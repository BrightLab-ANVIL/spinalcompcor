function [Eigenvectors,numEigenvectors,Eigenvectors_orig,coeff,score,explained,fittedNoise,denoised]=SpinalCompCor(basedir,spinemask,fMRIpath,prefix,options)
% notSpine performs regression of non-spinal cord noise from fMRI data
% spine mask should include CSF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ****************************** IMPORTANT! ******************************
% NOTE: NEED TO OPEN MATLAB THROUGH TERMINAL FOR FSL COMMANDS TO WORK
% USE 'open /Applications/MATLAB_R2019b.app/' with proper version inserted
% startup.m file should also include FSL Setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% USAGE: 
% notSpine('data/folder','mask.nii.gz','data/folder/fmri.nii.gz','C01_BH')
% 
% ARGUMENTS:
% basedir           > path to directory with spine mask 'data/folder'
% spinemask         > name of spine mask file 'mask.nii.gz'
% fMRIpath          > full path to fMRI data to denoise 'data/folder/fmri.nii.gz'
% prefix            > desired prefix for output files
% 'plot',1          > show pop up plots/images (1 or 0) default=0
% 'task',2          > breath-hold or hand-grip (1 or 2, respectively)
% 'CO2path','path'  > need if task=1 if plotting
% 'Lpath', 'path'   > (HG) need if task=2 if plotting
% 'Rpath', 'path'   > (HG) need if task=2 if plotting
% 'numPC', 'max'    > choose # PCs (max(default), min, mode, optimized)
% 'regression', 1   > perform regression (1 or 0) default=0
%
% OUTPUTS:
% Eigenvectors      > the chosen and demeaned PCs for each slice (cell)
% numEigenvectors   > number of PCs chosen for each slice (vector)
% Eigenvectors_orig > non-demeaned version of Eigenvectors (cell)
% coeff             > refer to mathworks.com/help/stats/pca.html
% score             > refer to mathworks.com/help/stats/pca.html
% explained         > refer to mathworks.com/help/stats/pca.html
% fittedNoise       > noise fitted for fMRI data (only use if regression: 1)
% denoised          > denoised fMRI data (only use if regression: 1)
% vox_dilate        > num voxels to dilate mask (default: 18)
% 
%
% Written by Kimberly J Hemmerling
% Based on Barry et al. 2016 (Reproducibility of resting state spinal cord
% networks in healthy volunteers at 7 Tesla)
arguments
    basedir (1,:) char
    spinemask (1,:) char
    fMRIpath (1,:) char
    prefix (1,:) char
    options.plot (1,1) {mustBeMember(options.plot,[0,1])} = 0
    options.task (1,1) {mustBeMember(options.task,[0,1,2])} = 0
    options.CO2path (1,:) char = "-"
    options.Lpath (1,:) char = "-"
    options.Rpath (1,:) char = "-"
    options.numPC (1,:) {mustBeMember(options.numPC,["max","min","mode","optimized","five"])} = "max"
    options.regression (1,1) {mustBeMember(options.regression,[0,1])} = 0
    options.vox_dilate (1,1) {mustBeInteger(options.vox_dilate)} = 18
end

%% Do checks
tic
fprintf('\nBeginning...\n')
if options.plot==1
    if options.task==0
        warning('No task provided.')
    elseif options.task==1
        if options.CO2path== "-"; error('No C02 path provided.'); end
    elseif options.task==2
        if options.Lpath=="-"; error('No left HG path provided.'); end
        if options.Rpath=="-"; error('No left HG path provided.'); end
    end
end
%% Load nifti file for regression
fprintf('\nLoading fMRI nifti file\n')
fMRIdata=load_untouch_nii(fMRIpath);

% TESTIGN section   -   scaling and interp ****************************************************
% This section applies the scaling parameters to the fMRI data and then
% updates the header information to be consistent... (04/28/22) IFF the 
% slope is not zero to begin with !!
if fMRIdata.hdr.dime.scl_slope ~=0
    fMRIdata.img=fMRIdata.img * fMRIdata.hdr.dime.scl_slope + fMRIdata.hdr.dime.scl_inter;
    fMRIdata.hdr.dime.scl_slope=1;
    fMRIdata.hdr.dime.scl_inter=0;
    fMRIdata.hdr.dime.glmax=max(max(max(max(fMRIdata.img)))); % tbd as to whether I should change this
    fMRIdata.hdr.dime.glmin=min(min(min(min(fMRIdata.img))));
end
% disp(fMRIdata.hdr.dime)
if options.plot==1
    figure; imagesc(fMRIdata.img(:,:,10,15)); axis image; colormap gray
    title({'Volume 15, Slice 10','(check if voxel values are correct)'})
end
% *********************************************************************************************

fMRIall = squeeze(fMRIdata.img(:,:,:,:));
dim.x=size(fMRIall,1); dim.y=size(fMRIall,2); dim.z=size(fMRIall,3); dim.t=size(fMRIall,4);
% Initialize output image
fMRI_denoised=fMRIdata;
notspine=fMRIdata;

%% Create a not spine mask by inverting the spine mask (FSL)
% % (mask-1) * -1
% fprintf('\nCreate not-spine mask\n')
% temp=unix(['fslmaths ' basedir '/' spinemask ' -sub 1 -abs ' basedir '/' prefix '_nostpinemask.nii.gz']);
% if temp==1; error('ERROR! Issue running fslmaths. Make sure MATLAB was opened from terminal.'); end
% temp=unix(['fslmeants -i ' fMRIpath ' -m ' basedir '/' prefix '_nostpinemask.nii.gz --showall > ' basedir '/' prefix '_notSpine_timeseries.txt']);
% if temp==1; error('ERROR! Issue running fslmeants. Make sure MATLAB was opened from terminal.'); end
% clear temp
%% Create a not spine mask by dilating the spine mask and subtracting (FSL/SCT)
fprintf('\nCreate not-spine mask\n')
% Dilate spine mask
vox=num2str(options.vox_dilate);
temp=unix(['sct_maths -i ' basedir '/' spinemask ' -o ' basedir '/' prefix ...
    '_dilate.nii.gz -dilate ' vox ' -shape disk -dim 2']);
if temp~=0; error('ERROR! Issue running sct_maths. Make sure MATLAB was opened from terminal.'); end
% Subtract spine from dilated mask
temp=unix(['fslmaths ' basedir '/' prefix '_dilate.nii.gz -sub ' basedir '/' spinemask ...
    ' ' basedir '/' prefix '_dilate_notspinemask_temp.nii.gz']);
if temp~=0; error('ERROR! Issue running fslmaths. Make sure MATLAB was opened from terminal.'); end
% Set outer voxels to zero using fslmaths -roi
temp=unix(['fslmaths ' basedir  '/' prefix '_dilate_notspinemask_temp.nii.gz' ...
    ' -roi 3 ' num2str(dim.x-6) ' 3 ' num2str(dim.y-6) ' 0 ' num2str(dim.z) ' 0 1 ' ...
    basedir '/' prefix '_dilate_notspinemask.nii.gz']);
if temp~=0; error('ERROR! Issue running fslmaths. Make sure MATLAB was opened from terminal.'); end
if options.plot==1
    temp=unix(['fsleyes ' fMRIpath ' ' basedir '/' prefix '_dilate_notspinemask.nii.gz' ...
        ' -cm blue-lightblue -a 30 &']);
    if temp~=0; error('ERROR! Issue opening fsleyes. Make sure MATLAB was opened from terminal.'); end
end
temp=unix(['fslmeants -i ' fMRIpath ' -m ' basedir '/' prefix '_dilate_notspinemask.nii.gz --showall > ' ...
    basedir '/' prefix '_notSpine_timeseries.txt']);
if temp==1; error('ERROR! Issue running fslmeants. Make sure MATLAB was opened from terminal.'); end
clear temp

%% Load data
fprintf('\nLoading fMRI timeseries\n')
timeseries_ns=load([basedir '/' prefix '_notSpine_timeseries.txt']);
if options.plot==1 && options.task==1
    CO2=load(options.CO2path);
elseif options.plot==1 && options.task==2
    left=load(options.Lpath);
    right=load(options.Rpath);
end
% Need voxel references to do slicewise PCA (column 3)
voxelDir=timeseries_ns(1:3,:);
timeseries_ns(1:3,:)=[];  % Delete voxel references
% DON'T Transpose timeseries_ns (01/23/2022)

%% SLICEWISE PCA
% n Rows should be observations/samples (TIME POINTS)
% p Columns should be predictor variables (VOXELS)
% timeseries_ns: n x p ---- TIME POINTS x VOXELS
% [coeff,score,latent,tsquared,explained,mu]=pca(timeseries_ns);
% coeff: loadings
% score: components

% SLICEWISE PCA: 
fprintf('\nPerforming slicewise PCA\n')
sliceVoxelDir=cell(dim.z,1); % To reassemble variance/"eigencord" images)
coeff=cell(dim.z,1);    score=cell(dim.z,1);
latent=cell(dim.z,1);   explained=cell(dim.z,1);
for slice_num=0:dim.z-1
    idx=find(voxelDir(3,:)==slice_num);
    sliceVoxelDir{slice_num+1}=voxelDir(:,idx);
    slice=timeseries_ns(:,idx);
    [coeff{slice_num+1},score{slice_num+1},latent{slice_num+1},~,explained{slice_num+1},~]=...
        pca(slice,'Algorithm','svd','Economy',true);
end

%% Visualize PCs and task vectors from a middle slice
% Each column of score contains principal components (ordered by variance)
% This corresponds to each "explained" cell and index of that ^ column
slice=round(dim.z/2); 
if options.plot==1
    figure
    subplot(611); title("Slice " + slice)
    if options.task==1
        plot(CO2); xlim([1 length(CO2)])
    elseif options.task==2
        hold on
        plot(left); plot(right); xlim([1 length(left)])
        hold off
    end
    subplot(612)
    plot(score{slice}(:,1)); title("Principal Component 1, " + explained{12}(1) + "% variance")
    xlim([1 dim.t])
    subplot(613)
    plot(score{slice}(:,2)); title("Principal Component 2, " + explained{12}(2) + "% variance")
    xlim([1 dim.t])
    subplot(614)
    plot(score{slice}(:,3)); title("Principal Component 3, " + explained{12}(3) + "% variance")
    xlim([1 dim.t])
    subplot(615)
    plot(score{slice}(:,4)); title("Principal Component 4, " + explained{12}(4) + "% variance")
    xlim([1 dim.t])
    subplot(616)
    plot(score{slice}(:,5)); title("Principal Component 5, " + explained{12}(5) + "% variance")
    xlim([1 dim.t])
    saveas(gcf,[basedir '/' prefix '_PCs.jpg'])
end

%% Choose PCs
% up to 80% of the slice-wise cumulative variance or until the difference 
% between two successive eigenvalues was less than 5%
fprintf('\nSelect PCs\n')
numEigenvectors=zeros(dim.z,1);
for slice=1:dim.z
    i=1;
    var=110; % Set arbitrarily just so first comparison always goes
    cumulative_var=explained{slice}(i);
    while (cumulative_var<80) && (var-explained{slice}(i+1)>=5)
        i=i+1;
        var=explained{slice}(i);
        cumulative_var=cumulative_var+var;
    end
    % Set eigenvectors to regress out
    numEigenvectors(slice)=i;
end
disp(numEigenvectors')
opt_num=numEigenvectors; 

% Select number of Eigenvectors/PCs to actually choose based on user input
if options.numPC=="optimized"
    % Do nothing
elseif options.numPC=="max"
    numEigenvectors(:)=max(numEigenvectors);
elseif options.numPC=="min"
    numEigenvectors(:)=min(numEigenvectors);
elseif options.numPC=="mode"
    numEigenvectors(:)=mode(numEigenvectors);
elseif options.numPC=="five"
    numEigenvectors(:)=5;
end

% output optimized if not requested %
% Eigenvectors_opt is a cell with Eigenvector matricies for each slice
% it uses the "optimized" choice, if this ISN'T requested by user
if options.numPC~="optimized"
    Eigenvectors_opt=cell(dim.z,1); 
    for slice=1:dim.z
    num=opt_num(slice);
    Eigenvectors_opt{slice}=score{slice}(:,1:num); 
        for PC=1:size(Eigenvectors_opt{slice},2)
            % Demean Eigenvector regressors
            Eigenvectors_opt{slice}(:,PC)=Eigenvectors_opt{slice}(:,PC)-mean(Eigenvectors_opt{slice}(:,PC));
        end
    end
end
clear opt_num

% DEMEAN EIGENVECTORS %
% Eigenvectors_orig is a cell with Eigenvector matrices for each slice 
% of non-demeaned, requested numEigenvectors
Eigenvectors_orig=cell(dim.z,1); Eigenvectors=cell(dim.z,1);
for slice=1:dim.z
    num=numEigenvectors(slice);
    Eigenvectors_orig{slice}=score{slice}(:,1:num); 
    for PC=1:size(Eigenvectors_orig{slice},2)
        % Demean Eigenvector regressors
        Eigenvectors{slice}(:,PC)=Eigenvectors_orig{slice}(:,PC)-mean(Eigenvectors_orig{slice}(:,PC));
    end
end


clear slice
%% TESTING SECTION  -  visualize eigenvectors *************************************************
if options.plot==1
    for i=1:3
        slice=round(dim.z/4)*i;
        figure; numPlots=numEigenvectors(slice);
        for plt=1:numPlots
            subplot(numPlots,1,plt)
            plot(Eigenvectors{slice}(:,plt))
            xlim([0 dim.t]); 
            if plt==1; title("Slice " + slice + " Demeaned Eigenvectors"); end
        end
    end
end
clear i j numPlots slice
% *********************************************************************************************
%% Make origin of variance image
% bluewhitered = [0 0 0.5; 0 0.5 1; 1 1 1; 1 0 0; 0.5 0 0];
bluewhitered = [linspace(0,1,128)' linspace(0,1,128)' linspace(0.5,1,128)'; ...
    linspace(1,0.5,128)' linspace(1,0,128)' linspace(1,0,128)'];

% do a regression from orig. timeseries and pull beta weights for each component??

% OR loadings (which can be understood as the weights for each original 
% variable when calculating the principal component):
varIm=cell(5,1);
if isfolder([basedir '/figures'])==0
    mkdir([basedir '/figures']); end
for PC=1:5
    varIm{PC}=zeros(dim.x,dim.y,dim.z);
    for slice=1:dim.z
        for voxel=1:size(sliceVoxelDir{slice},2)
            x=sliceVoxelDir{slice}(1,voxel);
            y=sliceVoxelDir{slice}(2,voxel);
            varIm{PC}( x,y,slice )=coeff{slice}( voxel,PC );
        end
%         figure(15); imagesc( imrotate(varIm{PC}(:,:,slice),90) ); 
%         title("SLICE: " + slice + " PC: " + PC)
%         colormap gray; axis image; colorbar
%         pause
    end
    varIm_img=figure; s=1; sgtitle("COMPONENT: " + PC)
    for im=round(1: dim.z/9 :dim.z)
%         subplot(3,3,s); imagesc( imrotate(varIm{PC}(50:200,50:200,im),90) ); s=s+1;
        subplot(3,3,s); imagesc( imrotate(varIm{PC}(:,:,im),90) ); s=s+1;
        title("SLICE: " + im); caxis([-max(max(varIm{PC}(:,:,1))) max(max(varIm{PC}(:,:,1)))])
        colormap(bluewhitered); axis image; colorbar
    end
    
    PC_img=figure; s=1; sgtitle("COMPONENT: " + PC)
    for im=round(1: dim.z/9 :dim.z)
        subplot(3,3,s); 
        plot(score{im}(:,PC)); xlim([0 size(score{im}(:,PC),1)]); s=s+1;
    end
    % Save and close figures
    saveas(varIm_img,[basedir '/figures/' prefix '_varIm_' num2str(PC) '.svg'])
    saveas(PC_img,[basedir '/figures/' prefix '_PCIm_' num2str(PC) '.svg'])
    close(varIm_img); close(PC_img);
end

clear slice x y z s im
%% Perform regression in ALL voxels (do what 3dTfitter does)
if options.regression==1
    % Slices may have differing #s of regressors (if numPC: optimized)
    fprintf('\nRegression with %d PCs\n', numEigenvectors)
    Bfits=zeros(dim.x,dim.y,dim.z,max(numEigenvectors)+1);
    fittedNoise=zeros(dim.x,dim.y,dim.z,dim.t);
    counter=0;
    for slice=1:dim.z
        X=[ones(size(timeseries_ns,1),1) Eigenvectors{slice,1}]; % design matrix
        pinvX = pinv(X); clear X; % Calculate pseudoinverse
        for x=1:dim.x
            for y=1:dim.y
                if fMRIall(x,y,slice,1)>0
                    % Fit the data in each voxel
                    Y = double(squeeze(fMRIall(x,y,slice,:))); % 1 voxel timeseries
                    B = pinvX*Y; % Fit the GLM
                    Bfits(x,y,slice,1:length(B))=B; 
                    % B length: # eigenvectors + 1 (for mean)
                    for i=1:length(B(2:end))
                        % Voxelwise multiplicaiton of regressor * beta fit
                        temp=Eigenvectors{slice,1}(:,i) * B(i+1); % i+1 (don't remove mean)             
                        fittedNoise(x,y,slice,:)=squeeze(fittedNoise(x,y,slice,:))+temp;
                    end
                else
                    counter=counter+1;
                end
            end
        end
        clear pinvX B Y x y slice temp
    end
    fprintf('Number of skipped voxels in regression: %d \n', counter)
    fprintf('(One sagittal slice is %d voxels) \n', dim.y*dim.z)

    % Calculate denoised image and save files
    denoised=double(fMRIall)-fittedNoise;
    fMRI_denoised.img=denoised;
    notspine.img=fittedNoise;
    denoise_path=[basedir '/' prefix '_denoised' '.nii.gz'];
    save_untouch_nii(fMRI_denoised,denoise_path)
    save_untouch_nii(notspine,[basedir '/' prefix '_fittednoise.nii.gz'])
    save([basedir '/' prefix '_Eigenvector_otherData.mat'],'basedir','fMRIpath','numEigenvectors','Eigenvectors_orig','coeff','explained','latent','score')
    save([basedir '/' prefix '_Eigenvectors.mat'],'Eigenvectors')
    if exist('Eigenvectors_opt','var')==1; save([basedir '/' prefix '_Eigenvectors_opt.mat'],'Eigenvectors_opt'); end
    % placeholder to save variance im!!!!
else
    % Save just the Eigenvectors file to output use in regression
    save([basedir '/' prefix '_Eigenvector_otherData.mat'],'basedir','fMRIpath','numEigenvectors','Eigenvectors_orig','coeff','explained','latent','score')
    save([basedir '/' prefix '_Eigenvectors.mat'],'Eigenvectors')
    if exist('Eigenvectors_opt','var')==1; save([basedir '/' prefix '_Eigenvectors_opt.mat'],'Eigenvectors_opt'); end
    % placeholder to save variance im!!!!
end
%% Calculate tSNR - probably don't need this
% % meanIm_path=[basedir '/' prefix '_denoised_mean.nii.gz'];
% % temp=unix(['fslmaths ' denoise_path ' -Tmean ' meanIm_path]);
% % if temp==1; error('1ERROR! Issue calculating tSNR. Make sure MATLAB was opened from terminal.'); end
% % stdIM_path=[basedir '/' prefix '_denoised_std.nii.gz'];
% % temp=unix(['fslmaths ' denoise_path ' -Tstd ' stdIM_path]);
% % if temp==1; error('2ERROR! Issue calculating tSNR. Make sure MATLAB was opened from terminal.'); end
% % tSNR_path=[basedir '/' prefix '_denoised_tsnr.nii.gz'];
% % temp=unix(['fslmaths ' meanIm_path ' -div ' stdIM_path ' ' tSNR_path]);
% % if temp==1; error('3ERROR! Issue calculating tSNR. Make sure MATLAB was opened from terminal.'); end
% % temp=unix(['fslstats ' tSNR_path ' -k ' basedir '/' spinemask ' -M >> ' basedir '/' prefix '_tSNR.txt']);
% % if temp==1; error('4ERROR! Issue calculating tSNR. Make sure MATLAB was opened from terminal.'); end
% % clear temp meanIm_path stdIM_path tSNR_path
fprintf('\n...done!\n\n')
time=toc; 
if time>60 
    mins=time/60; fprintf('\nElapsed time is %f minutes.\n',mins); 
else
    fprintf('\nElapsed time is %f seconds.\n',time); 
end