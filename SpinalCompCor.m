function [coeff,score,explained]=SpinalCompCor(basedir,spinemask,fMRIpath,prefix,options)
% SpinalCompCor performs PCA in spinal cord fMRI noise ROI and outputs
% components
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ****************************** IMPORTANT! ******************************
% NOTE: NEED TO OPEN MATLAB THROUGH TERMINAL FOR FSL & SCT COMMANDS TO WORK
% USE 'open /Applications/MATLAB_R2019b.app/' with proper version inserted
% startup.m file should also include FSL Setup
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% USAGE: 
% SpinalCompCor('data/folder','mask.nii.gz','data/folder/fmri.nii.gz','C01_BH')
% 
% ARGUMENTS:
% basedir           > path to directory with spine mask 'data/folder'
%                   > this is also the output directory
% spinemask         > name of spine mask file 'mask.nii.gz'
% fMRIpath          > full path to fMRI data 'data/folder/fmri.nii.gz'
% prefix            > desired prefix for output files
% 'plot',1          > show pop up spinemask image (1 or 0) default=0
% 'vox_dilate', 18  > number voxels to dilate spinemask (default: 18)
% 
% OUTPUTS:
% coeff             > refer to mathworks.com/help/stats/pca.html
% score             > refer to mathworks.com/help/stats/pca.html
% explained         > refer to mathworks.com/help/stats/pca.html
%
% OUTPUTS SAVED TO OUTPUT DIRECTORY:
% <prefix>_notSpine_timeseries.txt         > timeseries from noise ROI
% <prefix>_dilate_notspinemask.nii.gz      > the noise ROI mask
% <prefix>_dilate.nii.gz                   > dilated region
% <prefix>_dilate_notspinemask_temp.nii.gz > temp file
% <prefix>_Eigenvector_otherData.mat       > mat file with PCA and all other outputs
% 
% Code associated with:
% Hemmerling et al. Data-driven denoising in spinal cord fMRI with principal component analysis

arguments
    basedir (1,:) char
    spinemask (1,:) char
    fMRIpath (1,:) char
    prefix (1,:) char
    options.plot (1,1) {mustBeMember(options.plot,[0,1])} = 0
    options.numPC (1,:) {mustBeMember(options.numPC,["max","min","mode","optimized","five"])} = "max"
    options.vox_dilate (1,1) {mustBeInteger(options.vox_dilate)} = 18
end

%% Do checks
tic
fprintf('\nBeginning...\n')

if ~isfile(fullfile(basedir,spinemask))
    error(['Spine mask file not found: ' fullfile(basedir,spinemask)])
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
% if options.plot==1
%     figure; imagesc(fMRIdata.img(:,:,10,15)); axis image; colormap gray
%     title({'Volume 15, Slice 10','(check if voxel values are correct)'})
% end
% *********************************************************************************************

fMRIall = squeeze(fMRIdata.img(:,:,:,:));
dim.x=size(fMRIall,1); dim.y=size(fMRIall,2); dim.z=size(fMRIall,3); dim.t=size(fMRIall,4);
% Initialize output image
fMRI_denoised=fMRIdata;
notspine=fMRIdata;

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
% Need voxel references to do slicewise PCA (column 3)
voxelDir=timeseries_ns(1:3,:);
timeseries_ns(1:3,:)=[];  % Delete voxel references

%% SLICEWISE PCA
% n Rows should be observations/samples (TIME POINTS)
% p Columns should be predictor variables (VOXELS)
% timeseries_ns: n x p ---- TIME POINTS x VOXELS
% [coeff,score,latent,tsquared,explained,mu]=pca(timeseries_ns);
% coeff: loadings
% score: components

% SLICEWISE PCA: 
fprintf('\nPerforming slicewise PCA\n')
coeff=cell(dim.z,1);    score=cell(dim.z,1);
latent=cell(dim.z,1);   explained=cell(dim.z,1);
for slice_num=0:dim.z-1
    idx=find(voxelDir(3,:)==slice_num);
    slice=timeseries_ns(:,idx);
    [coeff{slice_num+1},score{slice_num+1},latent{slice_num+1},~,explained{slice_num+1},~]=...
        pca(slice,'Algorithm','svd','Economy',true);
end

%% Save just the Eigenvectors file to output use in regression
% Eigenvector_otherData
save([basedir '/' prefix '_Eigenvector_otherData.mat'],'basedir','fMRIpath',...
    'vox','coeff','explained','latent','score')

%%
fprintf('\n...done!\n\n')
time=toc; 
if time>60 
    mins=time/60; fprintf('\nElapsed time is %f minutes.\n',mins); 
else
    fprintf('\nElapsed time is %f seconds.\n',time); 
end