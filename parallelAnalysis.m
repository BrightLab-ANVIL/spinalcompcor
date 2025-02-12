function [numPCs, not_converged, score, sd_score, latent, sd_latent] = parallelAnalysis(prefix, basedir, datafolder, resultsfolder, slices, spv, perm, maxiter, truncate)
%PARALLELANALYSIS compares IAAFT generated surrogates to fMRI timeseries
%data to determine a number of principal components in a region outside the
%spinal cord that explain more variance than the similar surrogate timeseries

%INPUTS
%prefix: prefix of notSpine timeseries file (str)
%basedir: path to directory with timeseries and desired output folder (str)
%datafolder: folder/path within basedir containing notspine timeseries (str)
%resultsfolder: folder/path within basedir for outputs (str)
%slices: slice indices to analyze (vector)
%spv: number of surrogates to generate per voxel (value)
%perm: number of surrogate slices to generate by permutation (value)
%maxiter: max number iterations for IAAFT (value)
%truncate: number of timepoints (TRs) to truncate from end of timeseries

%OUTPUTS
%numPCs: number of PCs decided for each slice (vector)
%not_converged: # surrogates for which IAAFT stopped before convergence (value) 
%score: Matlab's pca() output score for each slice in timeseries_ns (cell)
%sd_score: Matlab's pca() output score for each surrogate slice (cell)
%latent: Matlab's pca() output latent for each slice in timeseries_ns (cell)
%sd_latent: Matlab's pca() output latent for each surrogate slice (cell)

%load timeseries
fprintf('\nLoading timeseries data.')
timeseries_ns = load([basedir '/' datafolder '/' prefix '_notSpine_timeseries.txt']);
if truncate>0
    % truncate TRs if requested from the END
    % (e.g., useful for data with different number of timepoints)
    timeseries_ns=timeseries_ns(1:end-truncate,:);
    fprintf('Truncated %d timepoints from END. New size: %d TRs.\n',truncate,size(timeseries_ns,1)-3)
elseif truncate<0
    % truncate TRs if requested from the START (e.g., for steady state)
    truncate=-truncate;
    timeseries_ns=[timeseries_ns(1:3,:); timeseries_ns(4+truncate:end,:)];
    fprintf('Truncated %d timepoints from START. New size: %d TRs.\n',truncate,size(timeseries_ns,1)-3)
end
%predefine variables
not_converged = 0; num_slices = length(slices);
numPCs = zeros(num_slices, 1); [rows, cols] = size(timeseries_ns);
score = cell(num_slices, 1); sd_score = cell(num_slices, perm);
latent = zeros((rows-4), num_slices); sd_latent = cell(num_slices, 1);
sd_latent_mean = zeros((rows-4), num_slices); 

%generate surrogates for each voxel
for s = slices
    fprintf('\nPerforming analysis on slice %d. ', s);
    slice_ts = [];
    for c = 1:cols
        if timeseries_ns(3, c)==s
            slice_ts = [slice_ts, timeseries_ns(4:rows, c)];
        elseif timeseries_ns(3, c)>s && timeseries_ns(3, (c-1))==s
            break
        end
    end
    
    %create database of spv surrogates per voxel
    fprintf('Generating surrogates for voxel timeseries. ');
    [timepoints, voxels] = size(slice_ts); sd_archive = cell(voxels, 1);
    for v = 1:voxels
        [sd_archive{v, 1}, nc] = IAAFT(slice_ts(:, v), spv, maxiter);
        not_converged = nc + not_converged;
    end
    
    %generate permuted surrogate slices and perform pca
    fprintf('Creating surrogate slices and performing slicewise PCA. ')
    sd_slice = zeros(size(slice_ts));
    sd_l = zeros(timepoints-1, perm);
    for p = 1:perm
        for v = 1:voxels
            sd_slice(:, v) = sd_archive{v, 1}(:, randi(spv));
        end
        [~,sd_score{s+1, p},sd_l(:, p)] = pca(sd_slice);
    end
    
    %perform pca on template and compare eigenvalues to determine numPCs
    [~, score{s+1, 1}, l] = pca(slice_ts);
    numPCs(s+1, 1) = sum(l>=mean(sd_l, 2));
    latent(:, s+1) = l; 
    sd_latent{s+1, 1} = sd_l;
    sd_latent_mean(:, s+1) = mean(sd_l, 2);
    
    % Uncomment below to save scree plot
%     %make scree plot (save without opening)
%     f=figure('visible','off'); 
%     plot([mean(sd_l, 2), l]); xline(numPCs(s+1, 1));
%     xlabel("Number of PCs"); ylabel("Eigenvalues"); title([prefix ' Slice ' num2str(s) ' Scree Plot '], 'Interpreter', 'none');
%     legend("Average Surrogate Eigenvalues", "Template Eigenvalues", "PC Cutoff");
%     saveas(f,[resultsfolder '/' prefix '_scree_slice' num2str(s) '.svg']);
end
fprintf('\n');

%save results
save([resultsfolder '/' prefix '_inputs.mat'], 'maxiter', 'perm', 'spv');
save([resultsfolder '/' prefix '_results.mat'], 'latent', 'not_converged', 'numPCs', 'score', 'sd_latent_mean');
% currently not saving sd_latent or sd_score because they are large and not
% necessary 
end