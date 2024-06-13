% function mlrSegmentNormalizeSPM12(studyDir,subjectNames,T1regexp)
%
%   Runs SPM12 segmentation and normalization preprocessing on T1-weighted volumes. For each participant,
%   the T1-weighted volume is assumed to be located in the "Anatomy" folder of an mrLoadRet folder structure
%   [Not sure how standard this is, maybe the freesurfer subject folder should be used by default].
%   The T1-weighted volume is segmented and normalized to the MNI152 average brain. The output (including
%   the segmented T1, the forward and inverse deformation fields and the normalized T1) is placed in a new
%   subfolder "SPMnormalize" in each participant's folder.
%   
%   Inputs: - studyDir:     path of study folder containg participant's data. Each participant's data folder
%                           is assumed to have an mrLoadRet folder structure
%           - subjectNames: name of all participants' folder (cell array of strings)
%           - T1baseName:   regular expression uniquely identifying each participant's T1-weigthed NIFTI file name
%                           across all participants, each assumed to be located in the participant's Anatomy subfolder
%

function mlrSegmentNormalizeSPM12(studyDir,subjectNames,T1regexp)

if nargin ~= 3
  help("mlrSegmentNormalizeSPM12");
  return;
end

if ischar(subjectNames)
  subjectNames = {subjectNames};
end

% check that SPM is installed
spmInstallPath = which('spm');
if isempty(spmInstallPath)
  mrErrorDlg('(mlrSegmentNormalizeSPM12) SPM is not installed or has not been added to Matlab''s path');
else
  spmInstallPath = fileparts(spmInstallPath);
  SPMversion = regexp(spmInstallPath, '\w+$', 'match', 'once');
  if ~strcmp(SPMversion,'spm12')
    mrWarnDlg(sprintf('(mlrSegmentNormalizeSPM12) This function has only been tested with SPM 12. You are using %s',SPMversion));
  end
end

% get full path of study folder because SPM does not deal well with ~ expansion
cwd = pwd; % easiest way to do this is to cd and pwd
cd(studyDir);
studyDir = pwd;
cd(cwd); % go back to current folder

% this loop will perform preprocessing steps for all subjects specified in the subjectNames list
startTime = tic;
for subject = subjectNames

    fprintf('\nStarting segmentation/normalization for %s', subject{1}) % add a print statement to tell you which subject is being processed

    anatDir = fullfile(studyDir, subject{1}, 'Anatomy'); % this combines the root with a specific subject directory to create the full path to the folder containing anatomical data
    % find the T1-weighted WH volume and make sure it's unique
    anatFileName = dir(anatDir);
    anatFileName = {anatFileName(~cellfun(@isempty,regexp({anatFileName.name},T1regexp))).name};
%     anatFileName = spm_select('List', anatDir, T1regexp) % not using spm_select because it cannot deal with a ~ in the path (although I'm now forcing the path to be fully expanded)
    if numel(anatFileName) > 1
      mrWarnDlg('(mlrSegmentNormalizeSPM12) Non-unique T1-weigthed file in Anatomy folder. Skipping this participant');
      continue
    end
    
    % copy T1w into new folder
    normDir = fullfile(studyDir, subject{1}, 'SPMnormalize'); % this combines the root with a specific subject directory to create the full path to the folder containing anatomical data
    if ~exist(normDir,'dir')
      mkdir(normDir);
    end
    copyfile(fullfile(anatDir,anatFileName{1}),normDir);
    anatFilePath = fullfile(normDir,anatFileName{1});

    % the following was created using the SPM batch processing GUI, saving the batch as a script (_job.m file)
    % and modified to work in this loop and be installation independent
    matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(anatFilePath); % modified to use the current participant's T1-weighted as input
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmInstallPath '/tpm/TPM.nii,1']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmInstallPath '/tpm/TPM.nii,2']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmInstallPath '/tpm/TPM.nii,3']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmInstallPath '/tpm/TPM.nii,4']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmInstallPath '/tpm/TPM.nii,5']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmInstallPath '/tpm/TPM.nii,6']}; % modified to use the current machine's SPM installation
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                  NaN NaN NaN];
    matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(anatFilePath); % modified to use the current participant's T1-weighted as base name for output
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

%     save preprocessing_batch matlabbatch % save the setup into a matfile called preprocessing_batch.mat
    spm_jobman('run',matlabbatch) % execute the batch
    clear matlabbatch % clear matlabbatch

    % delete original T1w file
    delete(anatFilePath)

    disp(['Completed preprocessing for ', subject{1}]) % add a print statement telling you which subject has been processed
    

end

fprintf('\n')
toc(startTime)
