% function mlrImportSPMnormalization(studyDir,subjectNames,overwrite)
%
%   Imports MNI non-linear registration information into mrLoadRet session (mrSession.mat) so it can be used
%   in mrLoadRet (e.g. to display MNI coordinates). This function can be run for group of participants (mrLoadRet
%   sessions) within a given study folder. It assumes that the non-linear registration was computed using SPM12
%   and that the linear and non-linear transforms were saved as NIFTI images/headers in folder "SPMnormalize"
%   within each participant's mrLoadRet Folder.
%   
%   Inputs: - studyDir:     path of study folder containg participant's data. Each participant's data folder
%                           is assumed to have an mrLoadRet folder structure
%           - subjectNames: name of all participants' folder (cell array of strings)
%           - overwrite:    whether to overwrite any existing MNI info for each participant (default: false)
%

function mlrImportSPMnormalization(studyDir,subjectNames,overwrite)

if ~ismember(nargin,[0 2 3])
  help('mlrImportSPMnormalization');
  return
end

if ieNotDefined('overwrite')
  overwrite = false;
end

if ieNotDefined('studyDir') || ieNotDefined('subjectNames') % if the function is called without input, use the current folder
  [studyDir,subjectNames] = fileparts(pwd);
end

if ischar(subjectNames)
  subjectNames = {subjectNames};
end

for subject = subjectNames

  fprintf(sprintf('(mlrImportSPMnormalization) Importing MNI normalization info for %s...\n',subject{1}));
  
  cd(fullfile(studyDir,subject{1}));
  
  thisView = newView; % create a new view (no need to load the currently save view, and this is much faster)
  if isempty(thisView)
    mrWarnDlg('(mlrImportSPMnormalization) No mrSession.mat file found. Are you sure this folder is an mrLoadRet participant/session folder?');
    continue
  end

  mniInfo = viewGet(thisView,'mniInfo');
  if ~isempty(mniInfo)
     mrWarnDlg('(mlrImportSPMnormalization) MNI non-linear registration is already defined for this participant.');
     if overwrite
       mrWarnDlg('(mlrImportSPMnormalization) MNI info will be overwritten.');
     else
       mrWarnDlg('(mlrImportSPMnormalization) Set ''overwrite'' input variable to ''true'' to overwrite.');
       mrQuit(0,thisView);
       continue;
     end
  end
  
  if exist(fullfile(studyDir,subject{1},'SPMnormalize'),'dir')
    mniCoordMapFilename = dir(fullfile(studyDir,subject{1},'SPMnormalize'));
    mniCoordMapFilename = {mniCoordMapFilename(~cellfun(@isempty,regexp({mniCoordMapFilename.name},'^iy_'))).name}; % find all filenames starting with 'iy_', which is the default name for deformation fields from T1w to MNI output by SPM12
  end
  
  if ~exist(fullfile(studyDir,subject{1},'SPMnormalize'),'dir') || isempty(mniCoordMapFilename) % deformation fields from MNI to T1w
    mrWarnDlg(sprintf('(mlrImportSPMnormalization) There is no MNI normalization deformation map (SPMnormalize/y_*.nii) for this subject (%s). You first need to run mlrSPMnormalization',subject{1}));
    mrQuit(0,thisView);
    continue
  elseif numel(mniCoordMapFilename)>2
    keyboard; % there are more than one deformation maps. Think what to do
  end
  
  %read deformation coord map
  mniInfo = struct();
  [mniInfo.T1w2mniCoordMap,hdrToMNI] = mlrImageReadNifti(fullfile(studyDir,subject{1},'SPMnormalize',mniCoordMapFilename{1})); % T1w2mniCoordMap is the (non-linear) deformation coordinates map going from the T1w volume coordinates to MNI coordinates
  mniInfo.mag2T1w = inv(hdrToMNI.sform44 * shiftOriginXform); % mag2T1w is the (linear) transformation matrix from magnet coordinates to T1w volume coordinates (this should be applied before the non-linear transform to go from magnet to MNI coordinates)
  if exist(fullfile(studyDir,subject{1},'SPMnormalize',mniCoordMapFilename{1}(2:end)), 'file')
    [mniInfo.mnivol2magCoordMap,hdrFromMNI] = mlrImageReadNifti(fullfile(studyDir,subject{1},'SPMnormalize',mniCoordMapFilename{1}(2:end)));  % mnivol2magCoordMap is the inverse (non-linear) deformation coordinates map going from MNI volume coordinates to magnet coordinates
    mniInfo.mni2mnivol = inv(hdrFromMNI.sform44 * shiftOriginXform); % mni2mnivol is the (linear) transformation matrix from MNI coordinates to the volume coordinates of the inverse deformation coordinates map (this should be applied before the non-linear transform to go from MNI to magnet coordinates)
  end
  
  viewSet(thisView,'mniInfo',mniInfo);
  saveSession; % save mrSession.nat, because this is where we store mniInfo
  mrQuit(0,thisView); % quit mrLoadRet without saveing the view: this removes MLR from the global scope without modifying the view currently saved in mrLastView.mat
  disp(['(mlrImportSPMnormalization) Imported MNI registration info for ', subject{1}]) % add a print statement saying which subject has been processed
  
end