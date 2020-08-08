%   function params = mlrSphericalNormGroup(params,<'justGetParams'>)
%
%   goal: Exports overlay data from multiple subjects located in the same study folder
%         into a common template space using Freesurfer's spherical normalization (keeping
%         only data located within the cortical sheet) and creates a new template MLR folder
%         (within the same folder) in which scans are concatenated overlays across all subjects.
%         Optionally concatenates subject data with their left-right-flipped counterpart. In this case
%         a twice left-right-flipped Freesurfer template surfer must be used (see surfRelaxFlipLR.m)
%         Optionally computes the group-average using mrLoadRet timeseries analysis.
%
%   usage:
%     params = mlrSphericalNormGroup(params,'justGetParams') %returns default parameters which need to be modified
%     mlrSphericalNormGroup(params) % run function with modified params
%
%   parameters:
%       params.studyDir:                      This is the folder in which subject-specific mrLoadRet folders are located (default: current folder)
%       params.mrLoadRetSubjectIDs            Names of the subjects (mrLoadRet folders) to include (default: all subfolders of current folder)
%       params.mrLoadRetSubjectLastView       Name of the last saved view in each subject's MLR folder (default: mrLastView.mat)
%       params.freesurferSubjectsFolder       Location of the Freesurfer subjects directory (default: from mrGetPref)
%       params.freesurferSubjectIDs           Freesurfer subject IDs corresponding to the mrLoadRet subject IDs
%       params.freesurferTemplateID           Freesurfer subject IF of the destination template (could be one of the subjects) (default: fsaverage)
%       params.mrLoadRetTemplateID            Name of the mrLoadRet directory where the normalized data will be located (deafult: 'Spherical Normalization')
%       params.mrLoadRetTemplateLastView      Name of the saved last view in the template MLR folder (default: mrLastView.mat)
%       params.subjectOverlayGroups           Group numbers of the overlays to normalize
%       params.subjectOverlayAnalyses         Analysis numbers of the overlays to normalize
%       params.subjectOverlays                Overlay numbers to normalize
%       params.templateOverlaysGroup          In what group of the group mrLoadRet folder the group-normalized overlays (concatenated in a scan) will be imported
%       params.combineLeftAndRight            If true, data will also be left-right flipped and combined across all hemispheres of all subjects (default: false)
%       params.lrFlippedFreesurferTemplateID  If combineLeftAndRight is true, name of the Freesurfer folder containing the left-right flipped surfaces of the left-right-flipped template
%       params.computeGroupAverage            If true, normalized overlays will be averaged across subjects (optionally hemispheres) (default: true)
%       params.templateOverlayNewNames        New names of the converted overlays
%       params.dryRun                         If true, just check whether the overlays and surfaces exist (default: false)
%
%   author: julien besle (28/07/2020)

function params = mlrSphericalNormGroup(params,varargin)

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end

if ieNotDefined('params')
  params = struct;
end

if fieldIsNotDefined(params,'studyDir')
  params.studyDir = pwd; % this is the folder in which subject-specific mrLoadRet folders are located
end
if fieldIsNotDefined(params,'mrLoadRetSubjectIDs')
  params.mrLoadRetSubjectIDs = dir(params.studyDir); % names of the subjects (mrLoadRet folders) to include
  params.mrLoadRetSubjectIDs = params.mrLoadRetSubjectIDs([params.mrLoadRetSubjectIDs(:).isdir] & ~ismember({params.mrLoadRetSubjectIDs(:).name},{'.','..','Spherical Normalization'}));
  params.mrLoadRetSubjectIDs = {params.mrLoadRetSubjectIDs(:).name};
end
if fieldIsNotDefined(params,'mrLoadRetSubjectLastView')
  params.mrLoadRetSubjectLastView = 'mrLastView.mat';
end
if fieldIsNotDefined(params,'freesurferSubjectsFolder')
  params.freesurferSubjectsFolder = mrGetPref('volumeDirectory'); %location of the Freesurfer subjects directory
end
if fieldIsNotDefined(params,'freesurferSubjectIDs')
  params.freesurferSubjectIDs = {''}; % Freesurfer subject IDs corresponding to the mrLoadRet subject IDs
end
if fieldIsNotDefined(params,'freesurferTemplateID')
  params.freesurferTemplateID = 'fsaverage'; % Freesurfer subject IF of the destination template (could be one of the subjects)
end
if fieldIsNotDefined(params,'mrLoadRetTemplateID')
  params.mrLoadRetTemplateID = 'Spherical Normalization'; % Name of the mrLoadRet directory where the normalized data will be located
end
if fieldIsNotDefined(params,'mrLoadRetTemplateLastView')
  params.mrLoadRetTemplateLastView = 'mrLastView.mat';
end
if fieldIsNotDefined(params,'subjectOverlayGroups')
  params.subjectOverlayGroups = 1; % group numbers of the overlays to normalize
end
if fieldIsNotDefined(params,'subjectOverlayAnalyses')
  params.subjectOverlayAnalyses = 1; % analysis numbers of the overlays to normalize
end
if fieldIsNotDefined(params,'subjectOverlays')
  params.subjectOverlays = 1; % overlay numbers to normalize
end
if fieldIsNotDefined(params,'templateOverlaysGroup')
  params.templateOverlaysGroup = params.freesurferTemplateID; % in what group of the group mrLoadRet folder the group-normalized overlays (concatenated in a scan) will be imported
end
if fieldIsNotDefined(params,'combineLeftAndRight')
  params.combineLeftAndRight = false; % if true, data will also be left-right flipped and combined across all hemispheres of all subjects
end
if fieldIsNotDefined(params,'lrFlippedFreesurferTemplateID')
  params.lrFlippedFreesurferTemplateID = ''; % combineLeftAndRight is true, name of the Freesurfer folder containing the left-right flipped surfaces of the left-right-flipped template
end
if fieldIsNotDefined(params,'computeGroupAverage')
  params.computeGroupAverage = true; % if true, normalized overlays will be averaged across subjects (optionally hemispheres)
end
if fieldIsNotDefined(params,'templateOverlayNewNames')
  params.templateOverlayNewNames = {}; % New names of the converted overlays
end
if fieldIsNotDefined(params,'dryRun')
  params.dryRun = false; % if true, just check whether the overlays and surfaces exist
end

if justGetParams, return; end


nSubjects = length(params.mrLoadRetSubjectIDs);
nOverlays = length(params.subjectOverlays);

if nSubjects ~= length(params.freesurferSubjectIDs)
  mrWarnDlg('(mlrSphericalNormGroup) There must be the same number of mrLoadRet and freesurfer subject IDs')
  return;
end
if nOverlays ~= length(params.subjectOverlayGroups)
  mrWarnDlg('(mlrSphericalNormGroup) Specify the same number of overlays and groups')
  return;
end
if nOverlays ~= length(params.subjectOverlayGroups)
  mrWarnDlg('(mlrSphericalNormGroup) Specify the same number of overlays and analyses')
  return;
end
if params.combineLeftAndRight && isempty(params.lrFlippedFreesurferTemplateID)
  mrWarnDlg('(mlrSphericalNormGroup) If combining left and right hemispheres, a left-right-flipped group Freesurfer ID must be specified')
  return;
end
if ismember(params.mrLoadRetTemplateID,params.mrLoadRetSubjectIDs)
  mrWarnDlg('(mlrSphericalNormGroup) The name of the mrLoadRet group folder cannot be the same as one of the subjects'' mrLoadRet ID')
  return;
end
if ~isempty(params.templateOverlayNewNames) &&  length(params.templateOverlayNewNames)~=nOverlays
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) There should be %d replacement overlay names',nOverlays));
  return;
end

badCharFixlist = {{'<','_lt_'},{'>','_gt_'},{':',';'},{'"',''''},{'/','_'},{'/','_'},{'|','_'},{'?','!'},{'*','%'}}; %characters that are not accepted in (Windows) file names

% create mrLoadRet group folder
mrLoadRetTemplateFolder = fullfile(params.studyDir,params.mrLoadRetTemplateID);
templateTseriesFolder = fullfile(mrLoadRetTemplateFolder,params.templateOverlaysGroup,'TSeries');
if ~params.dryRun
  if ~exist(mrLoadRetTemplateFolder,'dir')
    initMrLoadRetGroup = true;
    makeEmptyMLRDir(mrLoadRetTemplateFolder,'description=Spherical normalization folder','subject=group',...
                    'operator=Created by mlrSphericalNormGroup','defaultParams=1',sprintf('defaultGroup=%s', params.templateOverlaysGroup));
  else
    initMrLoadRetGroup = false;
    mkdir(mrLoadRetTemplateFolder,params.templateOverlaysGroup);
    mkdir(fullfile(mrLoadRetTemplateFolder,params.templateOverlaysGroup),'TSeries');
  end
end

mrQuit(0); % make sure there are no open MLR views

% export subject overlays to subject-specific scan spaces
for iSubj = 1:nSubjects
  fprintf('\n(mlrSphericalNormGroup) Exporting scan data for subject %s ...\n',params.mrLoadRetSubjectIDs{iSubj});
  % some OSs don't deal with files that have more than 259 characters (including path)
  maxNcharacters = 259 - (length(templateTseriesFolder)+1) - (length(params.mrLoadRetSubjectIDs{iSubj})+1) ...
                   - (max(params.combineLeftAndRight*(length(params.lrFlippedFreesurferTemplateID)),length(params.freesurferTemplateID))+1) - 7;
  
  cd(fullfile(params.studyDir,params.mrLoadRetSubjectIDs{iSubj}));
  thisView = mrLoadRet(params.mrLoadRetSubjectLastView,'No GUI');
  
  cOverlay = 0;
  exportNames{iSubj} = cell(0);
  convertNames{iSubj} = cell(0);
  convertFlippedNames{iSubj} = cell(0);
  for iOverlay = 1:nOverlays
    thisView = viewSet(thisView,'curGroup',params.subjectOverlayGroups(iOverlay));
    thisView = viewSet(thisView,'curAnalysis',params.subjectOverlayAnalyses(iOverlay));
    if params.subjectOverlays(iOverlay) <= viewGet(thisView,'nOverlays')
      thisView = viewSet(thisView,'curOverlay',params.subjectOverlays(iOverlay));
      fprintf('(mlrSphericalNormGroup) Will export Group %d, Analysis %d, overlay %d: %s\n',params.subjectOverlayGroups(iOverlay),params.subjectOverlayAnalyses(iOverlay),params.subjectOverlays(iOverlay),viewGet(thisView,'overlayName'));
      cOverlay = cOverlay + 1;
      overlayExists = true;
    else
      mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find overlay %d in group %d, analysis %d...',params.subjectOverlays(iOverlay),params.subjectOverlayGroups(iOverlay),params.subjectOverlayAnalyses(iOverlay)));
      if ~params.dryRun
        return;
      else
        overlayExists=false;
      end
    end
    
    if overlayExists
      groupName = viewGet(thisView,'groupName');
      if endsWith(groupName,'Volume')
        mrWarnDlg(sprintf('(mlrSphericalNormGrousprintf) Subject %s: group %s seems to be a volume version of a flat base and needs to be converted using flatVol2OriginalVolume',...
                          params.mrLoadRetSubjectIDs{iSubj},groupName));
      end
      overlayBaseName{iOverlay} = sprintf('%02d_%s',iOverlay,fixBadChars(viewGet(thisView,'overlayName'),badCharFixlist,[],maxNcharacters));
      if iSubj==1
        if isempty(params.templateOverlayNewNames)
          params.templateOverlayNewNames{iOverlay} = overlayBaseName{iOverlay};
        else
          params.templateOverlayNewNames{iOverlay} = sprintf('%02d_%s',iOverlay,fixBadChars(params.templateOverlayNewNames{iOverlay},badCharFixlist,[],maxNcharacters));
        end
      end
      % export overlays to NIFTI in scan space
      exportNames{iSubj}{cOverlay} = fullfile(templateTseriesFolder,sprintf('%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj}));
      convertNames{iSubj}{cOverlay} = fullfile(templateTseriesFolder,sprintf('%s_%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj},params.freesurferTemplateID));
      if params.combineLeftAndRight
        convertFlippedNames{iSubj}{cOverlay} = fullfile(templateTseriesFolder,sprintf('%s_%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj},params.lrFlippedFreesurferTemplateID));
      end
      if ~params.dryRun
        mrExport2SR(thisView,exportNames{iSubj}{cOverlay},0);
      end
    end
  end
  deleteView(thisView); % close view without saving (note that because there should be only one view in global variable MLR, this deletes MLR too
  
  % transform subject scans to group space
  fprintf('(mlrSphericalNormGroup) Converting data to %s template space ...\n',params.freesurferTemplateID);
  fsSphericalParams = freesurferSphericalNormalizationVolumes([],'justGetParams');
  fsSphericalParams.sourceVol = exportNames{iSubj};
  fsSphericalParams.fsSourceSubj = params.freesurferSubjectIDs{iSubj};
  fsSphericalParams.fsDestSubj = params.freesurferTemplateID;
  fsSphericalParams.destVol = convertNames{iSubj};
  fsSphericalParams.interpMethod = 'linear';
  fsSphericalParams.dryRun = params.dryRun;
  fsSphericalParamsOut = freesurferSphericalNormalizationVolumes(fsSphericalParams);
  fsSphericalParamsOut.destSurfRelaxVolume = fullfile(mrGetPref('volumeDirectory'),params.freesurferTemplateID,'surfRelax',[params.freesurferTemplateID '_mprage_pp.nii']); % DELETE
  if params.combineLeftAndRight
    fprintf('(mlrSphericalNormGroup) Converting data to %s template space ...\n',params.lrFlippedFreesurferTemplateID);
    fsSphericalParams.fsDestSubj = params.lrFlippedFreesurferTemplateID;
    fsSphericalParams.destVol = convertFlippedNames{iSubj};
    freesurferSphericalNormalizationVolumes(fsSphericalParams);
  end
end

if ~params.dryRun

  % delete subject-space exported NIFTI files
  for iSubj =1:nSubjects
    for iOverlay = 1:nOverlays
      delete(exportNames{iSubj}{iOverlay});
    end
  end

  % Concatenate normalized subject data
  cOverlay = 0;
  for iOverlay = 1:nOverlays
    cOverlay = cOverlay+1;
    templateOverlayName{cOverlay} = params.templateOverlayNewNames{iOverlay};
    fprintf('\n(mlrSphericalNormGroup) Concatenating subject data for overlay %s\n',templateOverlayName{cOverlay});
    for iSubj = 1:nSubjects
      [data(:,:,:,iSubj),hdr] = mlrImageLoad(convertNames{iSubj}{iOverlay});
    end
    hdr.dim = [hdr.dim nSubjects];
    hdr.filename = fullfile(templateTseriesFolder,[templateOverlayName{cOverlay} '.nii']);
    mlrImageSave(hdr.filename,data,hdr);
    if params.combineLeftAndRight
      cOverlay = cOverlay+1;
      templateOverlayName{cOverlay} = [params.templateOverlayNewNames{iOverlay} '_LRcombined'];
      fprintf('\n(mlrSphericalNormGroup) Concatenating hemisphere and subject data for overlay %s\n',templateOverlayName{cOverlay});
      for iSubj = 1:nSubjects
        [data(:,:,:,nSubjects+iSubj),hdr] = mlrImageLoad(convertFlippedNames{iSubj}{iOverlay});
      end
      hdr.dim = [hdr.dim nSubjects*2];
      hdr.filename = fullfile(templateTseriesFolder,[templateOverlayName{cOverlay} '.nii']);
      mlrImageSave(hdr.filename,data,hdr);
    end
    clear('data');
  end
  % delete temporary converted subject files
  for iSubj =1:nSubjects
    for iOverlay = 1:nOverlays
      delete(convertNames{iSubj}{iOverlay});
      if params.combineLeftAndRight
        delete(convertFlippedNames{iSubj}{iOverlay});
      end
    end
  end

  fprintf('(mlrSphericalNormGroup) Importing group data into MLR folder %s\n',params.mrLoadRetTemplateID);
  nTemplateOverlays = length(templateOverlayName);
  % initialize mrLoadRet view/import group
  cd(mrLoadRetTemplateFolder);
  if initMrLoadRetGroup
    subjectSessionParams = load(fullfile(params.studyDir,params.mrLoadRetSubjectIDs{1},'mrSession.mat'));
    [sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1');
    for iScan = 1:nTemplateOverlays
      groupParams.description{iScan} = templateOverlayName{iScan};
    end
    sessionParams.magnet = subjectSessionParams.session.magnet;
    sessionParams.coil = subjectSessionParams.session.coil;
    [sessionParams.pulseSequence,sessionParams.pulseSequenceText] = strtok(subjectSessionParams.session.protocol,':');
    mrInit(sessionParams,groupParams,'makeReadme=0','noPrompt=1');
    % load base anatomies
    % load whole-head MPRAGE anatomy
    thisView = mrLoadRet(params.mrLoadRetTemplateLastView,'No GUI');
    [fsPath,filename,ext] = fileparts(fsSphericalParamsOut.destSurfRelaxVolume);
    thisView = loadAnat(thisView,[filename ext],fsPath);
    thisView = viewSet(thisView,'basesliceindex',3); %set to axial view
    thisView = viewSet(thisView,'rotate',90);
    % import surfaces
    sides = {'left','right'};
    for iSide=1:2
      fprintf('(mlrSphericalNormGroup) Importing %s surface for %s\n',sides{iSide},params.freesurferTemplateID);
      importSurfParams.path = fsPath;
      importSurfParams.outerSurface = [params.freesurferTemplateID '_' sides{iSide} '_Inf.off']; % in order to view half-inflated surfaces, set the outer surface coordinates (outerCoords) to be the GM surface
      importSurfParams.outerCoords = [params.freesurferTemplateID '_' sides{iSide} '_GM.off']; % and the outer surface (outerSurface) to be the inflated surface
      importSurfParams.innerSurface = [params.freesurferTemplateID '_' sides{iSide} '_WM.off'];
      importSurfParams.innerCoords = [params.freesurferTemplateID '_' sides{iSide} '_WM.off'];
      importSurfParams.anatomy = [filename ext];
      importSurfParams.curv = [params.freesurferTemplateID '_' sides{iSide} '_Curv.vff'];
      base = importSurfaceOFF(importSurfParams);
      thisView = viewSet(thisView, 'newbase', base); %once the base has been imported into matlab, set it in the view
      thisView = viewSet(thisView,'corticalDepth',[0.2 0.8]); %set the range of of cortical depths used for displaying the overlays
    end

  else
    thisView = mrLoadRet(params.mrLoadRetTemplateLastView,'No GUI');
    thisView = viewSet(thisView,'newGroup',params.templateOverlaysGroup);
    [~,params] = importTSeries(thisView,[],'justGetParams=1','defaultParams=1');
    for iOverlay = 1:nOverlays
      params.pathname = fullfile(templateTseriesFolder,[overlayBaseName{iOverlay} '.nii']);
      [~,params] = importTSeries(thisView,params);
      if params.combineLeftAndRight
        params.pathname = fullfile(templateTseriesFolder,[overlayBaseName{iOverlay} '_LRcombined.nii']);
        [~,params] = importTSeries(thisView,params);
      end
    end
  end

  if params.computeGroupAverage
    statsParams.groupName = params.templateOverlaysGroup;
    statsParams.scanList = 1:nTemplateOverlays;
    thisView = timeSeriesStats(thisView,statsParams);
    % replace 0s by NaNs (because NaNs are written out in the Nifti files as 0; this is probably an MLR problem that can be fixed in mlrImageSave)
    [~,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
    params.combineFunction = '@(x)replaceByNaNs(x,0)';
    for iOverlay = 1:3 % only do the first 3 overlays (mean, median and std-deviation)
      params.outputName = ' ';
      params.overlayList = iOverlay;
      thisView = combineTransformOverlays(thisView,params);
    end
    thisView = viewSet(thisView,'deleteOverlay',1:6); % remove first 6 overlays
    thisView = viewSet(thisView,'curOverlay',1); % mean
    thisView = viewSet(thisView,'clipAcrossOverlays',false);

  end

  mrSaveView(thisView);
  deleteView(thisView);

end
