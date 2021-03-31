%   function [success,params] = mlrSphericalNormGroup(params,<'justGetParams'>)
%
%   goal: Exports overlay and ROI data from multiple subjects located in the same study folder
%         into a common template space using Freesurfer's spherical normalization (keeping
%         only data located within the cortical sheet) and creates a new template MLR folder
%         (within the same folder) in which scans are concatenated overlays and ROI masks across all subjects.
%         Optionally concatenates subject data with their left-right-flipped counterpart. In this case
%         a twice left-right-flipped Freesurfer template surfer must be used (see surfRelaxFlipLR.m)
%         Optionally computes group-average overlays and ROI probability maps.
%
%   usage:
%     params = mlrSphericalNormGroup(params,'justGetParams') %returns default parameters
%     ... % modify parameters
%     mlrSphericalNormGroup(params) % runs function with modified params
%
%   parameters:
%       params.studyDir                       This is the folder in which subject-specific mrLoadRet folders are located (default: current folder)
%       params.mrLoadRetSubjectIDs            Names of the subjects (mrLoadRet folders) to include (default: all subfolders of current folder)
%       params.mrLoadRetSubjectLastView       Name of the last saved view in each subject's MLR folder (default: mrLastView.mat)
%       params.freesurferSubjectsFolder       Location of the Freesurfer subjects directory (default: from mrGetPref)
%       params.freesurferSubjectIDs           Freesurfer subject IDs corresponding to the mrLoadRet subject IDs
%       params.fsSubjectSurfSuffix            Suffix(es) to add to the Freesurfer subject surface file names. Can be a single string or a cell array of strings (default: '')
%       params.freesurferTemplateID           Freesurfer subject ID of the destination template (could be one of the subjects) (default: fsaverage)
%       params.mrLoadRetTemplateID            Name of the mrLoadRet directory where the normalized data will be located (default: 'Spherical Normalization')
%       params.mrLoadRetTemplateLastView      Name of the saved last view in the template MLR folder (default: mrLastView.mat)
%       params.subjectOverlayGroups           Group names or numbers of the overlays to normalize
%       params.subjectOverlayAnalyses         Analysis names numbers of the overlays to normalize
%       params.subjectOverlays                Names or numbers or the overlay to normalize
%       params.subjectROIs                    Names or numbers or the ROIs to normalize
%       params.templateOverlayGroupNums       In what new group(s) of the template mrLoadRet folder the group-normalized overlays (concatenated in a scan) will be imported.
%                                             (must be an index into params.templateOverlayGroupNames)
%       params.templateOverlayGroupNames      Name(s) of the template group(s) (must not already exist)
%       params.combineLeftAndRight            If true, data will also be left-right flipped and combined across all hemispheres of all subjects (default: false)
%       params.lrFlippedFreesurferTemplateID  If combineLeftAndRight is true, name of the Freesurfer folder containing the left-right flipped surfaces of the left-right-flipped template
%       params.lrFlipTransformFunctions       Voxelwise transformation to apply to the overlay data when left-right flipping, useful for orientation sensitive measures
%       params.lrFlipTransform                Which overlays to transform. 0 for no transformation (default)
%       params.computeGroupAverage            If true, normalized overlays will be averaged across subjects (optionally hemispheres) (default: true)
%       params.templateOverlayNewNames        New names of the converted overlays
%       params.computeROIprobabilityMaps      If true, normalized overlays will be averaged across subjects (optionally hemispheres) (default: true)
%       params.templateROInewNames            New names of the converted ROIs
%       params.cropScans                      Whether to crop concatenated scans to smallest volume including all subject's data (default = true)
%       params.dryRun                         If true, just check whether the overlays and surfaces exist (default: false)
%
%   author: julien besle (28/07/2020)

function [success,params] = mlrSphericalNormGroup(params,varargin)

success = true;
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
if fieldIsNotDefined(params,'fsSubjectSurfSuffix')
  params.fsSubjectSurfSuffix = {''}; % Suffix(es) to add to the Freesurfer subject surfaces
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
if fieldIsNotDefined(params,'subjectOverlays')
  params.subjectOverlays = []; % overlay names or numbers to normalize
end
if fieldIsNotDefined(params,'subjectOverlayGroups')
  params.subjectOverlayGroups = ones(size(params.subjectOverlays)); % group numbers of the overlays to normalize
end
if fieldIsNotDefined(params,'subjectOverlayScans')
  params.subjectOverlayScans = ones(size(params.subjectOverlays)); % scan numbers of the overlays to normalize (within a given subjects, all scans should have identical xform)
end
if fieldIsNotDefined(params,'subjectOverlayAnalyses')
  params.subjectOverlayAnalyses = ones(size(params.subjectOverlays)); % analysis numbers of the overlays to normalize
end
if fieldIsNotDefined(params,'subjectROIs')
  params.subjectROIs = []; % ROI names or numbers to normalize
end
if fieldIsNotDefined(params,'subjectROIsides')
  params.subjectROIsides = []; % What hemisphere the ROI is in (1 = left, 2 = right). When averaging across left and right hemispheres, left and right ROIs should be matched and given in the same order.
end
if fieldIsNotDefined(params,'subjectROIbase')
  params.subjectROIbase = []; % export space of the ROI (by default, same as the overlay scan space)
end
if fieldIsNotDefined(params,'templateOverlayGroupNums')
  params.templateOverlayGroupNums = []; % in what group of the group mrLoadRet folder the group-normalized overlays (concatenated in a scan) will be imported
end
if fieldIsNotDefined(params,'templateOverlayGroupNames')
  params.templateOverlayGroupNames = params.freesurferTemplateID; % Names of the group of the template mrLoadRet folder the group-normalized overlays (concatenated in a scan) will be imported
end
if fieldIsNotDefined(params,'templateROIgroupName')
  params.templateROIgroupName = 'ROI group'; % Name of the group of the template mrLoadRet folder the group-normalized ROI masks (concatenated in a scan) will be imported
end
if fieldIsNotDefined(params,'combineLeftAndRight')
  params.combineLeftAndRight = false; % if true, data will also be left-right flipped and combined across all hemispheres of all subjects
end
if fieldIsNotDefined(params,'lrFlippedFreesurferTemplateID')
  params.lrFlippedFreesurferTemplateID = ''; % if combineLeftAndRight is true, name of the Freesurfer folder containing the left-right flipped surfaces of the left-right-flipped template
end
if fieldIsNotDefined(params,'lrFlipTransformFunctions')
  params.lrFlipTransformFunctions = {};
end
if fieldIsNotDefined(params,'lrFlipTransform')
  params.lrFlipTransform = 0;
end
if fieldIsNotDefined(params,'computeGroupAverage')
  params.computeGroupAverage = true; % if true, normalized overlays will be averaged across subjects (optionally hemispheres)
end
if fieldIsNotDefined(params,'templateOverlayNewNames')
  params.templateOverlayNewNames = {}; % New names of the converted overlays
end
if fieldIsNotDefined(params,'computeROIprobabilityMaps')
  params.computeROIprobabilityMaps = true; % if true, normalized ROI masks will be averaged across subjects (optionally hemispheres)
end
if fieldIsNotDefined(params,'templateROInewNames')
  params.templateROInewNames = {}; % New names of the converted ROIs. When combining across left and right hemispheres, the number of new ROI names should match the number of left/right ROI pairs and be given in the same order
end
if fieldIsNotDefined(params,'cropScans')
  params.cropScans = true;
end
if fieldIsNotDefined(params,'dryRun')
  params.dryRun = false; % if true, just check whether the overlays and surfaces exist
end

if justGetParams
  return;
else
  success = false;
end


if ischar(params.subjectOverlays)
  params.subjectOverlays = {params.subjectOverlays};
end
if ischar(params.subjectROIs)
  params.subjectROIs = {params.subjectROIs};
end
if ischar(params.subjectROIbase)
  params.subjectROIbase = {params.subjectROIbase};
end
if ischar(params.fsSubjectSurfSuffix)
  params.fsSubjectSurfSuffix = {params.fsSubjectSurfSuffix};
end


nSubjects = length(params.mrLoadRetSubjectIDs);
nOverlays = length(params.subjectOverlays);
nROIs = length(params.subjectROIs);
nTemplateGroups =  max(params.templateOverlayGroupNums);
if isempty(nTemplateGroups)
  nTemplateGroups = 0;
end

if numel(params.subjectOverlayGroups)==1
  params.subjectOverlayGroups = repmat(params.subjectOverlayGroups,1,nOverlays);
end
if numel(params.subjectOverlayScans)==1
  params.subjectOverlayScans = repmat(params.subjectOverlayScans,1,nOverlays);
end
if numel(params.subjectOverlayAnalyses)==1
  params.subjectOverlayAnalyses = repmat(params.subjectOverlayAnalyses,1,nOverlays);
end
if numel(params.templateOverlayGroupNums)==1
  params.templateOverlayGroupNums = repmat(params.templateOverlayGroupNums,1,nOverlays);
end
if numel(params.subjectROIbase)==1
  params.subjectROIbase = repmat(params.subjectROIbase,1,nSubjects);
end
if numel(params.fsSubjectSurfSuffix)==1
  params.fsSubjectSurfSuffix = repmat(params.fsSubjectSurfSuffix,1,nSubjects);
end
if numel(params.lrFlipTransform)==1
  params.lrFlipTransform = repmat(params.lrFlipTransform,1,nOverlays);
end
if numel(params.computeGroupAverage)==1
  params.computeGroupAverage = repmat(params.computeGroupAverage,1,nTemplateGroups);
end



if nSubjects ~= length(params.freesurferSubjectIDs)
  mrWarnDlg('(mlrSphericalNormGroup) There must be the same number of mrLoadRet and freesurfer subject IDs')
  return;
end
if nSubjects ~= length(params.fsSubjectSurfSuffix)
  mrWarnDlg('(mlrSphericalNormGroup) The subject surface suffix must be a single string or a cells of strings of length equal to the number of subjects')
  return;
end
if nOverlays ~= length(params.subjectOverlayGroups)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of subject groups must be a single value or a vector of same length as the number of overlays (%d)',nOverlays))
  return;
end
if nOverlays ~= length(params.subjectOverlayScans)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of subject scans must be a single value or a vector of same length as the number of overlays (%d)',nOverlays))
  return;
end
if nOverlays ~= length(params.subjectOverlayAnalyses)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of subject analyses must be a single value or a vector of same length as the number of overlays (%d)',nOverlays))
  return;
end
if nOverlays ~= length(params.templateOverlayGroupNums)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The template group number must be a scalar or a vector of same length as the number of subject overlays (%d)',nOverlays))
  return;
end
if nROIs ~= length(params.subjectROIsides)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of ROI sides must match the number of ROIs (%d)',nROIs))
  return;
end
if params.combineLeftAndRight
  if nnz(params.subjectROIsides==1)~=nnz(params.subjectROIsides==2)
    mrWarnDlg('(mlrSphericalNormGroup) When combining across left and right hemispheres, left and right ROIs should be matched (and be given in the same order)')
    return;
  end
  if length(params.templateROInewNames) ~= nROIs/2
    mrWarnDlg(sprintf('(mlrSphericalNormGroup) When combining ROIs across left and right hemispheres, provide new ROI names that do not include ''left'' or ''right''. Their number and order should match that of left/right ROI pairs (%d)',nROIs/2));
    return;
  end
else
  if ~isempty(params.templateROInewNames) && length(params.templateROInewNames) ~= nROIs
    mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of new ROI names should match the number of ROIs (%d)',nROIs));
    return;
  end
end
if isempty(params.subjectOverlays) && isempty(params.subjectROIbase)
  mrWarnDlg('(mlrSphericalNormGroup) If exporting only ROIs and no overlay, you must specify an ROI base space.')
  return;
end
if ~isempty(params.subjectROIbase) && nSubjects ~= length(params.subjectROIbase)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The number of ROI bases must be a scalar or a vector of same length as the number of subjects (%d)',nROIs))
  return;
end
if nOverlays ~= length(params.lrFlipTransform)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) lrFlipTransform must be a scalar or a vector of same length as the number of subject overlays and ROIs (%d)',nOverlays))
  return;
end
if ~isempty(params.templateOverlayGroupNums) && min(params.templateOverlayGroupNums)~=1 || (length(unique(params.templateOverlayGroupNums))>1 && any(diff(unique(params.templateOverlayGroupNums)))~=1)
  mrWarnDlg('(mlrSphericalNormGroup) Template group numbers must be consecutive and start at 1')
  return;
end
if nTemplateGroups > length(params.templateOverlayGroupNames)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The largest template group number (%d) cannot be larger than the number of template group names',nTemplateGroups,length(params.templateOverlayGroupNames)))
  return;
end
if strcmp(params.templateROIgroupName,'ROIs')
  mrWarnDlg('(mlrSphericalNormGroup) The template ROI group cannot be named ''ROIs''');
  return;
end
if nTemplateGroups > length(params.computeGroupAverage)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) computeGroupAverage must be a scalar or have the same number of elements as the number of template groups (%d)',nTemplateGroups))
  return;
end
if params.combineLeftAndRight && isempty(params.lrFlippedFreesurferTemplateID)
  mrWarnDlg('(mlrSphericalNormGroup) If combining left and right hemispheres, a left-right-flipped group Freesurfer ID must be specified')
  return;
end
if max(params.lrFlipTransform)>length(params.lrFlipTransformFunctions)
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) The largest LR-flip function  number (%d) cannot be larger than the number LR-flip transform functions (%d)',max(params.lrFlipTransform),length(params.lrFlipTransformFunctions)))
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
sides = {'left','right'};
if params.combineLeftAndRight
  whichROI = zeros(1,nROIs);
  whichROI(params.subjectROIsides==1) = 1:nROIs/2;
  whichROI(params.subjectROIsides==2) = 1:nROIs/2;
end

if params.dryRun
  fprintf('\n(mlrSphericalNormGroup) THIS IS A DRY RUN: no data will be exported or written to disc\n');
end

% create mrLoadRet group folder
mrLoadRetTemplateFolder = fullfile(params.studyDir,params.mrLoadRetTemplateID);
temporaryTseriesFolder = fullfile(params.studyDir,'tempTSeries');
if ~params.dryRun
  mkdir(temporaryTseriesFolder);
end
for iGroup = 1: nTemplateGroups + (nROIs>0)
  if iGroup >  nTemplateGroups
    templateGroupFolder{iGroup} =  params.templateROIgroupName;
  else
    templateGroupFolder{iGroup} = params.templateOverlayGroupNames{iGroup};
  end
  templateTseriesFolder{iGroup} = fullfile(mrLoadRetTemplateFolder,templateGroupFolder{iGroup},'TSeries');
end
if ~params.dryRun
  if ~exist(mrLoadRetTemplateFolder,'dir')
    initMrLoadRetGroup = true;
    makeEmptyMLRDir(mrLoadRetTemplateFolder,'description=Spherical normalization folder','subject=group',...
                    'operator=Created by mlrSphericalNormGroup','defaultParams=1',sprintf('defaultGroup=%s', templateGroupFolder{1}));
    for iGroup = 2: nTemplateGroups + (nROIs>0)
      mkdir(mrLoadRetTemplateFolder,templateGroupFolder{iGroup});
      mkdir(templateTseriesFolder{iGroup});
    end
  else
    initMrLoadRetGroup = false;
    cd(mrLoadRetTemplateFolder)
    thisView = mrLoadRet([],'No GUI');
    groupExists = false;
    for iGroup = 1: nTemplateGroups + (nROIs>0)
      if isempty(viewGet(thisView,'groupnum',templateGroupFolder{iGroup}))
        mkdir(mrLoadRetTemplateFolder,templateGroupFolder{iGroup});
        mkdir(templateTseriesFolder{iGroup});
      else
        mrWarnDlg(sprintf('(mlrSphericalNormGroup) Group ''%s'' already exists in template folder ''%s''',templateGroupFolder{iGroup},params.mrLoadRetTemplateID));
        groupExists = true;
      end
    end
    deleteView(thisView);
    if groupExists
      return;
    end
  end
end

mrQuit(0); % make sure there are no open MLR views

% export subject overlays to subject-specific scan spaces
cNiftiConcat = zeros(nTemplateGroups+(nROIs>0),3);
multipleScans = length(unique(params.subjectOverlayScans))>1;
for iSubj = 1:nSubjects
  fprintf('\n(mlrSphericalNormGroup) Exporting scan data for subject %s ...\n',params.mrLoadRetSubjectIDs{iSubj});
  % some OSs don't deal with files that have more than 259 characters (including path)
  maxNcharacters = 259 - (length(temporaryTseriesFolder)+1) - (length(params.mrLoadRetSubjectIDs{iSubj})+1) ...
                   - (max(params.combineLeftAndRight*(length(params.lrFlippedFreesurferTemplateID)),length(params.freesurferTemplateID))+1) ...
                   - multipleScans*6 - 4;
  cd(fullfile(params.studyDir,params.mrLoadRetSubjectIDs{iSubj}));
  thisView = mrLoadRet(params.mrLoadRetSubjectLastView,'No GUI');
  
  cOverlay = 0;
  exportNames = cell(0);
  convertNames = cell(0,1+params.combineLeftAndRight);
  overlayExists = false(1,nOverlays);
  roiExists = false(1,nROIs);
  roiExportList = [];
  for iOverlay = 1:nOverlays
    if iscell(params.subjectOverlayGroups)
      if ischar(params.subjectOverlayGroups{iOverlay})
        groupNum = viewGet(thisView,'groupNum',params.subjectOverlayGroups{iOverlay});
        if isempty(groupNum)
          mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find group ''%s''',params.subjectOverlayGroups{iOverlay}));
        end
      else
        groupNum = params.subjectOverlayGroups{iOverlay};
      end
    else
      groupNum = params.subjectOverlayGroups(iOverlay);
    end
    if ~isempty(groupNum)
      if groupNum > viewGet(thisView,'nGroups')
        mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find group %d',groupNum));
      else
        thisView = viewSet(thisView,'curGroup',groupNum);
        groupName = viewGet(thisView,'groupName');
        
        scanNum = viewGet(thisView,'nscans');
        if params.subjectOverlayScans(iOverlay) > scanNum
          mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find scan %d in group %d',params.subjectOverlayScans(iOverlay),groupNum));
        else
          thisView = viewSet(thisView,'curScan',scanNum);
          if iOverlay==1
            scanHdr = viewGet(thisView,'niftiHdr');
            thisScanHdr = scanHdr;
          else
            thisScanHdr = viewGet(thisView,'niftiHdr');
          end
          if any(any(abs(scanHdr.sform44 - thisScanHdr.sform44)>1e-4)) || ...
             any(abs(scanHdr.dim(2:4) - thisScanHdr.dim(2:4))>1e-4) || ...
             any(abs(scanHdr.pixdim(2:4) - thisScanHdr.pixdim(2:4))>1e-4)
            mrWarnDlg(sprintf('(mlrSphericalNormGroup) Header for scan %d in group %d differs from the rest',params.subjectOverlayScans(iOverlay),groupNum));
          else
            if iscell(params.subjectOverlayAnalyses)
              if ischar(params.subjectOverlayAnalyses{iOverlay})
                analysisNum = viewGet(thisView,'AnalysisNum',params.subjectOverlayAnalyses{iOverlay});
                if isempty(analysisNum)
                  mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find analysis ''%s'' in group ''%s''',params.subjectOverlayAnalyses{iOverlay},groupName));
                end
              else
                analysisNum = params.subjectOverlayAnalyses{iOverlay};
              end
            else
              analysisNum = params.subjectOverlayAnalyses(iOverlay);
            end
            if ~isempty(analysisNum)
              if analysisNum > viewGet(thisView,'nAnalyses')
                mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find analysis %d in group ''%s''',analysisNum,groupName));
              else
                thisView = viewSet(thisView,'curAnalysis',analysisNum);
                analysisName = viewGet(thisView,'analysisName');
                if iscell(params.subjectOverlays)
                  if ischar(params.subjectOverlays{iOverlay})
                    overlayNum = viewGet(thisView,'overlayNum',params.subjectOverlays{iOverlay});
                    if isempty(overlayNum)
                      mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find overlay ''%s'' in group ''%s'', analysis ''%s''',params.subjectOverlays{iOverlay},groupName,analysisName));
                    end
                  else
                    overlayNum = params.subjectOverlays(iOverlay);
                  end
                else
                  overlayNum = params.subjectOverlays(iOverlay);
                end
                if ~isempty(overlayNum)
                  if overlayNum > viewGet(thisView,'nOverlays')
                    mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find overlay %d in group ''%s'', analysis ''%s''...',overlayNum,groupName,analysisName));
                  else
                    thisView = viewSet(thisView,'curOverlay',overlayNum);
                    overlayName = viewGet(thisView,'overlayName');
                    overlayExists(iOverlay)=true;
                    fprintf('(mlrSphericalNormGroup) Will export Group %d, Analysis %d, overlay %d: ''%s''\n',groupNum,analysisNum,overlayNum,overlayName);
                    cOverlay = cOverlay + 1;
                  end
                end
              end
            end
          end
        end

        if overlayExists(iOverlay)
          if length(groupName) > 6 && isequal(groupName(end-5:end),'Volume')
            mrWarnDlg(sprintf('(mlrSphericalNormGrousprintf) Subject %s: group %s seems to be a volume version of a flat base and needs to be converted using flatVol2OriginalVolume',...
                              params.mrLoadRetSubjectIDs{iSubj},groupName));
          end
          overlayBaseName{iOverlay} = fixBadChars(viewGet(thisView,'overlayName'),badCharFixlist,[],maxNcharacters);
          if multipleScans
            overlayBaseName{iOverlay} = sprintf('%s_Scan%d',overlayBaseName{iOverlay},scanNum);
          end
          if iSubj==1
            if isempty(params.templateOverlayNewNames)
              templateOverlayNewNames{cOverlay} = overlayBaseName{iOverlay};
            else
              templateOverlayNewNames{cOverlay} = fixBadChars(params.templateOverlayNewNames{iOverlay},badCharFixlist,[],maxNcharacters);
            end
          end
          % export overlays to NIFTI in scan space
          exportNames{cOverlay} = fullfile(temporaryTseriesFolder,sprintf('%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj}));
          convertNames{cOverlay,1} = fullfile(temporaryTseriesFolder,sprintf('%s_%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj},params.freesurferTemplateID));
          if params.combineLeftAndRight
            convertNames{cOverlay,2} = fullfile(temporaryTseriesFolder,sprintf('%s_%s_%s.nii',overlayBaseName{iOverlay},params.mrLoadRetSubjectIDs{iSubj},params.lrFlippedFreesurferTemplateID));
          end
          if ~params.dryRun
            mrExport2SR(thisView,exportNames{cOverlay},0);
          end
        end
      end
    end
  end
  
  % export ROIs to masks in subject's anatomical space
  if ~isempty(params.subjectROIbase{iSubj}); % first make sure we know what space to export it to
    roiBaseNum = viewGet(thisView,'baseNum',params.subjectROIbase{iSubj});
    if isempty(roiBaseNum)
      mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find base %s...',params.subjectROIbase{iSubj}));
    end
  else
    roiBaseNum = [];
  end
  cROI = 0;
  for iRoi = 1:nROIs
    if iscell(params.subjectROIs)
      if ischar(params.subjectROIs{iRoi})
        roiNum = viewGet(thisView,'roiNum',params.subjectROIs{iRoi});
        if isempty(roiNum)
          mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find ROI ''%s''',params.subjectROIs{iRoi}));
        end
      else
        roiNum = params.subjectROIs(iOverlay);
      end
    else
      roiNum = params.subjectROIs(iOverlay);
    end
    if ~isempty(roiNum)
      if roiNum > viewGet(thisView,'nROIs')
        mrWarnDlg(sprintf('(mlrSphericalNormGroup) Cannot find ROI %d...',overlayNum));
      else
        thisView = viewSet(thisView,'curROI',roiNum);
        roiName = viewGet(thisView,'roiName');
        roiExists(iRoi)=true;
        fprintf('(mlrSphericalNormGroup) Will export ROI %d: ''%s''\n',roiNum,roiName);
        cROI = cROI + 1;
        roiExportList(cROI) = roiNum;
      end
    end
    if roiExists(iRoi)
      roiBaseName{iRoi} = fixBadChars(roiName,badCharFixlist,[],maxNcharacters);
      if iSubj==1
        if isempty(params.templateROInewNames)
          templateOverlayNewNames{cOverlay+cROI} = roiBaseName{iRoi};
        else
          if params.combineLeftAndRight
            templateOverlayNewNames{cOverlay+cROI} = params.templateROInewNames{whichROI(iRoi)};
          else
            templateOverlayNewNames{cOverlay+cROI} = params.templateROInewNames{iRoi};
          end
          templateOverlayNewNames{cOverlay+cROI} = fixBadChars(templateOverlayNewNames{cOverlay+cROI},badCharFixlist,[],maxNcharacters);
        end
      end
      % export overlays to NIFTI in scan space
      exportNames{cOverlay+cROI} = fullfile(temporaryTseriesFolder,sprintf('%s_%s.nii',roiBaseName{iRoi},params.mrLoadRetSubjectIDs{iSubj}));
      convertNames{cOverlay+cROI,1} = fullfile(temporaryTseriesFolder,sprintf('%s_%s_%s.nii',roiBaseName{iRoi},params.mrLoadRetSubjectIDs{iSubj},params.freesurferTemplateID));
      if params.combineLeftAndRight
        convertNames{cOverlay+cROI,2} = fullfile(temporaryTseriesFolder,sprintf('%s_%s_%s.nii',roiBaseName{iRoi},params.mrLoadRetSubjectIDs{iSubj},params.lrFlippedFreesurferTemplateID));
      end
      if ~params.dryRun 
        if isempty(roiBaseNum)
          mlrExportROI(thisView,exportNames{cOverlay+cROI},'scanNum',scanNum,'grouNnum',groupNum); % export to scan space
        else
          mlrExportROI(thisView,exportNames{cOverlay+cROI},'baseNum',roiBaseNum);
        end
      end
    end
    
  end
  
  deleteView(thisView); % close view without saving (note that because there should be only one view in global variable MLR), this deletes MLR too
  
  if all(overlayExists) && all(roiExists) || params.dryRun
    % transform subject scans to group space
    fprintf('\n(mlrSphericalNormGroup) Converting %s data to %s template space ...\n',params.mrLoadRetSubjectIDs{iSubj},params.freesurferTemplateID);
    fsSphericalParams = freesurferSphericalNormalizationVolumes([],'justGetParams');
    fsSphericalParams.sourceVol = exportNames;
    fsSphericalParams.fsSourceSubj = params.freesurferSubjectIDs{iSubj};
    fsSphericalParams.fsSourceSurfSuffix = params.fsSubjectSurfSuffix{iSubj};
    fsSphericalParams.fsDestSubj = params.freesurferTemplateID;
    fsSphericalParams.destVol = convertNames(:,1);
    fsSphericalParams.interpMethod = 'linear';
    fsSphericalParams.outputBinaryData = false;  % UNCLEAR IF THIS SHOULD BE TRUE OR FALSE (ONLY MATTERS FOR ROIs)
    fsSphericalParams.dryRun = params.dryRun;
    fsSphericalParamsOut = freesurferSphericalNormalizationVolumes(fsSphericalParams);
    if params.combineLeftAndRight
      fprintf('\n(mlrSphericalNormGroup) Converting %s data to %s template space ...\n',params.mrLoadRetSubjectIDs{iSubj},params.lrFlippedFreesurferTemplateID);
      fsSphericalParams.fsDestSubj = params.lrFlippedFreesurferTemplateID;
      fsSphericalParams.destVol = convertNames(:,2);
      freesurferSphericalNormalizationVolumes(fsSphericalParams);
    end
  else
    return
  end
  
  if ~params.dryRun

    % delete subject-space exported NIFTI files
    for iFile = 1:nOverlays+nROIs
      delete(exportNames{iFile});
    end

    % Concatenate normalized subject data
    fprintf('\n(mlrSphericalNormGroup) Concatenating transformed data for subject %s...\n',params.mrLoadRetSubjectIDs{iSubj});
    for iGroup = 1:nTemplateGroups + (nROIs>0)
      if iGroup > nTemplateGroups
        niftiFiles = nOverlays+ (1:nROIs);
        nScans(iGroup) = 1+2*params.combineLeftAndRight;
      else
        niftiFiles = find(params.templateOverlayGroupNums==iGroup);
        nScans(iGroup) = 1+params.combineLeftAndRight;
      end
      for iFile = niftiFiles
        cNiftiConcat(iGroup,1) = cNiftiConcat(iGroup,1)+1;
        if iSubj==1 && iFile == niftiFiles(1)
          if iGroup>nTemplateGroups
            templateScanName{iGroup,1} = sprintf('Concatenation of ROIs %s', mat2str(roiExportList));
          else
            templateScanName{iGroup,1} = sprintf('Concatenation of overlays-analyses-groups %s-%s-%s', ...
              mat2str(params.subjectOverlays(niftiFiles)),mat2str(params.subjectOverlayAnalyses(niftiFiles)),mat2str(params.subjectOverlayGroups(niftiFiles)));
          end
          templateScanName{iGroup,1} = fixBadChars(templateScanName{iGroup,1}); % replace spaces by underscores
          scanFileName{iGroup,1} = fullfile(templateTseriesFolder{iGroup},[templateScanName{iGroup,1} '.nii']);
          if params.combineLeftAndRight
            templateScanName{iGroup,2} = [templateScanName{iGroup,1} '_LRcombined'];
            if iGroup>nTemplateGroups
              templateScanName{iGroup,3} = [templateScanName{iGroup,2} 'OnRight'];
              templateScanName{iGroup,2} = [templateScanName{iGroup,2} 'OnLeft'];
              scanFileName{iGroup,3} = fullfile(templateTseriesFolder{iGroup},[templateScanName{iGroup,3} '.nii']);
            end
            scanFileName{iGroup,2} = fullfile(templateTseriesFolder{iGroup},[templateScanName{iGroup,2} '.nii']);
          end
        end
        [data,hdr] = cbiReadNifti(convertNames{iFile,1});
        hdr.time_units = 'subjects/conditions'; % I dont think this is is doing anything
        cbiWriteNifti(scanFileName{iGroup,1},data,hdr,'',{[],[],[],cNiftiConcat(iGroup,1)});
        % for log file
        if iGroup > nTemplateGroups && params.combineLeftAndRight % for the ROI group, add the side to the name if combining between left and right
          logfile{iGroup,1}.overlay{cNiftiConcat(iGroup,1),1} = [sides{params.subjectROIsides(iFile-nOverlays)} templateOverlayNewNames{iFile}];
        else
          logfile{iGroup,1}.overlay{cNiftiConcat(iGroup,1),1} = templateOverlayNewNames{iFile};
        end
        logfile{iGroup,1}.subject{cNiftiConcat(iGroup,1),1} = params.mrLoadRetSubjectIDs{iSubj};
        if iGroup > nTemplateGroups % for the ROI group, add a side variable
          logfile{iGroup,1}.hemisphere{cNiftiConcat(iGroup,1),1} = sides{params.subjectROIsides(iFile-nOverlays)};
        end

        if params.combineLeftAndRight
          if iGroup > nTemplateGroups % for the ROI group, split the data between left and right hemispheres
             if params.subjectROIsides(iFile-nOverlays)==1
               normalScanNum = 2; % put normally-oriented left-hemisphere ROIs in the 2nd scan
               reverseScanNum = 3; % but reverse-oriented right-hemisphere ROIs in a 3rd scan
             else
               normalScanNum = 3; % for right-hemisphere ROIs
               reverseScanNum = 2; % do the reverse
             end
          else
            normalScanNum = 2;
            reverseScanNum = 2;
          end
          cNiftiConcat(iGroup,normalScanNum) = cNiftiConcat(iGroup,normalScanNum)+1;
          cbiWriteNifti(scanFileName{iGroup,normalScanNum},data,hdr,'',{[],[],[],cNiftiConcat(iGroup,normalScanNum)});
          logfile{iGroup,normalScanNum}.overlay{cNiftiConcat(iGroup,normalScanNum),1} = templateOverlayNewNames{iFile};
          logfile{iGroup,normalScanNum}.leftRightDirection{cNiftiConcat(iGroup,normalScanNum),1} = 'normal';
          logfile{iGroup,normalScanNum}.subject{cNiftiConcat(iGroup,normalScanNum),1} = params.mrLoadRetSubjectIDs{iSubj};
          if iGroup > nTemplateGroups % for the ROI group, add the side
            logfile{iGroup,normalScanNum}.hemisphere{cNiftiConcat(iGroup,normalScanNum),1} = sides{params.subjectROIsides(iFile-nOverlays)};
          end
          cNiftiConcat(iGroup,reverseScanNum) = cNiftiConcat(iGroup,reverseScanNum)+1;
          [data,hdr] = cbiReadNifti(convertNames{iFile,2});
          if iGroup <= nTemplateGroups
            if params.lrFlipTransform(iFile)
              data = params.lrFlipTransformFunctions{params.lrFlipTransform(iFile)}(data);
            end
          end
          hdr.time_units = 'subjects/conditions'; % I dont think this is is doing anything
          cbiWriteNifti(scanFileName{iGroup,reverseScanNum},data,hdr,'',{[],[],[],cNiftiConcat(iGroup,reverseScanNum)});
          logfile{iGroup,reverseScanNum}.overlay{cNiftiConcat(iGroup,reverseScanNum),1} = templateOverlayNewNames{iFile};
          logfile{iGroup,reverseScanNum}.leftRightDirection{cNiftiConcat(iGroup,reverseScanNum),1} = 'reversed';
          logfile{iGroup,reverseScanNum}.subject{cNiftiConcat(iGroup,reverseScanNum),1} = params.mrLoadRetSubjectIDs{iSubj};
          if iGroup > nTemplateGroups % for the ROI group, add the side
            logfile{iGroup,reverseScanNum}.hemisphere{cNiftiConcat(iGroup,reverseScanNum),1} = sides{params.subjectROIsides(iFile-nOverlays)};
          end
        end
        % delete temporary converted subject files
        for iCombine = 1:1+params.combineLeftAndRight
          delete(convertNames{iFile,iCombine});
        end
      end
    end
  end
end

if ~params.dryRun
  
  if params.cropScans
    for iGroup = 1:nTemplateGroups+(nROIs>0)
      for iScan = 1:nScans(iGroup)
        fprintf('\n(mlrSphericalNormGroup) Cropping scan %d\n',iScan);
        cropBox = [inf -inf;inf -inf;inf -inf];
        croppedScanFileName = fullfile(templateTseriesFolder{iGroup},[templateScanName{iGroup,iScan} '_cropped.nii']);
        scanHdr = cbiReadNiftiHeader(scanFileName{iGroup,iScan});
        for iVolume = 1:scanHdr.dim(5)
          data = cbiReadNifti(scanFileName{iGroup,iScan},{[],[],[],iVolume});
          [X,Y,Z] = ind2sub(scanHdr.dim(2:4)',find(~isnan(data)&data~=0)); % can probably do better than this by finding main axes
          cropBox(:,1) = min(cropBox(:,1), floor([min(X);min(Y);min(Z)]));  % of coordinates by svd and applying rotation before cropping
          cropBox(:,2) = max(cropBox(:,2), floor([max(X);max(Y);max(Z)]));  % but this would involve resampling the data
        end
        scanHdr.dim(2:4) = diff(cropBox,[],2)+1;
        cropXform = eye(4);
        cropXform(1:3,4) = -cropBox(:,1)+1;
        scanHdr.qform44 = cropXform\scanHdr.qform44;
        scanHdr.sform44 = cropXform\scanHdr.sform44;
        for iVolume = 1:scanHdr.dim(5)
          data = cbiReadNifti(scanFileName{iGroup,iScan},{cropBox(1,:),cropBox(2,:),cropBox(3,:),iVolume});
          cbiWriteNifti(croppedScanFileName,data,scanHdr,'',{[],[],[],iVolume});
        end
        movefile(croppedScanFileName,scanFileName{iGroup,iScan});
      end
    end
  end
  
  % initialize mrLoadRet view/import group
  fprintf('\n(mlrSphericalNormGroup) Importing group data into MLR folder %s\n',params.mrLoadRetTemplateID);
  cd(mrLoadRetTemplateFolder);
  if initMrLoadRetGroup
    scanList{1} = 1:nScans(1);
    subjectSessionParams = load(fullfile(params.studyDir,params.mrLoadRetSubjectIDs{1},'mrSession.mat'));
    [sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1');
    for iScan = scanList{1}
      groupParams.description{iScan} = templateScanName{1,iScan};
    end
    sessionParams.magnet = subjectSessionParams.session.magnet;
    sessionParams.coil = subjectSessionParams.session.coil;
    [sessionParams.pulseSequence,sessionParams.pulseSequenceText] = strtok(subjectSessionParams.session.protocol,':');
    mrInit(sessionParams,groupParams,'makeReadme=0','noPrompt=1');
    templateGroupNums{1} = 1;
    % load base anatomies:
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
  end
  
  % import (other) groups
  for iGroup = 1+initMrLoadRetGroup:nTemplateGroups+(nROIs>0)
    thisView = viewSet(thisView,'newGroup',templateGroupFolder{iGroup});
    templateGroupNums{iGroup} = viewGet(thisView,'groupNum',templateGroupFolder{iGroup});
    thisView = viewSet(thisView,'curGroup',templateGroupNums{iGroup});
    for iScan = 1:nScans(iGroup)
      [~,importParams] = importTSeries(thisView,[],'justGetParams=1','defaultParams=1',['pathname=' scanFileName{iGroup,iScan}]);
      importParams.description = templateScanName{iGroup,iScan};
      importParams.overwrite = 1;
      thisView = importTSeries(thisView,importParams);
    end
    scanList{iGroup} = viewGet(thisView,'nScans')-(nScans(iGroup):-1:1)+1;
  end
  
  for iGroup = 1:nTemplateGroups+(nROIs>0)
    % create log files indicating which volumes correspond to which overlays and subjects
    for iScan = 1:nScans(iGroup)
      logfileStruct = logfile{iGroup,iScan};
      logFilename = [templateScanName{iGroup,iScan} '.mat'];
      save(fullfile(mrLoadRetTemplateFolder,'Etc',logFilename),'-struct','logfileStruct');
      % link log file to concatenated scan
      fprintf('(mlrSphericalNormGroup) Linking %s to scan %d\n',logFilename,scanList{iGroup}(iScan));
      viewSet(thisView,'stimfilename',logFilename, scanList{iGroup}(iScan),templateGroupNums{iGroup});
    end
    
    if ( iGroup>nTemplateGroups && params.computeROIprobabilityMaps ) || params.computeGroupAverage(iGroup)
      thisView = viewSet(thisView,'curGroup',templateGroupNums{iGroup});
      [~,averageParams] = mlrGroupAverage(thisView,[],'justGetParams');
      averageParams.factors = {'overlay'};
      averageParams.scanList = scanList{iGroup};
      thisView = mlrGroupAverage(thisView,averageParams);
      thisView = viewSet(thisView,'curOverlay',1);
      thisView = viewSet(thisView,'clipAcrossOverlays',false);
    end
  end

  mrSaveView(thisView);
  deleteView(thisView);

end

if params.dryRun
  success = true;
elseif length(dir(temporaryTseriesFolder))==2
  rmdir(temporaryTseriesFolder);
  success = true;
else
  mrWarnDlg(sprintf('(mlrSphericalNormGroup) Temporary files left in %s',temporaryTseriesFolder));
end
