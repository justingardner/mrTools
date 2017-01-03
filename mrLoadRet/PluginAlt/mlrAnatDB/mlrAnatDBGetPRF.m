% mlrAnatDBGetPRF.m
%
%        $Id:$ 
%      usage: [pRF sessionName pRFGroup pRFAnalName] = mlrAnatDBGetPRF(sid)
%         by: justin gardner
%       date: 01/03/17
%    purpose: get pRF for subject: sid (sid is number or view)
%             This function will also add to the pRF structure
%             the scan2mag xform and the scanNum for which scan
%             to use of the pRF
%
%       e.g.: pRF = mlrAnatDBGetPRF('s0389')
%
%             % without trying to pull directories
%             pRF = mlrAnatDBGetPRF('s0389','noPull=1')
%    
%             % set the name of the analysis folder to look for
%             pRF = mlrAnatDBGetPRF(389,'analysisName=pRFAnal');
%
function [pRF sessionName pRFGroup pRFAnalName] = mlrAnatDBGetPRF(sid,varargin)

% check arguments
if nargin <1
  help mlrAnatDBGetPRF
  return
end

% default returnarguments
pRF = [];sessionName = [];pRFGroup = [];pRFAnalName = [];

% parse rest of arguments
getArgs(varargin,{'noPull=0','analysisName=pRFAnal'});

% get subjectID
subjectID = mlrAnatDBSubjectID(sid);
if isempty(subjectID),return,end

% get repo
[localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID,'noPull',noPull);

% check that we got the localRepoLargeFiles
if isempty(localRepoLargeFiles)
  disp(sprintf('(mlrAnatDBGetPRF) Could not get mlrAnatDB repo with localizer for: %s',subjectID));
  return
end

% look for pRF localizers
localizerDir = fullfile(localRepoLargeFiles,'localizers');
if ~isdir(localizerDir)
  disp(sprintf('(mlrAnatDBGetPRF) Could not find localizers directory: %s',localizerDir));
  return
end

% find localizers in that directory
localizerList = dir(fullfile(localizerDir,sprintf('%s*',subjectID)));
if length(localizerList) == 0
  disp(sprintf('(mlrAnatDBGetPRF) Could not find any localizer directories starting with %s in %s',subjectID,localizerDir));
  return
end

% if there is more than one localizer then
% ask user which one they want to use
if length(localizerList) > 1
  paramsInfo = {{'localizerName',{localizerList(:).name},'Choose which localizer to get the pRF data from'}};
  params = mrParamsDialog(paramsInfo,'Choose which localizer to use');
  if isempty(params),return,end
  pRFDir = params.localizerName;
else
  pRFDir = localizerList(1).name;
end

% swap to pRF session
sessionName = fullfile(localizerDir,pRFDir);
currentSession = mlrSwapSession(sessionName);

% look for pRFAnal directory, use system command find
[status pRFAnalDirs] = system(sprintf('find %s -name %s -type d -print',sessionName,analysisName));

% check for error
if status ~= 0
  disp(sprintf('(mlrAnatDBGetPRF) Could not run system command: find'));
  mlrSwapSession(currentSession);
  return
end

% check for any pRF analyses
if isempty(pRFAnalDirs)
  disp(sprintf('(mlrAnatDBGetPRF) Could not find any pRF analyses in: %s',pRFDir));
  mlrSwapSession(currentSession);
  return
end

% check for all pRFANal directories
[pRFAnalDir pRFAnalDirs] = strtok(pRFAnalDirs);
% get the group where that directory is
pRFGroup = getLastDir(fileparts(pRFAnalDir));
% check if there is more than one
if ~isempty(pRFAnalDirs)
  pRFGroupList{1} = pRFGroup;
  while ~isempty(pRFAnalDirs)
    [pRFAnalDir pRFAnalDirs] = strtok(pRFAnalDirs);
    if ~isempty(pRFAnalDir)
      pRFGroupList{end+1} = getLastDir(fileparts(pRFAnalDir));
    end
  end
  if length(pRFGroupList) > 1
    paramsInfo = {{'pRFGroup',pRFGroupList,'Choose which group you want to get pRF analysis from'}};
    params = mrParamsDialog(paramsInfo,'Choose pRF Group');
    if isempty(params)
      mlrSwapSession(currentSession);
      return
    end
    pRFGroup = params.pRFGroup;
  else
    pRFGroup = pRFGroupList{1};
  end
end

% find all pRF analyses in that group
pRFAnalyses = dir(fullfile(sessionName,pRFGroup,analysisName,'*.mat'));

% if there is more than one analysis then have the user choose
if length(pRFAnalyses) > 1
  paramsInfo = {{'pRFAnalysis',{pRFAnalyses(:).name},'Choose pRF Analysis'}};
  params = mrParamsDialog(paramsInfo,'Choose which pRF Analysis to load');
  if isempty(params)
    mlrSwapSession(currentSession);
    return;
  end
  pRFAnalysis = params.pRFAnalysis;
else
  pRFAnalysis = pRFAnalyses.name;
end
pRFAnalName = fullfile(analysisName,pRFAnalysis);

% ok, ready to load now
v = newView;
v = viewSet(v,'curGroup',pRFGroup);
v = loadAnalysis(v,pRFAnalName);
pRF = viewGet(v,'analysis');

% check to see if it has been run on multiple scans
% if so, ask user which one
scanNum = pRF.params.scanNum;
if length(scanNum) > 1
  paramsInfo = {{'scanNum',num2cell(scanNum),'Multiple scans have been run for this pRF. Choose which scan number to use.'}};
  params = mrParamsDialog(paramsInfo,'Choose which scan to use');
  if isempty(params)
    deleteView(v);
    mlrSwapSession(currentSession);
  end
  scanNum = params.scanNum;
end
pRF.scanNum = scanNum;

% also get scan2mag and stick that on the pRF
pRF.scan2mag = viewGet(v,'scan2mag',scanNum);
pRF.scanDims = viewGet(v,'scanDims',scanNum);

% delete view
deleteView(v);

% and switch back
mlrSwapSession(currentSession);
