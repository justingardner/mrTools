% importGroupScans.m
%
%      usage: importGroupScans(params)
%         by: justin gardner
%       date: 04/11/07
%    purpose: import scans into a group in the current mrTools sesssion, from a group in a different (or the same) mrTools session.
%             params is a structure with the following fields:
%             - fromSession:  path to the source mrTools session
%             - fromGroup:    Name of the source group
%             - scanList:     Scan numbers of scans to import
%             - toGroup:      Name of the destination group
%             - linkFiles:    Whether to link scan files instead of copying them (not availale on Windows)
%             - hardLink:     Whether to use a hard link instead of a soft link (not availale on Windows)
%
function importGroupScans(params)

% check arguments
if ~ismember(nargin, [0, 1])
  help importGroupScans
  return
end

% new view
toView = newView;

if ieNotDefined('params')
  params = struct();
end

if fieldIsNotDefined(params,'fromSession')
  % go find the group that user wants to load here
  params.fromSession = uigetdir(viewGet(toView,'homeDir'),'Select session you want to import from');
  if (params.fromSession==0)
    deleteView(toView);
    return
  end
end


% now look for that sessions mrSession
mrSessionPath = fullfile(params.fromSession,'mrSession.mat');
if ~mlrIsFile(mrSessionPath)
  disp(sprintf('(importGroupScans) Could not find mrSession in %s',fileparts(params.fromSession)));
  disp(sprintf('                   Make sure you clicked on the directory'));
  disp(sprintf('                   with the mrSession.mat file (not the group'))
  disp(sprintf('                   directory you wanted to import)'))
  deleteView(toView);
  return
end

% check for MLR 4 session
mrSession = load(mrSessionPath);
if ~isfield(mrSession,'session') || ~isfield(mrSession,'groups')
  mrWarnDlg(sprintf('(importGroupScans) Unknown format for mrSession in %s',fileparts(params.fromSession)));
  deleteView(toView);
  return
end
clear mrSession

% get info from destination group
toHomeDir = viewGet(toView,'homeDir');
toGroupNames = viewGet(toView,'groupNames');

% now we are going to have to clear the MLR
% variable so that we can get info on the from groups
% we will then set back to the old MLR. Note that
% while we have switch the MLR session we cannot
% get info from the toView
switchSession(params.fromSession);
fromView = newView;

if fieldIsNotDefined(params,'fromGroup')

  % get the groups in the import session
  for gNum = 1:viewGet(fromView,'numGroups')
    fromGroups{gNum} = sprintf('%s:%s (%i scans)',getLastDir(params.fromSession),viewGet(fromView,'groupName',gNum),viewGet(fromView,'numScans',gNum));
  end

  % get from which and to which group we are doing
  paramsInfo = {...
      {'fromGroup',fromGroups,'type=popupmenu','The group to import from'},...
      {'toGroup',toGroupNames,'type=popupmenu','The group to import into'},...
      {'linkFiles',~ispc,'type=checkbox',sprintf('enable=%d',~ispc),'Link rather than copy the files (Mac/Linux only). This will make a soft link rather than copying the files which saves disk space.'},...
      {'hardLink',0,'type=checkbox',sprintf('enable=%d',~ispc),'contingent=linkFiles','(Mac/Linux only) Use hard links when linking files instead of soft links.'}};
  
  inputParams = params;
  params = mrParamsDialog(paramsInfo);
  if isempty(params)
    switchSession;
    deleteView(toView);
    return
  end
  
  params.fromSession = inputParams.fromSession;
  fromGroupNum = find(strcmp(params.fromGroup,fromGroups));
  
else
  fromGroupNum = viewGet(fromView,'groupNum',params.fromGroup);
end

% get whether to link or not
linkType = 0;
if params.linkFiles
  if ispc
    mrWarnDlg('(importGroupScans) Linking scan files is not implemented on Windows');
    switchSession;
    deleteView(toView);
    return
  end

  % for hard links, pass 2
  if params.hardLink
    linkType = 2;
  else
    linkType = 1;
  end
end

% now set up some variables
fromGroup = viewGet(fromView,'groupName',fromGroupNum);
toGroup = params.toGroup;
fromDir = fullfile(fullfile(params.fromSession,fromGroup),'TSeries');
if ~isdir(fromDir)
  mrWarnDlg(sprintf('(importGroupScans) Could not find directory %s',fromDir));
  switchSession;
  deleteView(toView);
  return
end
toDir = fullfile(fullfile(toHomeDir,toGroup),'TSeries');
if ~isdir(toDir)
  mrWarnDlg(sprintf('(importGroupScans) Could not find directory %s',toDir));
  % switch back to old one
  switchSession;
  deleteView(toView);
  return
end
fromName = getLastDir(params.fromSession);


% set the group
fromView = viewSet(fromView,'curGroup',fromGroupNum);

if fieldIsNotDefined(params,'scanList')
  % choose the scans to import
  params.scanList = selectInList(fromView,'scans','Choose scans to import');
end

if isempty(params.scanList)
  switchSession;
  deleteView(toView);
  return
end

% get the scan and aux paramters for the chosen scans
for i = 1:length(params.scanList)
  fromScanParams(i) = viewGet(fromView,'scanParams',params.scanList(i));
  fromAuxParams(i) = viewGet(fromView,'auxParams',params.scanList(i));
  % go through auxParams and get all fields
  if ~isempty(fromAuxParams(i))
    % get names of aux params
    fromAuxParamsNames{i} = fieldnames(fromAuxParams(i));
    % dont worry about stimFileName
    fromAuxParamsNames{i} = setdiff(fromAuxParamsNames{i},'stimFileName');
    % look for a few other ones
    fromAuxParamsNames{i} = {fromAuxParamsNames{i}{:} 'volTrigRatio' 'tSense' 'fidFilename'};
    % get values
    for iAuxParam = 1:length(fromAuxParamsNames{i})
      fromAuxParamsValues{i}{iAuxParam} = viewGet(fromView,'auxParam',fromAuxParamsNames{i}{iAuxParam});
    end
  end
end


% get the stimfiles for the selected scans
for scanNum = 1:length(params.scanList)
  stimFileName{scanNum} = viewGet(fromView,'stimFileName',params.scanList(scanNum));
end

% now switch back to old MLR session
toView = switchSession;

% set the group
toView = viewSet(toView,'currentGroup',toGroup);

% now cycle over all scans in group
mlrDispPercent(-inf,'Copying group scans');
r = 0;
for scanNum = 1:length(fromScanParams)
  startTime = clock;
  % get the fromFilename
  fromFilename = fullfile(fromDir,fromScanParams(scanNum).fileName);
  % get the scan params for this scan
  toScanParams = fromScanParams(scanNum);
  toScanParams.description = sprintf('%s:%s',fromName,fromScanParams(scanNum).description);
  % and remove any reference to originalFileName/GroupName
  toScanParams.originalFileName = [];
  toScanParams.originalGroupName = [];
  % clear filename, so that it gets a new unique filename
  toScanParams.fileName = [];
  % and now add the scan to our group
  saveNewTSeries(toView,fromFilename,toScanParams,[],linkType);
  toScanNum = viewGet(toView,'nScans');
  % copy the stimFile over
  toStimFileNames = {};
  for stimFileNum = 1:length(stimFileName{scanNum})
    % get the from and to stim file names
    fromStimFileName = fullfile(params.fromSession, 'Etc', getLastDir(stimFileName{scanNum}{stimFileNum}));
    toStimFileName = fullfile(viewGet(toView,'EtcDir'),getLastDir(stimFileName{scanNum}{stimFileNum}));
    % if it doesn't exist already, then copy it over
    if mlrIsFile(toStimFileName)
      disp(sprintf('(importGroupScans) Stimfile %s already exists',toStimFileName));
    else
      %      system(sprintf('cp %s %s',fromStimFileName,toStimFileName));
      if ~mlrIsFile(fromStimFileName)
	mrWarnDlg(sprintf('(importGroupScans) !!! Missing stimfile: %s !!!',fromStimFileName));
      else
	copyfile(fromStimFileName,toStimFileName);
      end
    end
    toStimFileNames{stimFileNum} = getLastDir(toStimFileName);
  end
  % and set the stimfiles in the view
  toView = viewSet(toView,'stimFileName',toStimFileNames,toScanNum);
  % now go through and add any non-empty auxparams
  for iAuxParams = 1:length(fromAuxParamsNames{scanNum})
    if ~isempty(fromAuxParamsValues{scanNum}{iAuxParams})
      toView = viewSet(toView,'auxParam',fromAuxParamsNames{scanNum}{iAuxParams},fromAuxParamsValues{scanNum}{iAuxParams},toScanNum);
    end
  end

  % pause for 1 second (so that we don't end up writing scans out with the same exact timestamp
  if etime(clock,startTime) < 1
    disp(sprintf('(importGroupScans) Pause for one second to avoid having same exact timestamps'));
    pause(1);
  end
  mlrDispPercent(scanNum/length(fromScanParams));
end
mlrDispPercent(inf);

deleteView(toView);
saveSession;


%%%%%%%%%%%%%%%%%%%%%%%
%%   switchSession   %%
%%%%%%%%%%%%%%%%%%%%%%%
function v = switchSession(pathStr)

% switch to the MLR session found at "inputParams.fromSession"
if (nargin == 1)
  % switch the path and globals
  mrGlobals;
  global oldMLR;
  oldMLR = MLR;
  clear global MLR;
  global oldPWD
  oldPWD = pwd;
  cd(pathStr);
  v = newView;
% switch back to old one
else
  % switch back to old one
  global oldPWD
  cd(oldPWD);
  clear global MLR;
  global MLR;
  global oldMLR;
  MLR = oldMLR;
  clear global oldMLR
end


