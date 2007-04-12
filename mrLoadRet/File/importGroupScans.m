% importGroupScans.m
%
%      usage: importGroupScans()
%         by: justin gardner
%       date: 04/11/07
%    purpose: 
%
function retval = importGroupScans()

% check arguments
if ~any(nargin == [0])
  help importGroupScans
  return
end

% new view
mrGlobals;
v = newView('Volume');

% go find the group that user wants to load here
pathStr = uigetdir(MLR.homeDir,'Select session you want to import from');
if (pathStr==0),return,end

% now look for that sessions mrSession
mrSessionPath = fullfile(pathStr,'mrSession.mat');
if ~isfile(mrSessionPath)
  disp(sprintf('(importGroupScans) Could not find mrSession in %s',getPath(pathStr)));
  return
end
mrSession = load(mrSessionPath);
if ~isfield(mrSession,'session') || ~isfield(mrSession,'groups')
  disp(sprintf('(importGroupScans) Unknown format for mrSession in %s',getPath(pathStr)));
  return
end

% get the groups in the import session
for gNum = 1:length(mrSession.groups)
  fromGroups{gNum} = sprintf('%s:%s (%i scans)',getLastDir(pathStr),mrSession.groups(gNum).name,length(mrSession.groups(gNum).scanParams));;
end

% get from which and to which group we are doing
paramsInfo = {...
    {'fromGroup',fromGroups,'The group to import from'},...
    {'toGroup',viewGet(v,'groupNames'),'The group to import into'}};
    
params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

% now set up some variables
fromGroupNum = find(strcmp(params.fromGroup,fromGroups));
fromGroup = mrSession.groups(fromGroupNum).name;
toGroupNum = viewGet(v,'groupNum',params.toGroup);
toGroup = params.toGroup;
fromDir = fullfile(fullfile(pathStr,fromGroup),'TSeries');
toDir = fullfile(fullfile(viewGet(v,'homeDir'),toGroup),'TSeries');
fromScanParams = mrSession.groups(fromGroupNum).scanParams;
fromAuxParams = mrSession.groups(fromGroupNum).auxParams;
fromName = getLastDir(pathStr);

% set the group
v = viewSet(v,'currentGroup',toGroupNum);

% now we are going to have to clear the MLR
% variable so that we can get info on the from groups
% we will then set back to the old MLR
oldMLR = MLR;
clear global MLR;
currentDir = pwd;
cd(pathStr);
fromView = newView('Volume');
for scanNum = 1:length(fromScanParams)
  stimFileName{scanNum} = viewGet(fromView,'stimFileName',scanNum);
end
% now switch back to old one
cd(pwd);
clear global MLR;
global MLR;
MLR = oldMLR;

% now cycle over all scans in group
for scanNum = 1:length(fromScanParams)
  % make sure there is a file where we think there should be
  fromFilename = fullfile(fromDir,fromScanParams(scanNum).fileName);
  toFilename = fullfile(toDir,fromScanParams(scanNum).fileName);
  if ~isfile(fromFilename)
    disp(sprintf('(importGroupScans) Could not find %s',fromFilename));
    return
  end
  % load up the file
  [data hdr] = cbiReadNifti(fromFilename);
  % write it to our group
  cbiWriteNifti(toFilename,data,hdr);
  % add where this scan is from into description
  toScanParams = fromScanParams(scanNum);
  toScanParams.description = sprintf('%s:%s',fromName,fromScanParams(scanNum).description);
  % and now add the scan to our group
  v = viewSet(v,'newScan',toScanParams);
  % copy the stimFile over
  for stimFileNum = 1:length(stimFileName{scanNum})
    % get the from and to stim file names
    fromStimFileName = stimFileName{scanNum}{stimFileNum};
    toStimFileName = fullfile(viewGet(v,'EtcDir'),getLastDir(stimFileName{scanNum}{stimFileNum}));
    % if it doesn't exist already, then copy it over
    if isfile(toStimFileName)
      disp(sprintf('(importGroupScans) File %s already exists',toStimFileName));
    else
      system(sprintf('cp %s %s',fromStimFileName,toStimFileName));
    end
    % and set it in the view
    v = viewSet(v,'stimFileName',getLastDir(toStimFileName),scanNum);
  end
end

