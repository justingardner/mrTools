% motionCompGUImrParams.m
%
%        $Id$
%      usage: motionCompGUImrParams(varargin)
%         by: justin gardner
%       date: 08/03/07
%    purpose: replace the old GUI with an mrParams dialog,
%             initial portion that parses input arguments taken
%             from motionComp.m
%
% params = motionCompGUI('groupName','Raw');
% params = motionCompGUI('params',params);
% params = motionCompGUI('groupName','Raw','params',params);
function mrParams = motionCompGUImrParams(varargin)

mrParams = [];
% check arguments
if ~any(nargin == [2 4])
  help motionCompGUImrParams
  return
end

% Parse varargin
for index = 1:2:length(varargin)
  field = varargin{index};
  val = varargin{index+1};
  switch field
    case 'groupName'
      groupName = val;
    case 'params'
      params = val;
    otherwise
      mrWarnDlg('Invalid initialization argument')
  end
end

% Error if neither params nor groupName specified
if ieNotDefined('params') & ieNotDefined('groupName')
  mrErrorDlg('Must initialize with either params or groupName.');
end

% If groupName not passed then set it according to params.groupName (which
% we now know must exist when groupName does not).
if ieNotDefined('groupName')
  if ~isfield(params,'groupName')
    mrErrorDlg('Must initialize with params.groupName.');
  end
  groupName = params.groupName;
end

% Group names, group number, and nScans
groupNames = viewGet([],'groupNames');
if ~any(strcmp('MotionComp',groupNames))
  groupNames{length(groupNames)+1} = 'MotionComp';
end
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
  mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
  params = motionCompReconcileParams(groupName);
else
  params = motionCompReconcileParams(groupName,params);
end

% the code below here is converted to work with mrParamsDialog
% instead of initializing handles it initializes the paramsInfo.
%
% Initialize handles. Store the parameters in the GUI until user
% clicks the ok or cancel buttons. Then return the params from
% motionCompGUI_OutputFcn.
paramsInfo = {};
%paramsInfo{end+1} = {'groupName',params.groupName};
%paramsInfo{end+1} = {'crop',params.crop};

% Initialize include (targetScans)
includeScans = cell(1,nScans);
for i = 1:length(params.targetScans)
  includeScans{params.targetScans(i)} = 1;
end
paramsInfo{11} = {'scanNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'Scan selector. Use this to choose which scans to include'};
paramsInfo{14} = {'include',includeScans,'type=checkbox','group=scanNum'};

% Initialize tseriesfiles
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end
paramsInfo{12} = {'tseriesfiles',tseriesfiles,'group=scanNum','type=String','editable=0','Filename of scan'};

% Initialize descriptions
descriptions = cell(1,nScans);
for scan = 1:nScans
  m = find(scan == params.targetScans);
  if m
    descriptions{scan} = params.descriptions{m};
  else
    descriptions{scan} = ['Motion compensation of ',groupName,' scan ',num2str(scan)];
  end
end
paramsInfo{13} = {'descriptions',descriptions,'group=scanNum','type=String','editable=0','Scan description'};

% Initialize group popup
groupNames = {'New',groupNames{:}};
motionCompGroupName = params.motionCompGroupName;
motionCompGroupNum = find(strcmp(motionCompGroupName,groupNames));
if isempty(motionCompGroupNum)
  groupNames{length(groupNames)+1} = motionCompGroupName;
  motionCompGroupNum = length(groupNames);
end
paramsInfo{1} = {'motionCompGroupName',putOnTopOfList(motionCompGroupName,groupNames),'Group name to put the motion compensated scans into'};

% Initialize base scan popup
paramsInfo{3} = {'baseScan',params.baseScan,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans)};

% Initialize base frame popup
baseFrameStrings = {'first','last','mean'};
paramsInfo{4} = {'baseFrame',putOnTopOfList(params.baseFrame,baseFrameStrings)};

% Initialize interp method popup
interpMethodStrings = {'nearest','linear','cubic','spline'};
paramsInfo{2} = {'interpMethod',putOnTopOfList(params.interpMethod,interpMethodStrings)};

% Initialize checkboxes
sliceTimeStrings = {'beginning of TR','middle of TR','end of TR'};
paramsInfo{7} = {'sliceTimeCorrection',params.sliceTimeCorrection,'type=checkbox'};
paramsInfo{8} = {'sliceTimeString',putOnTopOfList(params.sliceTimeString,sliceTimeStrings)};
paramsInfo{9} = {'robust',params.robust,'type=checkbox'};
paramsInfo{10} = {'correctIntensityContrast',params.correctIntensityContrast,'type=checkbox'};

% Initialize niters
paramsInfo{5} = {'niters',params.niters,'incdec=[-1 1]','minmax=[0 inf]'};

paramsInfo{6} = {'crop',params.crop,'callback',@thisSelectCropRegion,'buttonString=Set crop region','passParams=1','type=pushbutton'};

% put up dialog
mrParams = mrParamsDialog(paramsInfo,'Set motionComp parameters');

if ~isempty(mrParams)
  % do some clean up on the parameters
  % to make it return exactly what the old GUI returned
  mrParams.groupName = params.groupName;
  if isempty(mrParams.crop)
    mrParams.crop = [];
  end
  mrParams.targetScans = find(mrParams.include);
  tseriesfiles = {}; descriptions = {};
  for i = 1:length(mrParams.targetScans)
    tseriesfiles{i} = mrParams.tseriesfiles{mrParams.targetScans(i)};
    descriptions{i} = mrParams.descriptions{mrParams.targetScans(i)};
  end
  mrParams.tseriesfiles = tseriesfiles;
  mrParams.descriptions = descriptions;
  % remove extraneous fields
  mrParams = rmfield(mrParams,'include');
  mrParams = rmfield(mrParams,'paramInfo');
  mrParams = rmfield(mrParams,'scanNum');
  % order fields
  mrParams = orderfields(mrParams,{'groupName','baseScan','baseFrame','sliceTimeCorrection','sliceTimeString','robust','correctIntensityContrast','crop','niters','motionCompGroupName','interpMethod','targetScans','tseriesfiles','descriptions'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for setting crop region, gets the correct volume
% and passes it to selectCropRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crop = thisSelectCropRegion(params)

params.baseScan
view = newView('Volume');
% Load first frame of base scan
volume = loadTSeries(view,params.baseScan,'all',1);
% and get crop region
crop = selectCropRegion(volume);
