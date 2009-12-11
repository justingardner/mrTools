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
% params = motionCompGUImrParams('groupName','Raw');
% params = motionCompGUImrParams('params',params);
% params = motionCompGUImrParams('groupName','Raw','params',params);
function mrParams = motionCompGUImrParams(varargin)

mrParams = [];
% check arguments
if ~any(nargin == [2 4 6 8])
  help motionCompGUImrParams
  return
end

defaultParams = 0;

% Parse varargin
for index = 1:2:length(varargin)
  field = varargin{index};
  val = varargin{index+1};
  switch field
    case 'groupName'
      groupName = val;
    case 'params'
      params = val;
    case 'defaultParams'
      defaultParams = val;
    case 'scanList'
      scanList = val;
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

% Initialize include (targetScans)
includeScans = cell(1,nScans);
for i = 1:length(params.targetScans)
  includeScans{params.targetScans(i)} = 1;
end

% Initialize tseriesfiles
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end

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

% Initialize group popup
groupNames = {'New',groupNames{:}};
motionCompGroupName = params.motionCompGroupName;
motionCompGroupNum = find(strcmp(motionCompGroupName,groupNames));
if isempty(motionCompGroupNum)
  groupNames{length(groupNames)+1} = motionCompGroupName;
  motionCompGroupNum = length(groupNames);
end

% Initialize base frame popup
baseFrameStrings = {'first','last','mean'};
% Initialize interp method popup
interpMethodStrings = {'nearest','linear','cubic','spline'};
% Initialize sliceTime strings
sliceTimeStrings = {'beginning of TR','middle of TR','end of TR'};

% set up params
paramsInfo = {};
paramsInfo{1} = {'motionCompGroupName',putOnTopOfList(motionCompGroupName,groupNames),'Group name to put the motion compensated scans into'};
paramsInfo{2} = {'interpMethod',putOnTopOfList(params.interpMethod,interpMethodStrings), 'Interpolation method used to warp the images'};
paramsInfo{3} = {'baseScan',params.baseScan,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans), 'Specifies the scan that everything else will be aligned to'};
paramsInfo{4} = {'baseFrame',putOnTopOfList(params.baseFrame,baseFrameStrings), 'Specifies which frame (or mean) to align to.  First frame can be bad choice b/c of different T2 contrast'};
paramsInfo{5} = {'niters',params.niters,'incdec=[-1 1]','minmax=[0 inf]', 'How many iterations to estimate the optimal transform (use more for high-res images)'};
paramsInfo{6} = {'crop',params.crop,'callback',@thisSelectCropRegion,'buttonString=Set crop region','passParams=1','type=pushbutton', 'Crop the images.  This affects the intensity/contrast correction.  Important for high-res images'};
paramsInfo{7} = {'sliceTimeCorrection',params.sliceTimeCorrection,'type=checkbox', 'Apply slice time correction along with motion compensation.  Not appropriate for 3D scans.'};
paramsInfo{8} = {'sliceTimeString',putOnTopOfList(params.sliceTimeString,sliceTimeStrings),'Which point in time the slices should be aligned to. May loose first and last frames'};
paramsInfo{9} = {'driftCorrection', params.driftCorrection, 'type=checkbox', 'Correction for fluctuations in mean intensity over time. This divides each frame by the average over the (cropped) volume.'};
paramsInfo{10} = {'gradIntensityCorrection',params.gradIntensityCorrection,'type=checkbox', 'Compensation for intensity gradient before motion estimation'};
paramsInfo{11} = {'tSmooth', params.tSmooth, 'incdec=[-1 1]', 'minmax=[1 10]', 'How much temporal smoothing.  Only applied to estimate head motion, not to final time series'};
paramsInfo{12} = {'robust',params.robust,'type=checkbox', 'Robust contrast estimator, should be used if images are noisy with lots of outliers'};
paramsInfo{13} = {'include',includeScans,'type=checkbox','group=scanNum', 'Check this to include a particular scan, uncheck to skip'};
paramsInfo{14} = {'scanNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'Scan selector. Use this to choose which scans to include'};
paramsInfo{15} = {'tseriesfiles',tseriesfiles,'group=scanNum','type=String','editable=0','Filename of scan'};
paramsInfo{16} = {'descriptions',descriptions,'group=scanNum','type=String','editable=0','Scan description'};

% put up dialog
if defaultParams
  mrParams = mrParamsDefault(paramsInfo);
  % if a scan list has been specified then
  % only select the scans from the scan list
  if ~ieNotDefined('scanList') 
    if (min(scanList) < 1) || (max(scanList) > nScans)
      mrErrorDlg('(motionCompGUImrParams) scanList include scans that do not exist');
    end
    mrParams.include(:) = 0;
    mrParams.include(scanList) = 1;
  end
else
  mrParams = mrParamsDialog(paramsInfo,'Set motionComp parameters');
end

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
  mrParams = orderfields(mrParams);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for setting crop region, gets the correct volume
% and passes it to selectCropRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crop = thisSelectCropRegion(params)

view = newView;
% Load first frame of base scan
volume = loadTSeries(view,params.baseScan,'all',1);
% and get crop region
crop = selectCropRegion(volume);
