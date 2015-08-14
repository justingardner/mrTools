% pRF.m
%
%        $Id:$ 
%      usage: pRF(v,params,varargin)
%         by: justin gardner
%       date: 11/20/11
%    purpose: compute pRF analysis on MLR data
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = pRF(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = pRF(v,[],'justGetParams=1');
%
function [v d] = pRF(v,params,varargin)

% check arguments
if nargin < 1
  help pRF
  return
end

d = [];
% a version number in case we make major changes
pRFVersion = 1;

% params defaults to empty
if nargin < 2,params =[];end

% other arguments
justGetParams=[];defaultParams=[];scanList=[];
groupNum=[];
getArgs(varargin,{'justGetParams=0','defaultParams=0','scanList=[]','groupNum=[]'});

% first get parameters
if isempty(params)
  % get group
  if isempty(groupNum),groupNum = viewGet(v,'curGroup');end
  % put up the gui
  params = pRFGUI('v',v,'groupNum',groupNum,'defaultParams',defaultParams,'scanList',scanList);
end

% just return parameters
if justGetParams,d = params;return,end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

% Abort if params empty
if isempty(params),return,end

% set the group
v = viewSet(v,'curGroup',params.groupName);

% create the parameters for the r2 overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'pRF';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(v,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(v,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'pRFPlot';
r2.mergeFunction = 'defaultMergeParams';

% create the parameters for the polarAngle overlay
polarAngle = r2;
polarAngle.name = 'polarAngle';
polarAngle.range = [-pi pi];
polarAngle.clip = [-pi pi];
polarAngle.colormapType = 'normal';
polarAngle.colormap = hsv(256);

% create the parameters for the eccentricity overlay
eccentricity = r2;
eccentricity.name = 'eccentricity';
eccentricity.range = [0 15];
eccentricity.clip = [0 inf];
eccentricity.colormapType = 'normal';
eccentricity.colormap = copper(256);

% create the paramteres for the rfHalfWidth overlay
rfHalfWidth = r2;
rfHalfWidth.name = 'rfHalfWidth';
rfHalfWidth.range = [0 15];
rfHalfWidth.clip = [0 inf];
rfHalfWidth.colormapType = 'normal';
rfHalfWidth.colormap = pink(256);

% get number of workers 
nProcessors = mlrNumWorkers;

% code snippet for clearing precomputed prefit
%global gpRFFitStimImage;gpRFFitStimImage = [];

dispHeader
disp(sprintf('(pRF) Running on scans %s:%s (restrict %s)',params.groupName,num2str(params.scanNum,'%i '),params.restrict ));

for scanNum = params.scanNum
  % see how long it took
  tic;
  
  % get voxels that we are restricted to
  [x y z] = getVoxelRestriction(v,params,scanNum);
  if isempty(x)
    disp(sprintf('(pRF) No voxels to analyze with current restriction'));
    return
  end
  % load all time series, we do this to pass inot pRFFit. Generally
  % the purpose here is that if we run on distributed computing, we
  % can't load each voxel's time series one at a time. If this is
  % too large for memory then you can comment this out and not
  % pass it into pRFFit and pRFFit will load the tSeries itself
  loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
  loadROI.coords(1,:) = x;
  loadROI.coords(2,:) = y;
  loadROI.coords(3,:) = z;
  loadROI = loadROITSeries(v,loadROI,scanNum,params.groupName);
  x = loadROI.scanCoords(1,:)';
  y = loadROI.scanCoords(2,:)';
  z = loadROI.scanCoords(3,:)';
  n = length(x);

  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum);
  
  % init overlays
  r2.data{scanNum} = nan(scanDims);
  polarAngle.data{scanNum} = nan(scanDims);
  eccentricity.data{scanNum} = nan(scanDims);
  rfHalfWidth.data{scanNum} = nan(scanDims);

  % default all variables that will be returned
  % by pRFFIt, so that we can call it the
  % second time and save some time
  concatInfo = [];
  stim = [];
  
  % save pRF parameters
  pRFAnal.d{scanNum}.ver = pRFVersion;
  pRFAnal.d{scanNum}.linearCoords = sub2ind(scanDims,x,y,z);
  pRFAnal.d{scanNum}.params = [];

  % display information about what we are doing
  dispHeader;
  disp(sprintf('(pRF) Scan %s:%i (restrict %s) running on %i processor(s)',params.groupName,scanNum,params.restrict,nProcessors));

  % get some information from pRFFit that will be used again in
  % the fits, including concatInfo, stim, prefit, etc.
  fit = pRFFit(v,scanNum,[],[],[],'fitTypeParams',params.pRFFit,'returnPrefit',true);
  if isempty(fit),return,end
  stim = fit.stim;
  pRFAnal.d{scanNum}.stim = cellArray(stim);
  pRFAnal.d{scanNum}.stimX = fit.stimX;
  pRFAnal.d{scanNum}.stimY = fit.stimY;
  pRFAnal.d{scanNum}.stimT = fit.stimT;
  concatInfo = fit.concatInfo;
  pRFAnal.d{scanNum}.concatInfo = fit.concatInfo;
  prefit = fit.prefit;
  paramsInfo = fit.paramsInfo;
  pRFAnal.d{scanNum}.paramsInfo = paramsInfo;
  % grab all these fields and stick them onto a structure called paramsInfo
  % preallocate some space
  rawParams = nan(fit.nParams,n);
  r = nan(n,fit.concatInfo.n);
  thisr2 = nan(1,n);
  thisPolarAngle = nan(1,n);
  thisEccentricity = nan(1,n);
  thisRfHalfWidth = nan(1,n);
  % compute pRF for each voxel in the restriction
  if params.pRFFit.prefitOnly,algorithm='prefit-only';else,algorithm=params.pRFFit.algorithm;end
  disp(sprintf('(pRF) Computing %s fits using %s for %i voxels',params.pRFFit.rfType,algorithm,n));
  % get some info about the scan to pass in (which prevents
  % pRFFit from calling viewGet - which is problematic for distributed computing
  framePeriod = viewGet(v,'framePeriod');
  junkFrames = viewGet(v,'junkFrames',scanNum);

  % now loop over each voxel
  parfor i = 1:n
    fit = pRFFit(v,scanNum,x(i),y(i),z(i),'stim',stim,'concatInfo',concatInfo,'prefit',prefit,'fitTypeParams',params.pRFFit,'dispIndex',i,'dispN',n,'tSeries',loadROI.tSeries(i,:)','framePeriod',framePeriod,'junkFrames',junkFrames,'paramsInfo',paramsInfo);
    if ~isempty(fit)
      % keep data, note that we are keeping temporarily in
      % a vector here so that parfor won't complain
      % then afterwords we put it into the actual overlay struct
      thisr2(i) = fit.r2;
      thisPolarAngle(i) = fit.polarAngle;
      thisEccentricity(i) = fit.eccentricity;
      thisRfHalfWidth(i) = fit.std;
      % keep parameters
      rawParams(:,i) = fit.params(:);
      r(i,:) = fit.r;
    end
  end
      
  % set overlays
  for i = 1:n
    r2.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);
    polarAngle.data{scanNum}(x(i),y(i),z(i)) = thisPolarAngle(i);
    eccentricity.data{scanNum}(x(i),y(i),z(i)) = thisEccentricity(i);
    rfHalfWidth.data{scanNum}(x(i),y(i),z(i)) = thisRfHalfWidth(i);
  end
  pRFAnal.d{scanNum}.params = rawParams;
  pRFAnal.d{scanNum}.r = r;

  iScan = find(params.scanNum == scanNum);
  thisParams.scanNum = params.scanNum(iScan);
  r2.params{scanNum} = thisParams;
  polarAngle.params{scanNum} = thisParams;
  eccentricity.params{scanNum} = thisParams;
  rfHalfWidth.params{scanNum} = thisParams; 
  % display how long it took
  disp(sprintf('(pRF) Fitting for %s:%i took in total: %s',params.groupName,scanNum,mlrDispElapsedTime(toc)));
end

% install analysis
pRFAnal.name = params.saveName;
pRFAnal.type = 'pRFAnal';
pRFAnal.groupName = params.groupName;
pRFAnal.function = 'pRF';
pRFAnal.reconcileFunction = 'defaultReconcileParams';
pRFAnal.mergeFunction = 'defaultMergeParams';
pRFAnal.guiFunction = 'pRFGUI';
pRFAnal.params = params;
pRFAnal.overlays = [r2 polarAngle eccentricity rfHalfWidth];
pRFAnal.curOverlay = 1;
pRFAnal.date = dateString;
v = viewSet(v,'newAnalysis',pRFAnal);
if ~isempty(viewGet(v,'fignum'))
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

% Save it
saveAnalysis(v,pRFAnal.name);

%set(viewGet(v,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    pRFAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(pRFAnal.d) == 1
    d = pRFAnal.d{1};
  else
    d = pRFAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getVoxelRestriction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y z] = getVoxelRestriction(v,params,scanNum)

x = [];y = [];z = [];

if strncmp(params.restrict,'Base: ',6)
  % get the base name
  baseName = params.restrict(7:end);
  baseNum = viewGet(v,'baseNum',baseName);
  if isempty(baseNum)
    disp(sprintf('(pRF) Could not find base to restrict to: %s',params.restrict));
    return;
  end
  base = viewGet(v,'base',baseNum);
  if isempty(base)
    disp(sprintf('(pRF) Could not find base to restrict to: %s',params.restrict));
    return;
  end
  if any(base.type == [1 2])
    % get base coordinates from the coordMap
    baseCoords = base.coordMap.coords;
    % convert to 4xn array
    sizeBaseCoords = size(baseCoords);
    baseCoords = reshape(baseCoords,prod(sizeBaseCoords(1:3)),3)';
    baseCoords(4,:) = 1;
    % and convert to scan coordinates
    base2scan = viewGet(v,'base2scan',scanNum,params.groupName,baseNum);
    scanCoords = round(base2scan*baseCoords);
    % check against scandims
    scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
    scanCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
    % remove duplicates and nans
    scanCoords = scanCoords(~isnan(scanCoords));
    scanCoords = unique(scanCoords);
    % convert back to x,y,z coordinates
    [x y z] = ind2sub(scanDims,scanCoords);
  else
    keyboard
  end
elseif strncmp(params.restrict,'ROI: ',5)
  % get the roi name
  roiName = params.restrict(6:end);
  scanCoords = getROICoordinates(v,roiName,scanNum,params.groupName,'straightXform=1');
  if isempty(scanCoords),return,end
  x = scanCoords(1,:);y = scanCoords(2,:);z = scanCoords(3,:);
elseif strncmp(params.restrict,'None',4)
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  [x y z]  = ndgrid(1:scanDims(1),1:scanDims(2),1:scanDims(3));
  x = x(:);y = y(:);z = z(:);
else
  keyboard
end