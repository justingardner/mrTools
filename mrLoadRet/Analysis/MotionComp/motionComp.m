function [view  params] = motionComp(view,params,varargin)
%
% view = motionComp(view,[params])
%
% Robust 3D rigid body motion compensation. First, performs within scan
% motion compensation on the base scan. Then, computes the mean (over time)
% of the motion-compensated base scan. Finally, motion corrects each
% individual frame of the target scans to the mean motion-corrected base
% scan.
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure with all of the
%    following fields.
% targetScans: which scans to apply motion compensation.
%    Default all scans.
% baseScan: targetScans are registered to the motion-compensated baseScan.
%    Default scan 1.
% baseFrame: string ('first', last', or 'mean') specifying frame in
%    baseScan to which the rest of the baseScan is registered (ignoring
%    junk frames).
%    Default 'first'.
% robustFlag: use robust M-estimator for motion estimates.
%    Default 0.
% correctIntensityContrast: intensity and contrast normalization before
%    registration.
%    Default 0.
% crop specifies border size to crop/ignore around all sides of the volume.
%    Should be of the form [ymin xmin zmin; ymax xmax zmax]
%    Default [].
% niters: number of iterations in the motion estimation.
%    Default 3.
% sliceTimeCorrection: Performs slice time correction if True.
% sliceTimeString: Specifies slice time correction:
%    'beginning of TR'
%    'middle of TR'
%    'end of TR' (default)
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% tseriesfiles: cell array of strings the same length as targetScans,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
%
% If you change this function make parallel changes in:
%   MotionCompBetweenScans, motionCompWithScan
%
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
%
% Examples:
%
% params = motionCompGUImrParams('groupName','Raw');
% view = newView;
% view = motionComp(view,params);
%
% view = motionComp(view);
%
% You can also just get a parameters structure in the following
% ways:
%
%  v = newView;
%  [v params] = motionComp(v,[],'justGetParams=1');
%  [v params] = motionComp(v,[],'justGetParams=1','defaultParams=1');
%  [v params] = motionComp(v,[],'justGetParams=1','defaultParams=1','scanList=[1 2]');
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams
% Get analysis parameters from motionCompGUI.
% clear time series on each iteration to avoid running out of
% memory
%
% $Id$	
%
nScans = viewGet(view,'nScans');

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end

if (nScans == 0)
  mrWarnDlg('(motionComp) No scans in group');
  return
end

if ieNotDefined('params')
  % Initialize analysis parameters with default values
  %    params = motionCompGUI('groupName',viewGet(view,'groupName'));
  params = motionCompGUImrParams('groupName',viewGet(view,'groupName'),'defaultParams',defaultParams,'scanList',scanList);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields.
  params = motionCompReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('motion compensation cancelled');
  return
end

% if just getting params then return
if justGetParams,return,end
 
% Retrieve parameters
baseScan = params.baseScan;
baseFrame = params.baseFrame;
targetScans = params.targetScans;
sliceTimeCorrection = params.sliceTimeCorrection;
sliceTimeString = params.sliceTimeString;
robust = params.robust;
correctIntensityContrast = params.correctIntensityContrast;
crop = params.crop;
niters = params.niters;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
interpMethod = params.interpMethod;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;
tSmooth = params.tSmooth;

% temporal smoothing option.
if tSmooth ~= 0
  tSmooth = 2*fix(tSmooth/2) + 1;
else
  tSmooth = 1;
end

% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('motionComp: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Error if the number of slices is too small (the derivative computation
% discards the borders in z, 2 slices at the begining and 2 more at the
% end).
baseDataSize = viewGet(viewBase,'datasize',baseScan);
if baseDataSize(3) < 8
  mrErrorDlg('Motion compensation requires at least 8 slices');
end

% Ignore scans if data size is different from that for base scan.
for scanNum = targetScans
  datasize = viewGet(viewBase,'datasize',scanNum);
  if (datasize ~= baseDataSize)
    mrWarnDlg(['Ignoring scan ',num2str(scanNum),'. Motion compensation requires at the datasize to match the base scan']);
    targetScans = targetScans(find(targetScans ~= scanNum));
  end
end

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
viewMotionComp = newView(viewGet(view,'viewType'));
motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
if isempty(motionCompGroupNum)
  view = viewSet(view,'newgroup',motionCompGroupName);
  motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
end
viewMotionComp = viewSet(viewMotionComp,'currentGroup',motionCompGroupNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within scan motion compensation on base scan %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load tseries 
scanNum = baseScan;
tseries = loadTSeries(viewBase,scanNum,'all');
junkFrames = viewGet(viewBase,'junkframes',scanNum);
nFrames = viewGet(viewBase,'nFrames',scanNum);
totalFrames = viewGet(viewBase,'totalFrames',scanNum);

% Initialize the warped time series to zeros.
warpedTseries = zeros(size(tseries));

% Get slice times and replicate the last frame of tseries for slice time
% correction 
if sliceTimeCorrection
  sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
  tseries(:,:,:,end+1) = tseries(:,:,:,end);
  switch sliceTimeString
    case 'end of TR'
      sliceTimes = sliceTimes;
    case 'middle of TR'
      sliceTimes = sliceTimes - 0.5;
    case 'beginning of TR'
      sliceTimes = sliceTimes - 1;
    otherwise
      mrErrorDlg('Invalid slice times');
  end
else
  sliceTimes = [];
end

% Intensity/contrast correction
if correctIntensityContrast
  tseriesIC = intensityContrastCorrection(tseries,crop);
else
  tseriesIC = tseries;
end

% Get volume corresponding to base frame. Other frames will be motion
% compensated to this one.
switch baseFrame
  case 'first'
    baseF = junkFrames+1;
    baseVol = tseriesIC(:,:,:,baseF);
  case 'last'
    baseF = size(tseriesIC,4);
    baseVol = tseriesIC(:,:,:,baseF);
  case 'mean'
    baseF = 0;
    baseVol = nanmean(tseriesIC(:,:,:,junkFrames+1:junkFrames+nFrames),4);
  otherwise
    mrErrorDlg('Invalid base frame');
end

% Loop: computing motion estimates and warping the volumes to
% compensate for the motion in each temporal frame.
waitHandle = mrWaitBar(0,['Computing within scan motion compensation for base scan ',num2str(scanNum),'.  Please wait...']);
M = eye(4);
for frame = 1:totalFrames
  mrWaitBar(frame/totalFrames,waitHandle)
  if (frame == baseF)
    M = eye(4);
  else
    if (frame <= junkFrames)
      Minitial = eye(4);
    else
      Minitial = M;
    end
    % Compute rigid-body motion estimate
    if sliceTimeCorrection
      if strcmp(baseFrame,'mean')
        M = estMotionInterp3(baseVol,tseriesIC,1,frame,niters,Minitial,sliceTimes,1,robust,0,crop);
      else
        M = estMotionInterp3(tseriesIC,tseriesIC,baseF,frame,niters,Minitial,sliceTimes,1,robust,0,crop);
      end
    else
      frameMin = frame - fix(tSmooth/2);
      if frameMin < 1, frameMin = 1; end
      frameMax = frame + fix(tSmooth/2);
      if frameMax > totalFrames, frameMax = totalFrames; end
      vol = nanmean(tseriesIC(:,:,:,frameMin:frameMax), 4);
      M = estMotionIter3(baseVol,vol,niters,Minitial,1,robust,0,crop);
    end
  end
  % Warp the volume
  if sliceTimeCorrection
    warpedTseries(:,:,:,frame) = warpAffineInterp3(tseries,frame,M,sliceTimes,NaN,interpMethod);
  else
    warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod);
  end
end
mrCloseDlg(waitHandle);
% need to clear to make room for large datasets
clear tseries tseriesIC

% Finally, compute mean over time (ignoring junkFrames)
baseMean = nanmean(warpedTseries(:,:,:,junkFrames+1:junkFrames+nFrames),4);

if correctIntensityContrast
  baseMean = intensityContrastCorrection(baseMean,crop);
end

% clear temporary tseries
clear warpedTseries;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through target scans, perform motion estimation for each frame, %
% warp according to motion estimates and save new tseries.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(targetScans)
  scanNum = targetScans(s);
  
  % Load tseries
  tseries = loadTSeries(viewBase,scanNum,'all');
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  totalFrames = viewGet(viewBase,'totalFrames',scanNum);
  
  % Initialize the warped time series to zeros. Need to re-initialize this
  % for each scan because number of frames can differ. 
  warpedTseries = zeros(size(tseries));
  
  % Get slice times and replicate the last frame of tseries for slice time
  % correction
  if sliceTimeCorrection
    sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
    tseries(:,:,:,end+1) = tseries(:,:,:,end);
    switch sliceTimeString
      case 'end of TR'
        sliceTimes = sliceTimes;
      case 'middle of TR'
        sliceTimes = sliceTimes - 0.5;
      case 'beginning of TR'
        sliceTimes = sliceTimes - 1;
      otherwise
        mrErrorDlg('Invalid slice times');
    end
  else
    sliceTimes = [];
  end
  
  % Intensity/contrast correction
  if correctIntensityContrast
    tseriesIC = intensityContrastCorrection(tseries,crop);
  else
    tseriesIC = tseries;
  end
  
  % Loop through frames of target scan and estimate motion params
  waitHandle = mrWaitBar(0,['Computing motion estimates for scan ',num2str(scanNum),'.  Please wait...']);
  transforms = cell(1,totalFrames);
  M = eye(4);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    if (frame <= junkFrames)
      Minitial = eye(4);
    else
      Minitial = M;
    end
    % Compute rigid-body motion estimate with respect to baseMean
    if sliceTimeCorrection
      M = estMotionInterp3(baseVol,tseriesIC,1,frame,niters,Minitial,sliceTimes,1,robust,0,crop);
    else
      frameMin = frame - fix(tSmooth/2);
      if frameMin < 1, frameMin = 1; end
      frameMax = frame + fix(tSmooth/2);
      if frameMax > totalFrames, frameMax = totalFrames; end
      vol = nanmean(tseriesIC(:,:,:,frameMin:frameMax), 4);
      M = estMotionIter3(baseVol,vol,niters,Minitial,1,robust,0,crop);
    end
    % Collect the transform
    transforms{frame} = M;
  end
  clear tseriesIC;
  mrCloseDlg(waitHandle);
  
  % warp the images according to the motion estimates
  waitHandle = mrWaitBar(0,['Warping image volumes for scan ',num2str(scanNum),'.  Please wait...']);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    if sliceTimeCorrection
      warpedTseries(:,:,:,frame) = warpAffineInterp3(tseries,frame,transforms{frame},sliceTimes,NaN,interpMethod);
    else
      warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),transforms{frame},NaN,0,interpMethod);
    end
  end
  mrCloseDlg(waitHandle);
  clear tseries;
  
  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.description = ['Full ' descriptions{s}];
  scanParams.fileName = [];
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);

  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = motionComp(view,params);'];
  [pathstr,filename,ext,versn] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');

  % clear temporary tseries
  clear warpedTseries;

end


% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
