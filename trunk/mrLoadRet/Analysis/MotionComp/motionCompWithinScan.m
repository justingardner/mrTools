function view = motionCompWithinScan(view,params)
%
% view = motionCompWthinScan(view,[params])
%
% Robust 3D rigid body motion compensation
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure with all of the
%    following fields.
% targetScans: which scans to apply motion compensation.
%    Default all scans.
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
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% tseriesfiles: cell array of strings the same length as targetScans,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
%
% If you change this function make parallel changes in:
%   MotionCompBetweenScans, motionComp
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
% params = motionCompGUI('groupName','Raw');
% view = newView('Volume');
% view = motionCompWithinScan(view,params);
%
% view = motionCompWithinScan(view);
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams,
% clear time series on each iteration to avoid running out of memory

% Get analysis parameters from motionCompGUI.
nScans = viewGet(view,'nScans');
if (nScans == 0)
  mrWarnDlg('(motionCompWithinScans) No scans in group');
  return
end

if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = motionCompGUI('groupName',viewGet(view,'groupName'));
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

% Retrieve parameters
baseFrame = params.baseFrame;
targetScans = params.targetScans;
sliceTimeCorrection = params.sliceTimeCorrection;
robust = params.robust;
correctIntensityContrast = params.correctIntensityContrast;
crop = params.crop;
niters = params.niters;
interpMethod = params.interpMethod;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;

% *** Temporary warning until it is implemented and tested
if sliceTimeCorrection
  mrWarnDlg('Slice time correction not yet implemented');
end

% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('motionComp: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Ignore scans if the number of slices is too small (the derivative
% computation discards the borders in z, 2 slices at the begining and 2
% more at the end).
for scanNum = targetScans
  datasize = viewGet(viewBase,'datasize',scanNum);
  if (datasize(3) < 8)
    mrWarnDlg(['Ignoring scan ',num2str(scanNum),'. Motion compensation requires at least 8 slices']);
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

% Loop through target scans, perform motion estimation for each frame, warp
% according to motion estimates and save new tseries.
for s = 1:length(targetScans)
  scanNum = targetScans(s);
  % clear old tseries
  clear tseries warpedTseries;

  % Load tseries and dump junk frames
  tseries = loadTSeries(viewBase,scanNum,'all');
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
  
  % Intensity/contrast correction
  if correctIntensityContrast
    tseriesIC = intensityContrastCorrection(tseries,crop);
  else
    tseriesIC = tseries;
  end
  
  % Slice order for slice time correction
  if sliceTimeCorrection
    sliceOrder = viewGet(viewBase,'sliceOrder',scanNum);
  else
    sliceOrder = [];
  end

  % Get volume corresponding to base frame.  Other frames will be motion
  % compensated to this one.
  switch baseFrame
    case 'first'
      baseF = 1;
      baseVol = tseriesIC(:,:,:,baseF);
    case 'last'
      baseF = size(tseriesIC,4);
      baseVol = tseriesIC(:,:,:,baseF);
    case 'mean'
      baseF = 0;
      baseVol = nanmean(tseriesIC,4);
    otherwise
      mrErrorDlg('Invalid base frame');
  end

  % Initialize the warped time series to zeros.
  warpedTseries = zeros(size(tseries));

  % Loop: computing motion estimates and warping the volumes to
  % compensate for the motion in each temporal frame.
  waitHandle = mrWaitBar(0,['Computing motion compensation for scan ',num2str(scanNum),'.  Please wait...']);
  transforms = cell(1,nFrames);
  M = eye(4);
  for frame = 1:nFrames
    mrWaitBar(frame/nFrames,waitHandle)
    if (frame == baseF)
      M = eye(4);
    else
      % Compute rigid-body motion estimate
      vol = tseriesIC(:,:,:,frame);
      M = estMotionIter3(baseVol,vol,niters,M,1,robust,crop);
    end
    % Collect the transform
    transforms{frame} = M;
    % Warp the volume
    warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod);
  end
  mrCloseDlg(waitHandle);

  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.junkFrames = 0;
  scanParams.nFrames = nFrames;
  scanParams.description = ['Within ' descriptions{s}];
  scanParams.fileName = [];
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);

  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = motionCompWithinScan(view,params);'];
  [pathstr,filename,ext,versn] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');

end

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
