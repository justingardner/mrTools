function view = motionComp(view,params)
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
% params = motionCompGUI('groupName','Raw');
% view = newView('Volume');
% view = motionComp(view,params);
%
% view = motionComp(view);
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams
% Get analysis parameters from motionCompGUI.
% clear time series on each iteration to avoid running out of memory
nScans = viewGet(view,'nScans');

if (nScans == 0)
    mrWarnDlg('(motionComp) No scans in group');
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
baseScan = params.baseScan;
baseFrame = params.baseFrame;
targetScans = params.targetScans;
robust = params.robust;
correctIntensityContrast = params.correctIntensityContrast;
crop = params.crop;
niters = params.niters;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
interpMethod = params.interpMethod;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;

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

% Load tseries and dump junk frames
scanNum = baseScan;
tseries = loadTSeries(viewBase,scanNum,'all');
junkFrames = viewGet(viewBase,'junkframes',scanNum);
nFrames = viewGet(viewBase,'nFrames',scanNum);
tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
if correctIntensityContrast
    tseriesIC = intensityContrastCorrection(tseries,crop);
else
	tseriesIC = tseries;
end
	
% Get volume corresponding to base frame.  Other frames will be motion
% compensated to this one.
switch baseFrame
    case 'first'
		baseFrame = 1;
        baseVol = tseriesIC(:,:,:,baseFrame);
    case 'last'
		baseFrame = size(tseriesIC,4);
        baseVol = tseriesIC(:,:,:,baseFrame);
    case 'mean'
		baseFrame = 0;
        baseVol = nanmean(tseriesIC,4);
    otherwise
        mrErrorDlg('Invalid base frame');
end
	
% Initialize the warped time series to zeros.
warpedTseries = zeros(size(tseries));
	
% Loop: computing motion estimates and warping the volumes to
% compensate for the motion in each temporal frame.
waitHandle = mrWaitBar(0,['Computing within scan motion compensation for base scan ',num2str(scanNum),'.  Please wait...']);
M = eye(4);
for frame = 1:nFrames
    mrWaitBar(frame/nFrames,waitHandle)
    if (frame == baseFrame)
		M = eye(4);
	else
		% Compute rigid-body motion estimate
        vol = tseriesIC(:,:,:,frame);
        M = estMotionIter3(baseVol,vol,niters,M,1,robust,crop);
	end
	% Warp the volume
	warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod);
end
mrCloseDlg(waitHandle);

% Finally, compute mean over time
baseMean = nanmean(warpedTseries,4);
if correctIntensityContrast
    baseMean = intensityContrastCorrection(baseMean,crop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through target scans, perform motion estimation for each frame, %
% warp according to motion estimates and save new tseries.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(targetScans)
    scanNum = targetScans(s);
    % clear old tseries 
    clear tseries warpedTseries;    
	% Load tseries and dump junk frames
	tseries = loadTSeries(viewBase,scanNum,'all');
	junkFrames = viewGet(viewBase,'junkframes',scanNum);
	nFrames = viewGet(viewBase,'nFrames',scanNum);
	tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
    if correctIntensityContrast
        tseriesIC = intensityContrastCorrection(tseries,crop);
	else
		tseriesIC = tseries;
    end
		
	% Initialize the warped time series to zeros and fill the base frame.
	warpedTseries = zeros(size(tseries));
    
    % Loop through frames of target scan
	waitHandle = mrWaitBar(0,['Computing motion compensation for scan ',num2str(scanNum),'.  Please wait...']);
	transforms = cell(1,nFrames);
	M = eye(4);
	for frame = 1:nFrames
		mrWaitBar(frame/nFrames,waitHandle)
        vol = tseriesIC(:,:,:,frame);
        % Compute rigid-body motion estimate with respect to baseMean
        M = estMotionIter3(baseMean,vol,niters,M,1,robust,crop);
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
	scanParams.description = ['Full ' descriptions{s}];
	scanParams.fileName = [];
	scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
	scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
	try
	  % remember current verbose level, but set it temporarily to 0
	  currentVerbose = viewGet([],'pref','verbose');
	  if isempty(currentVerbose)
	    currentVerbose = 0;
	  end
	  viewSet(view,'pref','verbose',0);
	  % save the time series
	  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);
	catch
	  % reset verbose preference
	  viewSet(view,'pref','verbose',currentVerbose);
	  % display error and rethrow it
	  err = lasterror;
	  if isfield(err,'stack'),disp(sprintf('Error in ==> %s at %i',err.stack(1).name,err.stack(1).line)); end
	  rethrow(err);
	end
	% reset verbose preference
	viewSet([],'pref','verbose',currentVerbose);
	% Save evalstring for recomputing and params
	evalstr = ['view = newView(','''','Volume','''','); view = motionComp(view,params);'];
	[pathstr,filename,ext,versn] = fileparts(tseriesFileName);
	tseriesdir = viewGet(viewMotionComp,'tseriesdir');
	save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');

end 

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
