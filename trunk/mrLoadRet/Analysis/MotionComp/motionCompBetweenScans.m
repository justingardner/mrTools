function view = motionCompBetweenScans(view,params)
%
% view = motionCompBetweenScans(view,[params])
%
% Robust 3D rigid body motion compensation between different scans.
% Computes the mean over time of each scan, and then computes the
% registration between those mean volumes. Finally, rather than
% interpolating the time series, this function simply copies the tseries
% file, overwriting the sform of the nifti header.
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure with all of the
%    following fields.
% baseScan: targetScans are registered to the baseScan
%    Default scan 1.
% targetScans: which scans to apply motion compensation (normally including the
%    baseScan).
%    Default all scans
% robust: use robust M-estimator for motion estimates.
%    Default 0.
% correctIntensityContrast: intensity and contrast normalization before
%    registration.
%    Default 0.
% crop specifies border size to crop/ignore around all sides of the volume.
%    Should be of the form [ymin xmin zmin; ymax xmax zmax]
%    Default [].
% niters: number of iterations in the motion estimation.
%    Default 3.
% tseriesfiles: cell array of strings the same length as targetScans,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
%
% If you change this function make parallel changes in:
%    motionComp, motionCompWithinScan
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
% view = newView;
% view = motionCompBetweenScans(view,params);
%
% view = motionCompBetweenScans(view);
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams

% Get analysis parameters from motionCompGUI.
nScans = viewGet(view,'nScans');
if (nScans == 0)
    mrWarnDlg('(motionCompBetweenScans) No scans in group');
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
targetScans = params.targetScans;
robust = params.robust;
correctIntensityContrast = params.correctIntensityContrast;
crop = params.crop;
niters = params.niters;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;
interpMethod = params.interpMethod;

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

% Get base scanParams
scanParams = viewGet(viewBase,'scanParams',baseScan);
%baseSform = scanParams.niftiHdr.sform44;

% Open wait bar
waitHandle = mrWaitBar(0,'Computing between scan motion compensation.  Please wait...');

% Compute baseMean
% Load tseries and dump junk frames
tseries = loadTSeries(viewBase,baseScan,'all');
junkFrames = viewGet(viewBase,'junkframes',baseScan);
nFrames = viewGet(viewBase,'nFrames',baseScan);
% Compute mean of base scan
baseMean = nanmean(tseries(:,:,:,junkFrames+1:junkFrames+nFrames),4);
% Intensity and contrast correction
if correctIntensityContrast
	baseMean = intensityContrastCorrection(baseMean,crop);
end
	
% Loop through target scans, compute motion estimate transform, warp
% according to the transform, and save new tseries.
for s = 1:length(targetScans)
  scanNum = targetScans(s);
	
  % Update wait bar
  mrWaitBar(s/length(targetScans),waitHandle);
    
  % Load tseries 
  tseries = loadTSeries(viewBase,scanNum,'all');
	
  if (scanNum == baseScan)
    % Do not bother computing registration for baseScan (just copy it)
    transform = eye(4);
  else
    % Compute mean of target scan after dumping junk frames
    junkFrames = viewGet(viewBase,'junkframes',scanNum);
    nFrames = viewGet(viewBase,'nFrames',scanNum);
    tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
    targetMean = nanmean(tseries,4);
    % Intensity and contrast correction
    if correctIntensityContrast
      targetMean = intensityContrastCorrection(targetMean,crop);
    end
    % Estimate motion between mean volumes
    disp(['Computing motion estimates for scan ',num2str(scanNum),'...'])
    transform = estMotionIter3(baseMean,targetMean,niters,eye(4),1,robust,0,crop);
    disp(['Computing motion estimates for scan ',num2str(scanNum),'... done'])
    
    % jg: warp/transform each frame of the time series, rather than saving
    % the transformation in the sform
    if ~isequal(transform,eye(4))
      disp(['Warping images for scan ',num2str(scanNum)])
      for frame = 1:nFrames
	tseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),transform,NaN,0,interpMethod);
      end
      disp(['Warping images for scan ',num2str(scanNum),'... done'])
    end
  end
    
  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.description = ['Between ' descriptions{s}];
  scanParams.junkFrames = 0;
  scanParams.nFrames = nFrames;
  scanParams.fileName = [];
  %jg: removed the saving in the sform bit here, since we now
  %compute the transformation, also reset the junk frames and
  %nFrames fields so now we have new data
  %sform = baseSform * inv(transform);
  %scanParams.niftiHdr = cbiSetNiftiSform(scanParams.niftiHdr,sform);
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,tseries,scanParams,scanParams.niftiHdr);
  
  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = motionCompBetweenScans(view,params);'];
  [pathstr,filename,ext,versn] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transform','tseriesFileName');
  
end 
% Close waitbar
mrCloseDlg(waitHandle);

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
