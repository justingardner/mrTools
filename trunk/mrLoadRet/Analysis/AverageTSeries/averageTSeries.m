function view = averageTSeries(view,params)
%
% function averageTSeries(view,[params])
%
% Computes average of the time series. Time reverses and shifts as
% specified. Output is a new tSeries file added to the specified group, and
% that is in the coordinate frame of the base scan.
%
% Order of operations is: 1) dump junk frames, 2) time shift, 3) time
% reverse, 4) average (transforming to coordinate frame of the base scan).
%
% params: Optional initial parameters (default: user is prompted via
%    averageTSeriesGUI). Params must be a structure with all of the
%    following fields.
% scanList: vector of scan numbers to include in the average.
%    Default: all of the scans.
% tseriesfiles: cell array of strings the same length as scanList,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
% shiftList: vector of integers the same length as scanList, each entry of
%    which specifies the time shift (in frames, i.e., 1 means a shift of 1
%    TR) for each scan. 
%    Default: all zero (no shift).
% reverseList: vector of the same length as scanList with non-zero entries
%    indicating which scans to time reverse.
%    Default: all zero (no time reversal).
% baseScan: number specifying scan which is used to specify the coordinate
%    frame for the result.
%    Default: first scan
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% groupName: group of scans that will be averaged.
%    Default: current group of view.
% aveGroupName: the new average scan is added to the specified group. If
%    the group specified by aveGroupName does not exist, then a new group
%    is created with this name.
%    Default: 'Averages'.
% description: description string for the new average scan
%    Default: 'Average from <groupName>' of scans: <scanList>'
% fileName: the new tseries is saved in the specified filename.
%    Default: 'tseries-mmddyy-hhmmss.img' or 'tseries-mmddyy-hhmmss.nii'
%    where mmddyy = date and hhmmss = time. Calls niftiFileExtension to
%    choose between .img and .nii, based on mrLoadRet 'niftiFileExtension'
%    preference.
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
%
% Examples:
%
% params = averageTSeriesGUI('groupName','Raw');
% view = newView('Volume');
% view = averageTSeries(view,params);
%
% view = averageTSeries(view);
%
%
% djh, 7/2006

% Get analysis parameters from averageTSeriesGUI.
nScans = viewGet(view,'nScans');
if ieNotDefined('params')
	% Initialize analysis parameters with default values
	params = averageTSeriesGUI('groupName',viewGet(view,'groupName'));
else
    % Reconcile params with current status of group and ensure that it has
    % the required fields. 
    params = averageTSeriesReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
    mrMsgBox('averageTSeries cancelled');
    return
end

% Retrieve parameters
scanList = params.scanList;
reverseList = params.reverseList;
shiftList = params.shiftList;
baseScan = params.baseScan;
groupName = params.groupName;
aveGroupName = params.aveGroupName;
description  = params.description;
interpMethod = params.interpMethod;

% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
    mrErrorDlg('averageTSeries: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the average group name. Create the
% group if necessary.
viewAverage = newView(viewGet(view,'viewType'));
aveGroupNum = viewGet(viewAverage,'groupNum',aveGroupName);
if isempty(aveGroupNum)
	view = viewSet(view,'newgroup',aveGroupName);
	aveGroupNum = viewGet(viewAverage,'groupNum',aveGroupName);
end
viewAverage = viewSet(viewAverage,'currentGroup',aveGroupNum);

% Check that all scans in scanList have the same nframes, frameperiod,
% scanxform, scanvoxelsize, scandims
nFrames = viewGet(viewBase,'nFrames',baseScan);
framePeriod = viewGet(viewBase,'framePeriod',baseScan);
voxelSize = viewGet(viewBase,'scanvoxelsize',baseScan);
scanDims = viewGet(viewBase,'scandims',baseScan);
scanXform = viewGet(viewBase,'scanxform',baseScan);
for iscan = 1:length(scanList)
	if (viewGet(viewBase,'nFrames',scanList(iscan)) ~= nFrames)
		mrErrorDlg('Can not average these scans because they have different numFrames.');
	end
	if (viewGet(viewBase,'framePeriod',scanList(iscan)) ~= framePeriod)
		mrWarnDlg('These scans  have different frame periods.');
	end
	if find(viewGet(viewBase,'scanvoxelsize',scanList(iscan)) ~= voxelSize)
		mrErrorDlg('Can not average these scans because they have different voxel sizes.');
	end
	if find(viewGet(viewBase,'scandims',scanList(iscan)) ~= scanDims)
		mrErrorDlg('Can not average these scans because they have different sizes.');
	end
end

% Compute output volume
aveTSeries = zeros([scanDims(1) scanDims(2) scanDims(3) nFrames]);
waitHandle = mrWaitBar(0,'Computing average tSeries.  Please wait...');
for iscan = 1:length(scanList)
	scanNum = scanList(iscan);
	reverse = reverseList(iscan);
	shift = shiftList(iscan);
	
	% Load it
	tseries = loadTSeries(viewBase,scanNum,'all');
	
	% Dump junk frames
	junkFrames = viewGet(viewBase,'junkframes',scanNum);
	tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
	
	% Time shift
	if (shift > 0)
		tseries = cat(4, tseries(:,:,:,[nFrames-shift+1:nFrames]), tseries(:,:,:,[1:nFrames-shift]));
	end
	if (shift < 0)
		shift = -shift;
		tseries = cat(4, tseries(:,:,:,[shift+1:nFrames]), tseries(:,:,:,[1:shift]));
	end

	% Time reverse
	if reverse
		tseries = flipdim(tseries,4);
	end
	
	% Compute transform
	% *** Not fully tested yet ***
    baseXform = viewGet(view,'scanXform',baseScan,groupNum);
    scanXform = viewGet(view,'scanXform',scanNum,groupNum);
    % Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the
    % origin.
    shiftXform = shiftOriginXform;
    M = inv(shiftXform) * inv(scanXform) * baseXform * shiftXform;
	
	% Warp the frames
    for frame = 1:nFrames
        tseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod);
    end   

	% Add 'em up
	aveTSeries = aveTSeries + tseries;
    % Update waitbar
    mrWaitBar(iscan/length(scanList),waitHandle);
end
% Divide by number of scans in scanList
aveTSeries = aveTSeries / length(scanList);
mrCloseDlg(waitHandle);

% Save aveTSeries (using header of 1st scan on scanList as the template for
% the nifti header), and add it as a new scan.
scanParams.fileName = [];
scanParams.junkFrames = 0;
scanParams.nFrames = nFrames;
scanParams.description = description;
hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',baseScan));
[viewAverage,tseriesFileName] = saveNewTSeries(viewAverage,aveTSeries,scanParams,hdr);

% Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = averageTSeries(view,params);'];
[pathstr,filename,ext,versn] = fileparts(tseriesFileName);
tseriesdir = viewGet(viewAverage,'tseriesdir');
save(fullfile(tseriesdir,filename),'evalstr','params','tseriesFileName');

% Delete temporary viewBase and viewAverage
deleteView(viewBase);
deleteView(viewAverage);

return; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** Testing/debugging *** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.scanList = 1;
params.shiftList = 2;
params.reverseList = 0;
params.baseScan = 1;
params.groupName = 'Raw';
params.aveGroupName = 'Averages';
params.fileName = [];
params.description = 'test';

MLR.views{1} = averageTSeries(MLR.views{1},params);

MLR.views{1} = viewSet(MLR.views{1},'currentGroup',1);
ts1 = loadTSeries(MLR.views{1},1,6);
MLR.views{1} = viewSet(MLR.views{1},'currentGroup',2);
ts2 = loadTSeries(MLR.views{1},1,6);
figure; plot([ts1(9:168,2000) ts2(:,2000)]);



