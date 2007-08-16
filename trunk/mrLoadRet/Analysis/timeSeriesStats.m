function view = timeSeriesStats(view,params)
%
% view = timeSeriesStats(view,[params])
% 
% Loops throughs scans, loads corresponding tSeries, and computes time
% series statistics for each voxel:
% - mean 
% - median 
% - standard deviation
% - max frame-to-frame difference
% - max difference from median
%
% params: Optional initial parameters. Default: user is prompted via
%    GUI. Params must be a a structure with all of the following fields. 
% params.groupName: group of scans that will be analyzed.
%    Default: current group of view.
% params.scanList: vector specifying which scans to compute.
%    Default: all of the scans.
%
%
% Examples:
%
% paramss.groupName = 'Raw';
% n = viewGet([],'nScans',1)
% params.scanList = [1:n];
% view = timeSeriesStats(view,params);
%
% view = timeSeriesStats(view);
%
%
% djh, 7/2007
% $Id$	

if ~isview(view)
    help timeSeriesStats
    mrErrorDlg('(timeSeriesStats) Invalid view.')
end

% Get analysis parameters from timeSeriesStatsGUI
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  groupName = viewGet(view,'groupName');
  groupNum = viewGet(view,'groupNum',groupName);
  n = viewGet(view,'nScans',groupNum);
  params.groupName = groupName;
  params.scanList = selectScans(view);%[1:n];
end

% Reconcile params with current status of group and ensure that params
% has the required fields.
params = defaultReconcileParams(params.groupName,params);

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('timeSeriesStats cancelled',1);
  return
end

% Change group
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end

% Compute it
[tsMean,tsMedian,tsStd,tsMaxFrameDiff,tsMaxMedianDiff] = computeTimeSeriesStats(view,params);

% Make analysis structure
tsStats.name = 'timeSeriesStats';  % This can be reset by editAnalysisGUI
tsStats.type = 'timeSeriesStats';
tsStats.groupName = params.groupName;
tsStats.function = 'timeSeriesStats';
tsStats.guiFunction = 'timeSeriesStatsGUI';
tsStats.params = params;

% Install it in the view
view = viewSet(view,'newanalysis',tsStats);
view = viewSet(view,'newoverlay',tsMean);
view = viewSet(view,'newoverlay',tsMedian);
view = viewSet(view,'newoverlay',tsStd);
view = viewSet(view,'newoverlay',tsMaxFrameDiff);
view = viewSet(view,'newoverlay',tsMaxMedianDiff);

% Save it
saveAnalysis(view,tsStats.name);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tsMean,tsMedian,tsStd,tsMaxFrameDiff,tsMaxMedianDiff] = ...
  computeTimeSeriesStats(view,params)

% Get nScans from view and get scanList from params
scanList = params.scanList;
nScans = viewGet(view,'nScans');

% Intialize overlay structures

% mean
tsMean.name = 'mean';
tsMean.function = 'timeSeriesStats';
tsMean.data = cell(1,nScans);
tsMean.params = params;
tsMean.colormap = jet(256);
tsMean.groupName = params.groupName;

% median
tsMedian.name = 'median';
tsMedian.function = 'timeSeriesStats';
tsMedian.data = cell(1,nScans);
tsMedian.params = params;
tsMedian.colormap = jet(256);
tsMedian.groupName = params.groupName;

% std
tsStd.name = 'std';
tsStd.function = 'timeSeriesStats';
tsStd.data = cell(1,nScans);
tsStd.params = params;
tsStd.colormap = jet(256);
tsStd.groupName = params.groupName;

% max frame-to-frame diff
tsMaxFrameDiff.name = 'maxFrameDiff';
tsMaxFrameDiff.function = 'timeSeriesStats';
tsMaxFrameDiff.data = cell(1,nScans);
tsMaxFrameDiff.params = params;
tsMaxFrameDiff.colormap = jet(256);
tsMaxFrameDiff.groupName = params.groupName;

% max diff from median
tsMaxMedianDiff.name = 'maxMedianDiff';
tsMaxMedianDiff.function = 'timeSeriesStats';
tsMaxMedianDiff.data = cell(1,nScans);
tsMaxMedianDiff.params = params;
tsMaxMedianDiff.colormap = jet(256);
tsMaxMedianDiff.groupName = params.groupName;

disp('Computing time series statistics...');
waitHandle = mrWaitBar(0,'Computing statistics from the tSeries.  Please wait...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
    scanNum = scanList(scanIndex);
    disp(['Processing scan ', int2str(scanNum),'...']);
    
    % sliceDims: [ydim xdim] for single slice
    % volDims; [ydim xdim nslices] for single scan
    sliceDims = viewGet(view,'sliceDims',scanNum);
    volDims = viewGet(view,'dims',scanNum);
    
    % Initialize data with NaNs
    tsMean.data{scanNum} = NaN*ones(volDims);
    tsMedian.data{scanNum} = NaN*ones(volDims);
    tsStd.data{scanNum} = NaN*ones(volDims);
    tsMaxFrameDiff.data{scanNum} = NaN*ones(volDims);
    tsMaxMedianDiff.data{scanNum} = NaN*ones(volDims);
    
    nslices = viewGet(view,'nslices',scanNum);
    for sliceNum = 1:nslices
        [tsMeanSeries,tsMedianSeries,tsStdSeries,tsMaxFrameDiffSeries,tsMaxMedianDiffSeries] = ...
          computeTimeSeriesStatsSeries(view,scanNum,sliceNum);
        tsMean.data{scanNum}(:,:,sliceNum) = reshape(tsMeanSeries,sliceDims);
        tsMedian.data{scanNum}(:,:,sliceNum) = reshape(tsMedianSeries,sliceDims);
        tsStd.data{scanNum}(:,:,sliceNum) = reshape(tsStdSeries,sliceDims);
        tsMaxFrameDiff.data{scanNum}(:,:,sliceNum) = reshape(tsMaxFrameDiffSeries,sliceDims);
        tsMaxMedianDiff.data{scanNum}(:,:,sliceNum) = reshape(tsMaxMedianDiffSeries,sliceDims);
    end
    % Update waitbar
    mrWaitBar(scanIndex/length(scanList),waitHandle);
end

mrCloseDlg(waitHandle);

% Fill range fields 
tsMean.range = findRange(tsMean.data);
tsMedian.range = findRange(tsMean.data);
tsStd.range = findRange(tsMean.data);
tsMaxFrameDiff.range = findRange(tsMean.data);
tsMaxMedianDiff.range = findRange(tsMean.data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tsMeanSeries,tsMedianSeries,tsStdSeries,tsMaxFrameDiffSeries,tsMaxMedianDiffSeries] = ...
  computeTimeSeriesStatsSeries(view,scan,slice)

% Get junk frames and nframes
junkframes = viewGet(view,'junkframes',scan);
nframes = viewGet(view,'nframes',scan);

% Load tSeries
tSeries = loadTSeries(view, scan, slice);

% Reshape the tSeries
tSeries = reshapeTSeries(tSeries);

% Remove junkFrames
tSeries = tSeries(junkframes+1:junkframes+nframes,:);

tsMeanSeries = mean(tSeries);
tsMedianSeries = median(tSeries);
tsStdSeries = std(tSeries);
tsMaxFrameDiffSeries = max(abs(tSeries(2:end,:)-tSeries(1:end-1,:)));
tsMaxMedianDiffSeries = max(abs(tSeries - repmat(tsMedianSeries,[nframes,1])));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
    ampMin = min([ampMin min(data{scan}(:))]);
    ampMax = max([ampMax max(data{scan}(:))]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end

