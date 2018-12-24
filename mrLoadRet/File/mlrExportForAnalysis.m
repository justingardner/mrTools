% mlrExportForAnalysis.m
%
%        $Id:$ 
%      usage: output = mlrExportForAnalysis(v)
%         by: justin gardner
%       date: 12/06/18
%    purpose: Export a structure for data analysis by pyton etc
%
function output = mlrExportForAnalysis(varargin)

% check arguments
if ~any(nargin == [2])
  help mlrExportForAnalysis
  return
end

% check if first argument is not view, that
% means that we have been called from the gui
% so we got an hObject
if ishandle(varargin{1})
  % get view
  v = viewGet(getfield(guidata(varargin{1}),'viewNum'),'view');
elseif isview(varargin{1})
  v = varargin{1};
else
  disp(sprintf('(mlrExportForAnalysis) Call with a view'));
  return
end

% start to create paramsInfo for choosing what to save
paramsInfo = {{'fullTimeSeries',0,'type=checkbox','Save full time series'}};

% get what ROIs are loaded
if viewGet(v,'nROIs') > 0
  paramsInfo{end+1} = {'roi',1,'type=checkbox','Save roi'};
  paramsInfo{end+1} = {'roiFullTimeSeries',1,'type=checkbox','contingent=roi','Save full time series for all the voxels in the rois. If unclicked, will still the average roi time series'};
else
  paramsInfo{end+1} = {'roi',0,'type=checkbox','visible=0'};
  paramsInfo{end+1} = {'roiTimeSeries',0,'type=checkbox','visible=0'};
end

% check if there is an analysis
analysis = viewGet(v,'analysis');

% check for event-related
if ~isempty(analysis) &&  strcmp(analysis.type,'erAnal')
  paramsInfo{end+1} = {'eventRelatedAnalysis',1,'type=checkbox','Save event-related analysis'};
else
  paramsInfo{end+1} = {'eventRelatedAnalysis',0,'visible=0'};
end
  
% put up dialog for choice on what to save
params = mrParamsDialog(paramsInfo,'Choose export options');
if isempty(params),return,end

% get save name
filename = sprintf('%s_%s_%i.mat',getLastDir(viewGet(v,'homeDir')),viewGet(v,'groupName'),viewGet(v,'curScan'));
[filename,pathname] = uiputfile({'*.mat','Matlab file'},'Save as',filename);
if isequal(filename,0) || isequal(pathname,0)
  return
end

% save experiment info
output.experimentName = getLastDir(viewGet(v,'homeDir'));
output.groupName = viewGet(v,'groupName');
output.scanNum = viewGet(v,'curScan');
output.description = viewGet(v,'description');

% load the time series
if params.fullTimeSeries
  disppercent(-inf,'(mlrExportForAnalysis) Loading time series');
  output.tSeries = loadTSeries(v);
  disppercent(inf);
else
  output.tSeries = [];
end

% set dimensions
output.dims = viewGet(v,'scanDims');
output.dims(end+1) = viewGet(v,'totalFrames');
output.framePeriod = viewGet(v,'framePeriod');
output.concatInfo = viewGet(v,'concatInfo');

% load the rois
if params.roi
  % get how many rois
  output.nroi = viewGet(v,'nrois');
  % and read in each one
  for iROI = 1:output.nroi
    % load the roi
    output.roi(iROI) = loadROITSeries(v,iROI,output.scanNum,output.groupName);
    % compute mean of time series is we are not saving full
    if ~params.roiFullTimeSeries
      output.roi(iROI).tSeries = mean(output.roi(iROI).tSeries);
    end
  end
else
  output.nroi = 0;
  output.roi = [];
end

% save event related analysis if it is there
if params.eventRelatedAnalysis
  output.eventRelated = analysis.d{output.scanNum};
end

% load overlays
nOverlays = viewGet(v,'nOverlays');
for iOverlay = 1:nOverlays
  % get the overlay
  output.overlay(iOverlay) = viewGet(v,'overlay',iOverlay);
  % extract only the data for this scan
  output.overlay(iOverlay).data = viewGet(v,'overlayData',output.scanNum,iOverlay);
end

% save file
save(fullfile(pathname,filename),output,'-v7.3');

