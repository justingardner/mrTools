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

  if viewGet(v,'nROIs') > 0
    % check to see if we have the function getInstances
    if ~isempty(which('getInstances'))
      % put up params for doing instances
      paramsInfo{end+1} = {'getInstances',1,'type=checkbox','Compute and save instances that can be used for doing classification and other analyses '};
      % parameters for average instances computation
      paramsInfo{end+1} = {'computeAverageInstances',1,'type=checkbox','contingent=getInstances','Compute instances using averaging method - fast and dirty method that just takes the average using the below parameter of timepoints after each stimulus presentation. See getInstances for more info'};
      paramsInfo{end+1} = {'startLag',0,'type=numeric','incdec=[-1 1]','minmax=[0 inf]','round=1','contingent=computeAverageInstances','Starglag in volumes to start after timulus presentation when computing averages. 0 is to use default. See getInstances for more info'};
      paramsInfo{end+1} = {'blockLen',0,'type=numeric','incdec=[-1 1]','minmax=[0 inf]','round=1','contingent=computeAverageInstances','Blocklen in volumes to compute average over. 0 is to use default. See getInstances for more info'};
      paramsInfo{end+1} = {'groupTrials',1,'type=numeric','incdec=[-1 1]','minmax=[0 inf]','round=1','contingent=computeAverageInstances','Group trials into k trials and average. Default is 1. See getInstances for more info.'};
      paramsInfo{end+1} = {'minResponseLen',0,'type=numeric','incdec=[-1 1]','minmax=[0 inf]','round=1','contingent=computeAverageInstances','Minimum length a trial must be to be used in averages. Use 0 to default. See getInstances for more info.'};
      % parameters for deconv instnaces
      paramsInfo{end+1} = {'computeDeconvInstances',1,'type=checkbox','contingent=getInstances','Compute instances using deconvolution method. See getInstances for more info.'};
      paramsInfo{end+1} = {'canonicalType',{'allfit2','allfit1','all'},'type=popupmenu','contingent=computeDeconvInstances','How to compute the canonical response for making deconv instances. All is the average deconvolution across all stimulus types and all voxels. Allfit1 fits that with a single gamma and Allfit2 fits that with a difference of gammas. See getCanonical for more info.'};
      paramsInfo{end+1} = {'canonicalR2cutoff',0,'type=numeric','minmax=[0 1]','incdec=[-0.01 0.01]','contingent=computeDeconvInstances','r2 cutoff for voxels to be used in computing canonical'};
    else
      % put up params for doing instances
      paramsInfo{end+1} = {'getInstances',0,'type=checkbox','enable=0','Compute and save instances that can be used for doing classification and other analyses. Not available - download from https://github.com/justingardner/gru.git'};
      disp(sprintf('(mlrExportForAnalysis) Git repo wigh getInstances not available. If you wish to compute instances, download from: https://github.com/justingardner/gru.git'));
    end
  else
    % put up params for doing instances
    paramsInfo{end+1} = {'getInstances',0,'type=checkbox','enable=0','No ROIs are loaded, so it is not possible to compute instances. '};
    disp(sprintf('(mlrExportForAnalysis) No ROIs found. Skipping getInstances'));
  end
  
else
  paramsInfo{end+1} = {'eventRelatedAnalysis',0,'visible=0'};
end

% prf analysis
if ~isempty(analysis) &&  strcmp(analysis.type,'pRFAnal')
  paramsInfo{end+1} = {'pRFAnalysis',1,'type=checkbox','Save pRF analysis'};
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
    output.roi{iROI} = loadROITSeries(v,iROI,output.scanNum,output.groupName);
    % compute mean of time series is we are not saving full
    if ~params.roiFullTimeSeries
      output.roi{iROI}.tSeries = mean(output.roi(iROI).tSeries);
    end
  end
else
  output.nroi = 0;
  output.roi = {};
end

% save event related analysis if it is there
if params.eventRelatedAnalysis
  % get event related
  if length(analysis.d) >= output.scanNum
    output.eventRelated = analysis.d{output.scanNum};
  else
    disp(sprintf('(mlrExportForAnalysis) No event related analysis for %i',output.scanNum));
    output.eventRelated = [];
  end
end

% save pRF Analysis
if params.pRFAnalysis
  % get pRF 
  if length(analysis.d) >= output.scanNum
    output.pRF = analysis.d{output.scanNum};
  else
    disp(sprintf('(mlrExportForAnalysis) No pRF analysis for %i',output.scanNum));
    output.pRF = [];
  end
end

% load overlays
nOverlays = viewGet(v,'nOverlays');
for iOverlay = 1:nOverlays
  % get the overlay
  output.overlay(iOverlay) = viewGet(v,'overlay',iOverlay);
  % extract only the data for this scan
  output.overlay(iOverlay).data = viewGet(v,'overlayData',output.scanNum,iOverlay);
end

% now run instances
if params.eventRelatedAnalysis
  if isequal(params.getInstances,1) && (params.computeAverageInstances || params.computeDeconvInstances)
    % get r2
    r2 = viewGet(v,'overlayData',viewGet(v,'overlayNum','r2'));
    % get sort index, sets r2 for the roi and can be used for sorting voxel order
    output.roi = getSortIndex(v,output.roi,r2);

    % compute average instnaces
    if params.computeAverageInstances
      % set default arguments
      if params.startLag == 0, params.startLag = [];end
      if params.blockLen == 0, params.blockLen = [];end
      if params.minResponseLen == 0, params.minResponseLen = [];end
      % getInstances using average method      
      output.roi = getInstances(v,output.roi,output.eventRelated.stimvol,'startLag',params.startLag,'blockLen',params.blockLen,'groupTrials',params.groupTrials,'minResponseLen',params.minResponseLen,'fieldName=averageInstances','n=inf');
    end

    % compute deconv instnaces
    if params.computeDeconvInstances
      % set default agruments
      if params.canonicalR2cutoff == 0,params.canonicalR2cutoff = [];end
      % getInstances using deconvmethod      
      output.roi = getInstances(v,output.roi,output.eventRelated.stimvol,'type=deconv','canonicalType',params.canonicalType,'r2cutoff',params.canonicalR2cutoff,'fieldName=deconvInstances','n=inf');
    end
  end
end

% save file
save(fullfile(pathname,filename),'output','-v7.3');

