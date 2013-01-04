function view = mrOpenWindow(viewType,mrLastView)
%
%  view = openWindow(viewType)
%
% djh, 6/2004
%        $Id$

if ieNotDefined('viewType'),viewType = 'Volume';end
% note we don't use ieNotDefined here, because
% if mrLastView is empty then the user doesn't
% want to ignore mrLastView
if ~exist('mrLastView','var')
  mrLastView = 'mrLastView.mat';
end

mrGlobals;

% startup a view
view = newView(viewType);

% view is empty if it failed to initialize
if ~isempty(view)
  fig = mrLoadRetGUI('viewNum',view.viewNum);
  set(fig,'CloseRequestFcn',@mrQuit);
  view = viewSet(view,'figure',fig);
else
  return
end
% set the location of the figure
figloc = mrGetFigLoc('mrLoadRetGUI');
if ~isempty(figloc)
  %deal with multiple monitors
  [whichMonitor,figloc]=getMonitorNumber(figloc,getMonitorPositions);
  set(fig,'Position',figloc);
end

set(fig,'Renderer','painters')
% set the keyoard accelerator
%mrAcceleratorKeys('init',view.viewNum);

% Initialize the scan slider
nScans = viewGet(view,'nScans');
mlrGuiSet(view,'nScans',nScans);
mlrGuiSet(view,'scan',min(1,nScans));
% Initialize the slice slider
mlrGuiSet(view,'nSlices',0);
% init showROIs to all perimeter
view = viewSet(view,'showROIs','all perimeter');
view = viewSet(view,'labelROIs',1);

% Add plugins
if ~isempty(which('mlrPlugin')), view = mlrPlugin(view);end

baseLoaded = 0;
if ~isempty(mrLastView) && isfile(sprintf('%s.mat',stripext(mrLastView)))
  disppercent(-inf,sprintf('(mrOpenWindow) Loading %s',mrLastView));
  mrLastView=load(mrLastView);
  disppercent(inf);
  % if the old one exists, then set up fields
%   disppercent(-inf,'(mrOpenWindow) Restoring last view');
  if isfield(mrLastView,'view')
    % open up base anatomy from last session
    if isfield(mrLastView.view,'baseVolumes')
      disppercent(-inf,sprintf('(mrOpenWindow) installing Base Anatomies'));
      if length(mrLastView.view.baseVolumes) >= 1
        baseLoaded = 1;
        % Add it to the list of base volumes and select it
	for i = 1:length(mrLastView.view.baseVolumes)
	  % make sure sliceOrientation is not 0
	  if ~isfield(mrLastView.view.baseVolumes(i),'sliceOrientation') ...
		|| isequal(mrLastView.view.baseVolumes(i).sliceOrientation,0)
	    mrLastView.view.baseVolumes(i).sliceOrientation = 1;
	  end
	  % install the base
	  view = viewSet(view,'newBase',mrLastView.view.baseVolumes(i));
	end
      else
        %try to load 
  [view,baseLoaded] = loadAnatomy(view);
      end
      disppercent(inf);
    end
    % change group
    view = viewSet(view,'curGroup',mrLastView.view.curGroup);
    nScans = viewGet(view,'nScans');
    mlrGuiSet(view,'nScans',nScans);
    if baseLoaded  && isfield(mrLastView,'viewSettings')
      % slice orientation from last run
      view = viewSet(view,'curBase',mrLastView.view.curBase);
      % change scan
      view = viewSet(view,'curScan',mrLastView.view.curScan);
      % change slice/corticalDepth
      if viewGet(view,'baseType') && isfield(mrLastView.view.curslice,'corticalDepth')
        view = viewSet(view,'corticalDepth',mrLastView.view.curslice.corticalDepth);
      elseif isfield(mrLastView.view.curslice,'sliceNum')
        view = viewSet(view,'curSlice',mrLastView.view.curslice.sliceNum);
      end
      if isfield(mrLastView.viewSettings,'corticalDepth')
	view = viewSet(view,'corticalDepth',mrLastView.viewSettings.corticalDepth);
      end
    end
    % read analyses
    if isfield(mrLastView.view,'analyses')
      for anum = 1:length(mrLastView.view.analyses)
        view = viewSet(view,'newAnalysis',mrLastView.view.analyses{anum});
%         disppercent(anum /length(mrLastView.view.analyses));
      end
      view = viewSet(view,'curAnalysis',mrLastView.view.curAnalysis);
    end
    % read loaded analyses
    if isfield(mrLastView.view,'loadedAnalyses')
      for g = 1:length(mrLastView.view.loadedAnalyses)
	view = viewSet(view,'loadedAnalyses', mrLastView.view.loadedAnalyses{g},g);
      end
    end
    % read which scan we were on in for each group
    if isfield(mrLastView.view,'groupScanNum')
      for g = 1:length(mrLastView.view.groupScanNum)
	view = viewSet(view,'groupScanNum', mrLastView.view.groupScanNum(g),g);
      end
    end
    
    % read ROIs into current view
    if isfield(mrLastView.view,'ROIs')
      disppercent(-inf,sprintf('(mrOpenWindow) installing ROIs'));
      for roinum = 1:length(mrLastView.view.ROIs)
        view = viewSet(view,'newROI',mrLastView.view.ROIs(roinum));
      end
      view = viewSet(view,'currentROI',mrLastView.view.curROI);
      if ~fieldIsNotDefined(mrLastView.view,'showROIs')
	view = viewSet(view,'showROIs',mrLastView.view.showROIs);
      end
      if ~fieldIsNotDefined(mrLastView.view,'labelROIs')
	view = viewSet(view,'labelROIs',mrLastView.view.labelROIs);
      end
      if ~fieldIsNotDefined(mrLastView.view,'roiGroup')
	view = viewSet(view,'roiGroup',mrLastView.view.roiGroup);
      end
      disppercent(inf);
    end
    
    % add here, to load more info...
    % and refresh
    disppercent(-inf,sprintf('(mrOpenWindow) Refreshing MLR display'));
    refreshMLRDisplay(view.viewNum);
    disppercent(inf);
  end

else
  [view,baseLoaded] = loadAnatomy(view);

  if baseLoaded
    refreshMLRDisplay(view.viewNum);
  end
end

% reset some preferences
mrSetPref('importROIPath','');

function [view,baseLoaded] = loadAnatomy(view)
  % for when there is no mrLastView
  % open an anatomy, if there is one
  anatdir = dir('Anatomy/*.img');
  if ~isempty(anatdir)
    baseLoaded = 1;
    % load the first anatomy in the list
    view = loadAnat(view,anatdir(1).name);
    % if it is a regular anatomy
    if viewGet(view,'baseType') == 0
      view = viewSet(view,'sliceOrientation','coronal');
      % set to display a middle slice
      baseDims = viewGet(view,'baseDims');
      baseSliceIndex = viewGet(view,'baseSliceIndex');
      view = viewSet(view,'curSlice',floor(baseDims(baseSliceIndex)/2));
    end
    % change group to last in list
    view = viewSet(view,'curGroup',viewGet(view,'numberOfGroups'));
    % and refresh
  else
    baseLoaded = 0;
  end

