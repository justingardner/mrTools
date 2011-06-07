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
    end
    % change group
    if isfield(mrLastView,'viewSettings')
      view = viewSet(view,'curGroup',mrLastView.viewSettings.curGroup);
      mlrGuiSet(view.viewNum,'group',mrLastView.viewSettings.curGroup);
    end
    nScans = viewGet(view,'nScans');
    mlrGuiSet(view,'nScans',nScans);
    if baseLoaded && isfield(mrLastView,'viewSettings')
      % slice orientation from last run
      if isfield(mrLastView.viewSettings,'curBase')
	view = viewSet(view,'curBase',mrLastView.viewSettings.curBase);
      end
      view = viewSet(view,'sliceOrientation',mrLastView.viewSettings.sliceOrientation);
      % rotate
      mlrGuiSet(view.viewNum,'rotate',mrLastView.viewSettings.rotate);
      % change scan
      view = viewSet(view,'curScan',mrLastView.viewSettings.curScan);
      % change slice
      view = viewSet(view,'curSlice',mrLastView.viewSettings.curSlice);
    end
    % read analyses
    if isfield(mrLastView.view,'analyses')
      for anum = 1:length(mrLastView.view.analyses)
        view = viewSet(view,'newAnalysis',mrLastView.view.analyses{anum});
%         disppercent(anum /length(mrLastView.view.analyses));
      end
      if anum >= 1
        % overlay settings
        if isfield(mrLastView.viewSettings,'overlayMin')
          mlrGuiSet(view.viewNum,'overlayMin',mrLastView.viewSettings.overlayMin);
          mlrGuiSet(view.viewNum,'overlayMax',mrLastView.viewSettings.overlayMax);
          mlrGuiSet(view.viewNum,'alpha',mrLastView.viewSettings.alpha);
        end
      end
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
    %drawnow        %this doesn't seem useful because it will be done by refreshMLRDisplay
    % read ROIs into current view
    if isfield(mrLastView.view,'ROIs')
      for roinum = 1:length(mrLastView.view.ROIs)
        view = viewSet(view,'newROI',mrLastView.view.ROIs(roinum));
      end
      view = viewSet(view,'currentROI',mrLastView.view.curROI);
      if isfield(mrLastView.viewSettings,'showROIs')
	view = viewSet(view,'showROIs',mrLastView.viewSettings.showROIs);
      end
      if isfield(mrLastView.viewSettings,'labelROIs')
	view = viewSet(view,'labelROIs',mrLastView.viewSettings.labelROIs);
      end
      if isfield(mrLastView.viewSettings,'roiGroup')
	view = viewSet(view,'roiGroup',mrLastView.viewSettings.roiGroup);
      end
    end
    % Add plugins
    if ~isempty(which('mlrPlugin')), view = mlrPlugin(view);end
    
    % add here, to load more info...
    % and refresh
    refreshMLRDisplay(view.viewNum);
  end
%   disppercent(inf);

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
    view = viewSet(view,'sliceOrientation','coronal');
    % rotate 270
    mlrGuiSet(view.viewNum,'rotate',270);
    % change group to last in list
    view = viewSet(view,'curGroup',viewGet(view,'numberOfGroups'));baseLoaded
    % and refresh
  else
    baseLoaded = 0;
  end

