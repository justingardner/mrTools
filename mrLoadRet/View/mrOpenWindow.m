function view = mrOpenWindow(viewType,mrLastView,noGUI)
%
%  view = openWindow(viewType)
%
% djh, 6/2004
%        $Id: mrOpenWindow.m 2838 2013-08-12 12:52:20Z julien $

if ieNotDefined('viewType'),viewType = 'Volume';end
if ieNotDefined('noGUI'), noGUI = false;end

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
if isempty(view)
  return
end

if noGUI
  % skip all the GUI-related stuff
else
  
  fig = mrLoadRetGUI('viewNum',view.viewNum);
  set(fig,'CloseRequestFcn',@mrQuit);
  view = viewSet(view,'figure',fig);
  
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

  % set the position of the main base viewer
  gui = guidata(fig);
  gui.marginSize = 0.01;
  gui.anatPosition = [0.3 0.2+gui.marginSize 1-0.3-gui.marginSize 1-0.2-2*gui.marginSize];

  % create 3 axis for display all three orientations at once
  % set up the position of each of the 3 axis. Start them
  % in an arbitrary position, being careful not to overlap
  gui.sliceAxis(1) = subplot('Position',[0 0 0.01 0.01],'Parent',fig);
  axis(gui.sliceAxis(1),'off');
  set(gui.sliceAxis(1),'HandleVisibility','off');
  gui.sliceAxis(2) = subplot('Position',[0.02 0 0.01 0.01],'Parent',fig);
  axis(gui.sliceAxis(2),'off');
  set(gui.sliceAxis(2),'HandleVisibility','off');
  gui.sliceAxis(3) = subplot('Position',[0.02 0.02 0.01 0.01],'Parent',fig);
  axis(gui.sliceAxis(3),'off');
  set(gui.sliceAxis(3),'HandleVisibility','off');

  % save the axis handles
  guidata(fig,gui);

  % add controls for multiAxis
  mlrAdjustGUI(view,'add','control','axisSingle','style','radio','value',0,'position',  [0.152    0.725   0.1    0.025],'String','Single','Callback',@multiAxisCallback);
  mlrAdjustGUI(view,'add','control','axisMulti','style','radio','value',0,'position',  [0.152    0.695   0.1    0.025],'String','Multi','Callback',@multiAxisCallback);
  mlrAdjustGUI(view,'add','control','axis3D','style','radio','value',0,'position',  [0.152    0.665   0.1    0.025],'String','3D','Callback',@multiAxisCallback);

  % Initialize the scan slider
  nScans = viewGet(view,'nScans');
  mlrGuiSet(view,'nScans',nScans);
  mlrGuiSet(view,'scan',min(1,nScans));
  % Initialize the slice slider
  mlrGuiSet(view,'nSlices',0);
  % init showROIs to all perimeter
  view = viewSet(view,'showROIs','all perimeter');
  view = viewSet(view,'labelROIs',1);

  %get colormaps in the colormapFunctions directory
  colorMapsFolder = [fileparts(which('mrLoadRet')) '/colormapFunctions/'];
  colorMapFiles =  dir([colorMapsFolder '*.m']);
  if ~isempty(colorMapFiles)
    colorMapList = cell(1,length(colorMapFiles));
    for iFile=1:length(colorMapFiles)
       colorMapList{iFile} = stripext(colorMapFiles(iFile).name);
    end
    % install default colormaps
    % that will show up when you do /Edit/Overlay
    mlrAdjustGUI(view,'add','colormap',colorMapList);
  else
    disp(['(mrOpenWindow) No colormap found in folder ' colorMapsFolder]);
  end

  % add a menu item to export analysis struct
  mlrAdjustGUI(view,'add','menu','Export for Analysis','/File/Export/','Callback',@mlrExportForAnalysis);

  % add a menu item to export surfaces to wavefront off
  mlrAdjustGUI(view,'add','menu','Export surface','/File/Base anatomy/Use current scan','Callback',@mlrExportSurface,'Separator','on');
  mlrAdjustGUI(view,'set','Export surface','Enable','off');


  % Add plugins
  if ~isempty(which('mlrPlugin')), view = mlrPlugin(view);end
end

baseLoaded = 0;
if ~isempty(mrLastView) && mlrIsFile(sprintf('%s.mat',stripext(mrLastView)))
  mlrDispPercent(-inf,sprintf('(mrOpenWindow) Loading %s',mrLastView));
  [mrLastView, lastViewSettings]=mlrLoadLastView(mrLastView);
  mlrDispPercent(inf);
  % if the old one exists, then set up fields
%   mlrDispPercent(-inf,'(mrOpenWindow) Restoring last view');
  if ~isempty(mrLastView)
    %Add any missing field to make sure things don't crash
    [~,mrLastView,unknownFields]=isview(mrLastView);
    if ~isempty(unknownFields) %warn if there are any unknown fields
      fprintf('(mrOpenWindow) Unknown view field(s):')
      for i = unknownFields'
        fprintf(' %s',i{1});
      end
      fprintf('\n');
    end
    % open up base anatomy from last session
    if isfield(mrLastView,'baseVolumes')
      mlrDispPercent(-inf,sprintf('(mrOpenWindow) installing Base Anatomies'));
      if length(mrLastView.baseVolumes) >= 1
        baseLoaded = 1;
        % Add it to the list of base volumes and select it
	for i = 1:length(mrLastView.baseVolumes)
	  % make sure sliceOrientation is not 0
	  if ~isfield(mrLastView.baseVolumes(i),'sliceOrientation') ...
		|| isequal(mrLastView.baseVolumes(i).sliceOrientation,0)
	    mrLastView.baseVolumes(i).sliceOrientation = 1;
	  end
	  % install the base
	  view = viewSet(view,'newBase',mrLastView.baseVolumes(i));
	end
      else
        %try to load 
	[view,baseLoaded] = loadAnatomy(view);
      end
      mlrDispPercent(inf);
    end
    % change group
    view = viewSet(view,'curGroup',mrLastView.curGroup);
    if ~noGUI
      nScans = viewGet(view,'nScans');
      mlrGuiSet(view,'nScans',nScans); % JB: I dont' think this call to mlrGuiSet belongs in this function.
      % It's probably not even necessary since the number of scans should have been set by viewSet(v,'curGroup')
    end
    if baseLoaded
      % slice orientation from last run
      view = viewSet(view,'curBase',mrLastView.curBase);
      % change scan
      view = viewSet(view,'curScan',mrLastView.curScan);
      % set displayGyrusSulcusBoundary on/off
      view = viewSet(view,'displayGyrusSulcusBoundary',viewGet(mrLastView,'displayGyrusSulcusBoundary'));
    end
    % read analyses
    if isfield(mrLastView,'analyses')
      for anum = 1:length(mrLastView.analyses)
        view = viewSet(view,'newAnalysis',mrLastView.analyses{anum});
%         mlrDispPercent(anum /length(mrLastView.analyses));
      end
      view = viewSet(view,'curAnalysis',mrLastView.curAnalysis);
    end
    % read loaded analyses
    if isfield(mrLastView,'loadedAnalyses')
      for g = 1:length(mrLastView.loadedAnalyses)
	view = viewSet(view,'loadedAnalyses', mrLastView.loadedAnalyses{g},g);
      end
    end
    % read which scan we were on in for each group
    if isfield(mrLastView,'groupScanNum')
      for g = 1:length(mrLastView.groupScanNum)
	view = viewSet(view,'groupScanNum', mrLastView.groupScanNum(g),g);
      end
    end
    
    % read ROIs into current view
    if isfield(mrLastView,'ROIs')
      mlrDispPercent(-inf,sprintf('(mrOpenWindow) installing ROIs'));
      for roinum = 1:length(mrLastView.ROIs)
        view = viewSet(view,'newROI',mrLastView.ROIs(roinum));
      end
      view = viewSet(view,'currentROI',mrLastView.curROI);
      if ~fieldIsNotDefined(mrLastView,'showROIs')
	view = viewSet(view,'showROIs',mrLastView.showROIs);
      end
      if ~fieldIsNotDefined(mrLastView,'labelROIs')
	view = viewSet(view,'labelROIs',mrLastView.labelROIs);
      end
      if ~fieldIsNotDefined(mrLastView,'roiGroup')
	view = viewSet(view,'roiGroup',mrLastView.roiGroup);
      end
      mlrDispPercent(inf);
    end
    
    if noGUI
      % skip GUI stuff
    else
      % check panels that need to be hidden
      if ~isempty(lastViewSettings) && isfield(lastViewSettings,'panels')
        for iPanel = 1:length(lastViewSettings.panels)
          % if it is not displaying then turn it off
          if ~lastViewSettings.panels{iPanel}{4}
            panelName = lastViewSettings.panels{iPanel}{1};
            mlrGuiSet(view,'hidePanel',panelName);
            % turn off check
            mlrAdjustGUI(view,'set',panelName,'Checked','off');
          end
        end
      end
      % add here, to load more info...
      % and refresh
      mlrDispPercent(-inf,sprintf('(mrOpenWindow) Refreshing MLR display'));
      refreshMLRDisplay(view.viewNum);
      mlrDispPercent(inf);
    end
  end

else
  [view,baseLoaded] = loadAnatomy(view);

  if ~noGUI && baseLoaded
    refreshMLRDisplay(view.viewNum);
  end
end

% reset some preferences
mrSetPref('importROIPath','');

%%%%%%%%%%%%%%%%%%%%%
%    loadAnatomy    %
%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiAxisCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function multiAxisCallback(hObject,eventdata)

% get gui data and v
gui = guidata(get(hObject,'Parent'));
v = viewGet(gui.viewNum,'view');

% this function gets called when single/Multi/3D radio button
% get called. 

% First, determine who was hit
controlNames = {'Single','Multi','3D'};
multiAxisVal = find(strcmp(controlNames,get(hObject,'String')));

% now set the curent base val
viewSet(v,'baseMultiAxis',multiAxisVal-1);

% redraw
refreshMLRDisplay(v.viewNum);



  
