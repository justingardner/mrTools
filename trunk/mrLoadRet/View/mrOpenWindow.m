function view = mrOpenWindow(viewType)
%
%  view = openWindow(viewType)
%
% djh, 6/2004

if ~exist('viewType','var')
    viewType = 'Volume';
end

mrGlobals;

view = newView(viewType);
% view is empty if it failed to initialize
if ~isempty(view)
  fig = mrLoadRetGUI('viewNum',view.viewNum);
  view = viewSet(view,'figure',fig);
else
  return
end
% set the location of the figure
if isfield(MLR.figloc,'mrLoadRetGUI')
  set(fig,'Position',MLR.figloc.mrLoadRetGUI);
end

% Initialize the scan slider
nScans = viewGet(view,'nScans');
mlrGuiSet(view,'nScans',nScans);
mlrGuiSet(view,'scan',min(1,nScans));
% Initialize the slice slider
mlrGuiSet(view,'nSlices',0);
% open an anatomy, if there is one
anatdir = dir('Anatomy/*.img');
if ~isempty(anatdir)
  % load the first anatomy in the list
  view = loadAnat(view,anatdir(1).name);
  view = viewSet(view,'sliceOrientation','coronal');
  % rotate 270
  mlrGuiSet(view.viewNum,'rotate',270);
  % change group to last in list
  view = viewSet(view,'curGroup',viewGet(view,'numberOfGroups'));
  % and refresh
  refreshMLRDisplay(view.viewNum);
end

if isfile('mrLastView.mat')
  mrLastView=load('mrLastView');
  % if the old one exists, then set up fields
  if isfield(mrLastView,'view')
    % read ROIs into current view
    if isfield(mrLastView.view,'ROIs')
      for roinum = 1:length(mrLastView.view.ROIs)
	view = viewSet(view,'newROI',mrLastView.view.ROIs(roinum));
      end
    end
    % read analyses
%    if isfield(mrLastView.view,'analyses')
%      for anum = 1:length(mrLastView.view.analyses)
%	view = viewSet(view,'newAnalysis',mrLastView.view.analyses{anum});
%      end
%    end
    % add here, to load more info...
  end
end

