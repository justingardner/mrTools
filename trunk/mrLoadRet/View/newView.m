function view = newView(viewType)
%
%  v = newView
%
% djh, 6/2004

if nargin == 0
  viewType = 'Volume';
end

% Define and initialize global variable MLR.
mrGlobals
if isempty(MLR.session),return,end

viewNum = length(MLR.views) + 1;
view.viewNum = viewNum;
view.viewType = viewType;

% Initialize anat
view.baseVolumes = struct([]);
view.curBase = [];

% Initialize analysis list
view.analyses = {};
view.curAnalysis = [];

% Initialize ROIs
view.ROIs = struct([]);
view.curROI = [];
view.prevROIcoords = [];
view.showROIs = 'all';

% Initialize curGroup
view.curGroup = 1;

% Figure handle
view.figure = [];

% Coordinates corresponding to the slice currently displayed
view.curslice.baseCoords = [];
view.curslice.overlayCoords = [];

% Add the new view to the list of views
MLR.views{viewNum} = view;

% validate view (add any optional fields);
[tf view] = isview(view);

% add it back to the globals, in case it has changed
MLR.views{viewNum} = view;

% add the caches
MLR.caches{viewNum} = [];

% init the caches
view = viewSet(view,'roiCache','init');
view = viewSet(view,'overlayCache','init');
view = viewSet(view,'baseCache','init');


return;
