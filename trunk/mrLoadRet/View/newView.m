function view = newView(viewType)
%
%  view = newView(viewType)
%
% djh, 6/2004

% Define and initialize global variable MLR.
mrGlobals

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

return;
