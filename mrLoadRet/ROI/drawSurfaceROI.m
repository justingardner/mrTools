% drawSurfaceROI.m
%
%        $Id$ 
%      usage: drawSurfaceROI()
%         by: Luke Robot Hospadaruk <hospada1@msu.edu>
%       date: 11/28/09
%    purpose: draw an ROI on 3D surface
%
function coords = drawSurfaceROI(view)

%use modified getpts function to grab the correct screen locations
%TODO figure out how to draw stuff on the screen

oneTimeWarning('drawROIOnSurface','(drawROI) Click on each boundary vertex of ROI. To finish making the ROI double click on a point that you want inside the ROI boundary',-1)

hold(get(view.figure, 'CurrentAxes'), 'on');
[xi, yi] = mlr_getpts_3d(view.figure);
hold(get(view.figure, 'CurrentAxes'), 'off');

disp('(drawSurfaceROI) Calculating key vertices');

cmap = view.baseVolumes(view.curBase).coordMap;
roiKeyVertices = [];
%calculate vertex index of each clicked location
for i=1:(numel(xi))
	[x y s xBase yBase sBase xT yT sT vi] = getMouseCoords(view.viewNum, [xi(i) yi(i)]);
	roiKeyVertices = [roiKeyVertices vi];
end
%make sure there are at least 3 perimiter and 1 interior vertices
if(numel(roiKeyVertices) < 4)
	mrWarnDlg('(drawSurfaceROI) Not enough vertices selected, you need at least 4');
  coords = [];
  return
end

disp('(drawSurfaceROI) Setting up boundary search');
%We need to generate a list of connections
nodes = [(1:size(cmap.coords, 2))' squeeze(cmap.coords(1,:,1,1))' squeeze(cmap.coords(1,:,1,2))' squeeze(cmap.coords(1,:,1,2))'];
edges = [(1:size(cmap.tris, 1))' cmap.tris(:,1) cmap.tris(:,2)];
edges = [edges;(1:size(cmap.tris, 1))'+size(cmap.tris, 1) cmap.tris(:,2) cmap.tris(:,3)];
edges = [edges;(1:size(cmap.tris, 1))'+size(cmap.tris, 1)+size(cmap.tris, 1) cmap.tris(:,1) cmap.tris(:,3)];
[tr loc1] = ismember(edges(:,2), nodes(:,1));
[tr loc2] = ismember(edges(:,3), nodes(:,1));
nodes1 = nodes(loc1,2:end);
nodes2 = nodes(loc2,2:end);
distances = sqrt(sum([...
	(nodes1(:,1)-nodes2(:,1)).^2 ...
	(nodes1(:,2)-nodes2(:,2)).^2 ...
	(nodes1(:,3)-nodes2(:,3)).^2], 2));

%make the edges go in both directions, from a-b and b-a to guarantee
%bi-directional edges
edgeDistances = sparse([edges(:,2);edges(:,3)], [edges(:,3);edges(:,2)], [distances;distances], size(cmap.coords, 2), size(cmap.coords, 2));


disp('(drawSurfaceROI) Searching for boundary vertices');

[dist pred] = dijkstrap(edgeDistances, roiKeyVertices);
pred(pred == -1) = 0;

boundaryVertices = [];
for i=1:(numel(roiKeyVertices)-1)
	v1 = roiKeyVertices(i);
	v2 = roiKeyVertices(i+1);
	%if we're at the last boundary point, connect the loop
	if i == (numel(roiKeyVertices)-1)
		v2 = roiKeyVertices(1);
	end
	boundaryVertices = [boundaryVertices pred2path(pred(i,:),v1,v2)];
end
%make sure there are no duplicates
boundaryVertices = unique(boundaryVertices);
seed = roiKeyVertices(end);
%grow the ROI, this will return all vertices, including the boundary
roiVertices = growROI(edgeDistances, boundaryVertices, seed);

baseCoordMap = viewGet(view,'baseCoordMap');
%we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
corticalSlice = ceil(mean(viewGet(view,'corticalDepth'))*size(baseCoordMap.coords,5));
pos = round(squeeze(baseCoordMap.coords(1,roiVertices,1,:,corticalSlice)));
xBase = pos(:,1);yBase = pos(:,2);sBase = pos(:,3);

base2scan = viewGet(view,'base2scan');
if isempty(base2scan), x = nan; y=nan; s=nan; return,  end
transformed = round(base2scan*[xBase yBase sBase ones(length(xBase),1)]');

coords = [xBase yBase sBase ones(length(xBase),1)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%given the distance matrix passed into dijkstra, the boundary vertices
%(which include any previously visited vertices) grow the roi and return
%all vertices.  edge matrix must be bi-directional;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiVertices = growROI(edges, boundaryVertices, seed)
%add the seed and the boundary to the list of visited nodes
roiVertices = [seed boundaryVertices];

%get all vertices connected to seed
seedEdges = find(edges(seed,:));

%select only unvisited, non-boundary edges
seedEdges = seedEdges(~ismember(seedEdges, roiVertices));

disp('(drawSurfaceROI) Growing ROI from interior point');
disppercent(-inf,'(drawSurfaceROI) Growing ROI from interior point. If this is taking a long time, you may have put your interior point on the outside of the boundary. In this case abort with ctrl-c. growROI percentage of total edges: ');
%loop until we can't find any more unvisited nodes
while numel(seedEdges) > 0
  disppercent(length(roiVertices)/length(edges));
	%add the current list of seeds to the return
	roiVertices = [roiVertices seedEdges];
	newSeedEdges = [];
	%for each edge, find all it's connecting edges
	for i=1:numel(seedEdges)
		newSeedEdges = [newSeedEdges find(edges(seedEdges(i),:))];
	end
	%select only those edges which are unvisited
	seedEdges = newSeedEdges(~ismember(newSeedEdges, roiVertices));
	%throw out non-unique nodes
	seedEdges = unique(seedEdges);
	
end
disppercent(inf);
%make sure we arn't returning duplicates
roiVertices = unique(roiVertices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s xBase yBase sBase xTal yTal zTal vi] = getMouseCoords(viewNum, point)

% get the view
mrGlobals;
view = viewGet([],'view',viewNum);

% no bases
if viewGet(view,'numBase') == 0
	x = nan;y= nan;s = nan;xBase = nan;yBase = nan;sBase = nan;xTal = nan;yTal = nan;zTal = nan;
	return
end
baseType = viewGet(view,'baseType');


% get location of pointer
pointerLoc = point;
%set it as the current point so select3d will get the right point
set(view.figure, 'CurrentPoint', point);
%pointerLoc = get(MLR.interrogator{viewNum}.axesnum,'CurrentPoint');

if baseType <= 1
	mouseY = round(pointerLoc(1,1));
	mouseX = round(pointerLoc(1,2));
	
	% get base coordinates
	baseCoords = viewGet(view,'cursliceBaseCoords');
	% convert mouse to baseCoords
	if (mouseX>0) && (mouseX<=size(baseCoords,1)) && (mouseY>0) && (mouseY<=size(baseCoords,2))
		xBase = baseCoords(mouseX,mouseY,1);
		yBase = baseCoords(mouseX,mouseY,2);
		sBase = baseCoords(mouseX,mouseY,3);
	else
		x = nan;y = nan; s = nan;
		xBase = nan;yBase = nan; sBase = nan;
		xTal = nan; yTal = nan; zTal = nan;
		return
	end
else
	% handle getting coordinates for surface
	baseSurface = viewGet(view,'baseSurface');
	baseDims = viewGet(view,'baseSurfaceDims');
	pos = [];xBase = nan; yBase = nan; sBase = nan;
	% check mouse bounding box coords against baseDims
	% for a quick check to see if we are in the volume
	if all(pointerLoc(1,:) >= 0)
		% then use select3d which is slooow, but accurate
		hobj = get(get(view.figure, 'CurrentAxes'),'Children');
		% make sure we are using the correct object (should be the 3D
		% brain). Use end here because with searchForVoxel we plot a
		% point on the image and the brain object always seems to be
		% last. But if this is not always the case, then we may need
		% to do a little more work here to find the correct object
		[pos v vi] = select3d(hobj(end));
		% convert the index to the coordinates
		if ~isempty(pos)
			baseCoordMap = viewGet(view,'baseCoordMap');
      %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
      corticalSlice = ceil(mean(viewGet(view,'corticalDepth'))*size(baseCoordMap.coords,5));
			pos = round(squeeze(baseCoordMap.coords(1,vi,1,:,corticalSlice)));
			xBase = pos(1);yBase = pos(2);sBase = pos(3);
		end
	end
end

% transform from base coordinates to talairach coordinates, if the base has a talairach transform defined
base2tal = viewGet(view,'base2tal'); % keyboard
if(~isempty(base2tal))
	talCoords = round(base2tal * [xBase yBase sBase 1]');
	xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
else
	xTal = nan; yTal = nan; zTal = nan;
end


% transform from base coordinates into scan coordinates
base2scan = viewGet(view,'base2scan');
if isempty(base2scan), x = nan; y=nan; s=nan; return,  end
transformed = round(base2scan*[xBase yBase sBase 1]');

x = transformed(1);
y = transformed(2);
s = transformed(3);

% get the scan dims to make sure we haven't jumped off end
% of scan
scanDims = viewGet(view,'scanDims',viewGet(view,'curScan'));
if ((x < 1) || (x > scanDims(1)) || ...
		(y < 1) || (y > scanDims(2)) || ...
		(s < 1) || (s > scanDims(3)))
	x = nan;y = nan;s = nan;
end


function [x,y] = mlr_getpts_3d(varargin)
%Modified version to get points for getMouseCoords from 3d surface
%for ROI creation on surfaces
%GETPTS Select points with mouse.
%   [X,Y] = GETPTS(FIG) lets you choose a set of points in the
%   current axes of figure FIG using the mouse. Coordinates of
%   the selected points are returned in the vectors X and Y. Use
%   normal button clicks to add points.  A shift-, right-, or
%   double-click adds a final point and ends the selection.
%   Pressing RETURN or ENTER ends the selection without adding
%   a final point.  Pressing BACKSPACE or DELETE removes the
%   previously selected point.
%
%   [X,Y] = GETPTS(AX) lets you choose points in the axes
%   specified by the handle AX.
%
%   [X,Y] = GETPTS is the same as [X,Y] = GETPTS(GCF).
%
%   Example
%   --------
%       imshow('moon.tif')
%       [x,y] = getpts_luke
%
%   See also GETRECT, GETLINE.

%   Callback syntaxes:
%       getpts_luke('KeyPress')
%       getpts_luke('FirstButtonDown')
%       getpts_luke('NextButtonDown')

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1533 $  $Date: 2009-10-07 02:26:45 +0900 (Wed, 07 Oct 2009) $
%	Modified by: Luke Hospadaruk

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
global GETPTS_PT1 line_points
if ((nargin >= 1) && (ischar(varargin{1})))
	% Callback invocation: 'KeyPress', 'FirstButtonDown', or
	% 'NextButtonDown'.
	feval(varargin{:});
	return;
end

if (nargin < 1)
	GETPTS_AX = gca;
	GETPTS_FIG = ancestor(GETPTS_AX, 'figure');
else
	if (~ishandle(varargin{1}))
		eid = 'Images:getpts:expectedHandle';
		error(eid, '%s', 'First argument is not a valid handle');
	end
	
	switch get(varargin{1}, 'Type')
		case 'figure'
			GETPTS_FIG = varargin{1};
			GETPTS_AX = get(GETPTS_FIG, 'CurrentAxes');
			if (isempty(GETPTS_AX))
				GETPTS_AX = axes('Parent', GETPTS_FIG);
			end
			
		case 'axes'
			GETPTS_AX = varargin{1};
			GETPTS_FIG = ancestor(GETPTS_AX, 'figure');
			
		otherwise
			eid = 'Images:getpts:expectedFigureOrAxesHandle';
			error(eid, '%s', 'First argument should be a figure or axes handle');
			
	end
end

% Bring target figure forward
figure(GETPTS_FIG);

% Remember initial figure state
state = uisuspend(GETPTS_FIG);

% Set up initial callbacks for initial stage
[pointerShape, pointerHotSpot] = CreatePointer;
set(GETPTS_FIG, 'WindowButtonDownFcn', 'mlr_getpts_3d(''FirstButtonDown'');', ...
	'KeyPressFcn', 'mlr_getpts_3d(''KeyPress'');', ...
	'Pointer', 'custom', ...
	'PointerShapeCData', pointerShape, ...
	'PointerShapeHotSpot', pointerHotSpot);

% Initialize the lines to be used for the drag
markerSize = 9;
GETPTS_H1 = line('Parent', GETPTS_AX, ...
	'XData', [], ...
	'YData', [], ...
	'Visible', 'off', ...
	'Clipping', 'off', ...
	'Color', 'c', ...
	'LineStyle', 'none', ...
	'Marker', '+', ...
	'MarkerSize', markerSize);

GETPTS_H2 = line('Parent', GETPTS_AX, ...
	'XData', [], ...
	'YData', [], ...
	'Visible', 'off', ...
	'Clipping', 'off', ...
	'Color', 'm', ...
	'LineStyle', 'none', ...
	'Marker', 'x', ...
	'MarkerSize', markerSize);

% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in try-catch so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
try
	waitfor(GETPTS_H1, 'UserData', 'Completed');
catch
	errCatch=1;
end

% After the waitfor, if GETPTS_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
	errStatus = 'trap';
	
elseif (~ishandle(GETPTS_H1) || ...
		~strcmp(get(GETPTS_H1, 'UserData'), 'Completed'))
	errStatus = 'unknown';
	
else
	errStatus = 'ok';
	x = get(GETPTS_H1, 'XData');
	y = get(GETPTS_H1, 'YData');
	x = x(:);
	y = y(:);
	% If no points were selected, return rectangular empties.
	% This makes it easier to handle degenerate cases in
	% functions that call getpts_luke.
	if (isempty(x))
		x = zeros(0,1);
	end
	if (isempty(y))
		y = zeros(0,1);
	end
end

% Delete the animation objects
if (ishandle(GETPTS_H1))
	delete(GETPTS_H1);
end
if (ishandle(GETPTS_H2))
	delete(GETPTS_H2);
end

% Restore the figure state
if (ishandle(GETPTS_FIG))
	uirestore(state);
end

% Clean up the global workspace
clear global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
clear global GETPTS_PT1 line_points

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
	case 'ok'
		% No action needed.
		
	case 'trap'
		% An error was trapped during the waitfor
		eid = 'Images:getpts:interruptedMouseSelection';
		error(eid, '%s', 'Interruption during mouse point selection.');
		
	case 'unknown'
		% User did something to cause the point selection to
		% terminate abnormally.  For example, we would get here
		% if the user closed the figure in the middle of the selection.
		eid = 'Images:getpts:interruptedMouseSelection';
		error(eid, '%s', 'Interruption during mouse point selection.');
end


%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
global GETPTS_PT1 line_points

key = get(GETPTS_FIG, 'CurrentCharacter');
switch key
	case {char(8), char(127)}  % delete and backspace keys
		x = get(GETPTS_H1, 'XData');
		y = get(GETPTS_H1, 'YData');
		switch length(x)
			case 0
				% nothing to do
			case 1
				% remove point and start over
				set([GETPTS_H1 GETPTS_H2], ...
					'XData', [], ...
					'YData', []);
				set(GETPTS_FIG, 'WindowButtonDownFcn', ...
					'getpts_luke(''FirstButtonDown'');');
			otherwise
				% remove last point
				set([GETPTS_H1 GETPTS_H2], ...
					'XData', x(1:end-1), ...
					'YData', y(1:end-1));
		end
		
	case {char(13), char(3)}   % enter and return keys
		% return control to line after waitfor
		set(GETPTS_H1, 'UserData', 'Completed');
		
end

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 line_points

%[x,y] = getcurpt(GETPTS_AX);
pos = get(GETPTS_FIG, 'CurrentPoint');
x = pos(1);
y = pos(2);
view = getMLRView;
% then use select3d which is slooow, but accurate
hobj = get(get(view.figure, 'CurrentAxes'),'Children');
% make sure we are using the correct object (should be the 3D
% brain). Use end here because with searchForVoxel we plot a
% point on the image and the brain object always seems to be
% last. But if this is not always the case, then we may need
% to do a little more work here to find the correct object
[pos v vi] = select3d(hobj(end));
% convert the index to the coordinates

if ~isempty(pos)
	line_points =[line_points pos];
	if(numel(line_points) > 3)
		plot3(line_points(1,:), line_points(2,:), line_points(3,:),'Parent', get(view.figure, 'CurrentAxes'));
	end
end
set([GETPTS_H1 GETPTS_H2], ...
	'XData', x, ...
	'YData', y, ...
	'Visible', 'on');

if (~strcmp(get(GETPTS_FIG, 'SelectionType'), 'normal'))
	% We're done!
	set(GETPTS_H1, 'UserData', 'Completed');
else
	set(GETPTS_FIG, 'WindowButtonDownFcn', 'mlr_getpts_3d(''NextButtonDown'');');
end

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 line_points

selectionType = get(GETPTS_FIG, 'SelectionType');
if (~strcmp(selectionType, 'open'))
	% We don't want to add a point on the second click
	% of a double-click
	
	%[newx, newy] = getcurpt(GETPTS_AX);
	pos = get(GETPTS_FIG, 'CurrentPoint');
	newx = pos(1);
	newy = pos(2);
	view = getMLRView;
	% then use select3d which is slooow, but accurate
	hobj = get(get(view.figure, 'CurrentAxes'),'Children');
	% make sure we are using the correct object (should be the 3D
	% brain). Use end here because with searchForVoxel we plot a
	% point on the image and the brain object always seems to be
	% last. But if this is not always the case, then we may need
	% to do a little more work here to find the correct object
	[pos v vi] = select3d(hobj(end));
	% convert the index to the coordinates
	
	if ~isempty(pos)
		line_points =[line_points pos];
		if(numel(line_points) > 3)
			plot3(line_points(1,:), line_points(2,:), line_points(3,:), 'Parent', get(view.figure, 'CurrentAxes'));
		end
	end
	
	
	
	
	x = get(GETPTS_H1, 'XData');
	y = get(GETPTS_H2, 'YData');
	
	set([GETPTS_H1 GETPTS_H2], 'XData', [x newx], ...
		'YData', [y newy]);
	
end

if (~strcmp(get(GETPTS_FIG, 'SelectionType'), 'normal'))
	% We're done!
	set(GETPTS_H1, 'UserData', 'Completed');
end



%----------------------------------------------------
% Subfunction CreatePointer
%----------------------------------------------------
function [pointerShape, pointerHotSpot] = CreatePointer

pointerHotSpot = [8 8];
pointerShape = [ ...
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
	2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
	2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
	1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
