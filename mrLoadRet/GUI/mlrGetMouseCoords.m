% mlrGetMouseCoords.m
%
%        $Id:$ 
%      usage: coords = mlrGetMouseCoords(v)
%         by: justin gardner
%       date: 12/20/14
%    purpose: Gets where the mouse is on the view. Returns a coords structures
%             with fields for base, scan and tal coordinates and what
%             axis number was clicked. Used by mrInterrogator and mrLoadRetGUI
%
function coords = mlrGetMouseCoords(v)

% init return structure
coords.inAxis = -1;
coords.scan = [];
coords.base = [];
coords.tal = [];
coords.mni = [];

% get the viewNum, globals and test for fig
mrGlobals;
if isnumeric(v)
  % passed in viewNum, get view
  viewNum = v;
  v = viewGet(v,'view');
else
  viewNum = viewGet(v,'viewNum');
end

% check fig
figNum = viewGet(v,'figNum');
if isempty(figNum),return,end

% get figure handles which store information about
% axis and other gui items
handles = guidata(figNum);

% no bases
if viewGet(v,'numBase') == 0
  x = nan;y= nan;s = nan;xBase = nan;yBase = nan;sBase = nan;xTal = nan;yTal = nan;zTal = nan;
  return
end
baseType = viewGet(v,'baseType');

% get location of pointer
pointerLoc = get(handles.axis,'CurrentPoint');

if (baseType == 0) && (viewGet(v,'baseMultiAxis')>=1)
  % for multi axis, we have to check which axis the
  % user clicked on first

  for iAxis = 1:3
    % get current mouse position
    currentPoint = get(handles.sliceAxis(iAxis),'CurrentPoint');
    % get limits on axis
    xLim = get(handles.sliceAxis(iAxis),'xLim');
    yLim = get(handles.sliceAxis(iAxis),'yLim');
    % get the x,y point
    mouseX = round(currentPoint(1,1));
    mouseY = round(currentPoint(1,2));
    % see if we are in bounds
    if (mouseX>0) && (mouseX<=xLim(2)) && (mouseY>0) && (mouseY<=yLim(2))
      % if we are here, then change coords appropriately
      % to redisplay base
      if iAxis==1
	coords.base.x = handles.coords(1);
	coords.base.z = mouseX;
	coords.base.y = mouseY;
	%       handles.coords([3 2]) = [mouseX mouseY];
      elseif iAxis==2
	coords.base.y = handles.coords(2);
	coords.base.z = mouseX;
	coords.base.x = mouseY;
	%       handles.coords([3 1]) = [mouseX mouseY];
      elseif iAxis==3
	coords.base.z = handles.coords(3);
	coords.base.y = mouseX;
	coords.base.x = mouseY;
	%       handles.coords([2 1]) = [mouseX mouseY];
      end
      coords.inAxis = iAxis;
      % found the point, so break.
      break;
    end
  end
elseif baseType <= 1
  mouseY = round(pointerLoc(1,1));
  mouseX = round(pointerLoc(1,2));

  % get base coordinates
  baseCoords = viewGet(v,'cursliceBaseCoords');
  % convert mouse to baseCoords
  if (mouseX>0) && (mouseX<=size(baseCoords,1)) && (mouseY>0) && (mouseY<=size(baseCoords,2))
    coords.base.x = baseCoords(mouseX,mouseY,1);
    coords.base.y = baseCoords(mouseX,mouseY,2);
    coords.base.z = baseCoords(mouseX,mouseY,3);
    % in main axis
    coords.inAxis = 0;
  else
    % out of bounds, so just return with nans set
    return
  end
% get coordinates for surface
else
  % handle getting coordinates for surface
  baseSurface = viewGet(v,'baseSurface');
  baseDims = viewGet(v,'baseSurfaceDims');
  pos = [];
  % check mouse bounding box coords against baseDims
  % for a quick check to see if we are in the volume
  if all(pointerLoc(1,:) >= 0)
    % then use select3d which is slooow, but accurate
    if isfield(MLR.interrogator{viewNum},'axesnum')
      if verLessThan('matlab','8.4')
        hobj = get(MLR.interrogator{viewNum}.axesnum,'Children');
        hobj = hobj(end);
      % make sure we are using the correct object (should be the 3D
      % brain). Use end here because with searchForVoxel we plot a
      % point on the image and the brain object always seems to be
      % last. But if this is not always the case, then we may need
      % to do a little more work here to find the correct object
      else %from  matlab version 8.5a, the handle to the base has been 
        %stored in the view, so use that
        hobj = viewGet(v,'baseHandle');
      end
      [pos, vertex, vertexIndex] = select3d(hobj);
      % convert the index to the coordinates
      if ~isempty(pos)
	baseCoordMap = viewGet(v,'baseCoordMap');
	%we'll take the coordinates of the middle of whatever range of cortical depth is currently selected
	corticalSlice = max(1,ceil(mean(viewGet(v,'corticalDepth'))*size(baseCoordMap.coords,5)));
	pos = round(squeeze(baseCoordMap.coords(1,vertexIndex,1,:,corticalSlice)));
	%pos = round(squeeze(baseCoordMap.coords(1,vertexIndex,1,:)));
	coords.base.x = pos(1);coords.base.y = pos(2);coords.base.z = pos(3);
      end
    end
  end      
end

% if no base coords, then return
if isempty(coords.base),return,end

% transform from base coordinates to talairach coordinates, if the base has a talairach transform defined
base2tal = viewGet(v,'base2tal'); % keyboard
if(~isempty(base2tal))
  talCoords = round(base2tal * [coords.base.x coords.base.y coords.base.z 1]');
  coords.tal.x = talCoords(1); coords.tal.y = talCoords(2); coords.tal.z = talCoords(3);
end

% transform from base coordinates to MNI coordinates. This uses deformation maps computed using SPM
mniInfo = viewGet(v,'mniInfo');
if ~isempty(mniInfo)
  % convert coords from base to T1w volume that was use to compute the deformation coords map
  base2mag = viewGet(v,'base2mag');
  coordsT1w = mniInfo.mag2T1w * base2mag * [coords.base.x coords.base.y coords.base.z 1]';
  % non-linear registration from T1w volume space to mni
  [T1wGridX,T1wGridY,T1wGridZ] = ndgrid(1:size(mniInfo.T1w2mniCoordMap,1),1:size(mniInfo.T1w2mniCoordMap,2),1:size(mniInfo.T1w2mniCoordMap,3));
  coords.mni.x = interpn(T1wGridX,T1wGridY,T1wGridZ, mniInfo.T1w2mniCoordMap(:,:,:,1), coordsT1w(1), coordsT1w(2), coordsT1w(3));
  coords.mni.y = interpn(T1wGridX,T1wGridY,T1wGridZ, mniInfo.T1w2mniCoordMap(:,:,:,2), coordsT1w(1), coordsT1w(2), coordsT1w(3));
  coords.mni.z = interpn(T1wGridX,T1wGridY,T1wGridZ, mniInfo.T1w2mniCoordMap(:,:,:,3), coordsT1w(1), coordsT1w(2), coordsT1w(3));
%   coords.mni.x = mniInfo.T1w2mniCoordMap(round(coordsT1w(1)),round(coordsT1w(2)),round(coordsT1w(3)),1); %
%   coords.mni.y = mniInfo.T1w2mniCoordMap(round(coordsT1w(1)),round(coordsT1w(2)),round(coordsT1w(3)),2); % less accurate but faster?
%   coords.mni.z = mniInfo.T1w2mniCoordMap(round(coordsT1w(1)),round(coordsT1w(2)),round(coordsT1w(3)),3); %
end

% transform from base coordinates into scan coordinates
base2scan = viewGet(v,'base2scan');
if ~isempty(base2scan)
  transformed = round(base2scan*[coords.base.x coords.base.y coords.base.z 1]');
  coords.scan.x = transformed(1);
  coords.scan.y = transformed(2);
  coords.scan.z = transformed(3);

  % get the scan dims to make sure we haven't jumped off end of scan
  scanDims = viewGet(v,'scanDims',viewGet(v,'curScan'));
  if ((coords.scan.x < 1) || (coords.scan.x > scanDims(1)) || ...
      (coords.scan.y < 1) || (coords.scan.y > scanDims(2)) || ...
      (coords.scan.z < 1) || (coords.scan.z > scanDims(3)))
    coords.scan = [];
  end
end

