function view = drawROI(view,descriptor,sgn)

% view = drawROI(view,[descriptor],[sgn])
%
% Adds/deletes to/from the current ROI based on user interaction.
%
% descriptor: option for how the new coordinates are to be specified.
% Current options are:
%    'rectangle'[default]
%
% sgn: If sgn~=0 [default, adds user-specified coordinates to selected ROI
% in current slice. If sgn==0, removes those coordinates from the ROI.
%
% $Id$	
% djh, 8/2005 (modified from mrLoadRet-3.1)

% Error if no current ROI
curROInum = viewGet(view,'currentROI');
if isempty(curROInum)
  mrWarnDlg('No ROI currently selected.');
  return
end

if ieNotDefined('descriptor')
  descriptor = 'rectangle'
end
if ieNotDefined('sgn')
  sgn = 1;
end

% Select main axes of view figure for user input
fig = viewGet(view,'figNum');
gui = guidata(fig);
set(fig,'CurrentAxes',gui.axis);

% baseCoords contains the mapping from pixels in the displayed slice to
% voxels in the current base volume.
baseCoords = viewGet(view,'cursliceBaseCoords');
baseSliceDims = [size(baseCoords,1),size(baseCoords,2)];
if isempty(baseCoords)
  mrErrorDlg('Load base anatomy before drawing an ROI');
end

% unhook existing button down/up/move functions -- to keep mrInterrogator from
% running when we define ROIs, but remember the callbacks so we can put them back on
windowButtonDownFcn = get(fig,'WindowButtonDownFcn');
windowButtonMotionFcn = get(fig,'WindowButtonMotionFcn');
windowButtonUpFcn = get(fig,'WindowButtonUpFcn');
set(fig,'WindowButtonDownFcn','');
set(fig,'WindowButtonMotionFcn','');
set(fig,'WindowButtonUpFcn','');

switch descriptor
  case 'rectangle'
    % Get region from user.
    region = round(ginput(2));
    
    % Note: ginput hands them back in x, y order (1st col is x and 2nd col is
    % y). But we use them in the opposite order (y then x), so flip 'em.
    region = fliplr(region);
    % Check if outside image
    if (min(region(:,1))< 1 | max(region(:,1))>baseSliceDims(1) | ...
        min(region(:,2))< 1 | max(region(:,2))>baseSliceDims(2))
      mrWarnDlg('Must choose rect entirely within image boundaries');
      % reset button down/up/move functions
      set(fig,'WindowButtonDownFcn',windowButtonDownFcn);
      set(fig,'WindowButtonMotionFcn',windowButtonMotionFcn);
      set(fig,'WindowButtonUpFcn',windowButtonUpFcn);
      return;
    end
    % Make sure 2nd value is larger than the 1st.
    for i=1:2
      if region(2,i) < region(1,i)
        region(:,i)=flipud(region(:,i));
      end
    end
    
    % Extract coordinates in base reference frame
    baseX = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],1);
    baseY = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],2);
    baseZ = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],3);
    coords = [baseX(:)'; baseY(:)'; baseZ(:)'; ones(1,prod(size(baseX)))];
    
 case 'polygon'
    roiPolygonMethod = mrGetPref('roiPolygonMethod');
    % Get polygon region using matlab's roipoly function
    if strcmp(roiPolygonMethod,'roipoly')
      % this is sometimes very slow if you have a lot
      % of lines already drawn on the figure.
      % i.e. if you have rois already being displayed
      polyIm = roipoly;
    else
      % this doesn't have to redraw lines all the time
      % so it is faster
      % but has the disadvantage that you don't get
      % to see the lines connecting the points.
      [x y a] = getimage;
      if strcmp(roiPolygonMethod,'getptsNoDoubleClick')
	[xi yi] = getptsNoDoubleClick;
      else
	[xi yi] = getpts;
      end
      % draw the lines temporarily
      if ~isempty(xi)
        %	    line([xi;xi(1)],[yi;yi(1)]);
        %	    drawnow;
      end
      polyIm = roipoly(a,xi,yi);
    end
    % Extract coordinates in base reference frame
    baseX = baseCoords(:,:,1);
    baseY = baseCoords(:,:,2);
    baseZ = baseCoords(:,:,3);
    polyImIndices = find(polyIm);
    x = baseX(polyImIndices);
    y = baseY(polyImIndices);
    z = baseZ(polyImIndices);
    coords = [x'; y'; z'; ones(1,length(x))];

  case 'line'
    % grab two points from the image;
    [xi yi] = getpts;

    xii=[]; yii=[];
    for p=1:length(xi)-1
      [tmpX, tmpY] = findLinePoints([xi(p) yi(p)], [xi(p+1) yi(p+1)]); 
      xii = [xii tmpX];
      yii = [yii tmpY];
    end

    if ~isempty(xii)
      line(xii, yii);
      drawnow;
    end

    baseX = baseCoords(:,:,1);
    baseY = baseCoords(:,:,2);
    baseZ = baseCoords(:,:,3);
    lineInd = sub2ind(size(baseX), round(yii), round(xii));
    x = baseX(lineInd);
    y = baseY(lineInd);
    z = baseZ(lineInd);
    
    coords = [x; y; z; ones(1,length(x))];

  otherwise
    mrErrorDlg(['Invalid descriptor: ',descriptor]);
end

% reset button down/up/move functions
set(fig,'WindowButtonDownFcn',windowButtonDownFcn);
set(fig,'WindowButtonMotionFcn',windowButtonMotionFcn);
set(fig,'WindowButtonUpFcn',windowButtonUpFcn);

% Modify ROI
baseNum = viewGet(view,'currentBase');
base2roi = viewGet(view,'base2roi',[],baseNum);
voxelSize = viewGet(view,'baseVoxelSize',baseNum);
view = modifyROI(view,coords,base2roi,voxelSize,sgn);


