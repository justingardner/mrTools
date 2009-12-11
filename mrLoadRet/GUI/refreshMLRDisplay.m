function [img base roi overlay] = refreshMLRDisplay(viewNum)
%	$Id$

mrGlobals
img = [];base = [];roi = [];overlay=[];
% for debugging/performance tests
% set to 0 for no info
% set to 1 for info on caching
% set to 2 for all info
% if you want to comment out all of this, replace the string
% "if verbose" to "%if verbose"
verbose = 0;
if verbose,tic,end

% Get current view and baseNum.
% Get interp preferences.
% Get slice, scan, alpha, rotate, and sliceIndex from the gui.
if verbose>1,disppercent(-inf,'viewGet');,end
view = viewGet(viewNum,'view');
fig = viewGet(view,'figNum');
slice = viewGet(view,'curslice');
rotate = viewGet(view,'rotate');
baseNum = viewGet(view,'currentBase');
sliceIndex = viewGet(view,'baseSliceIndex',baseNum);
baseType = viewGet(view,'baseType',baseNum);
baseGamma = viewGet(view,'baseGamma',baseNum);
% if no base then clear axis and return
if isempty(baseNum)
  fig = viewGet(view,'figNum');
  if ~isempty(fig)
    gui = guidata(fig);
    cla(gui.axis,'reset');
    set(fig,'CurrentAxes',gui.colorbar);
    cla(gui.colorbar,'reset');
    axis(gui.axis,'off');
    axis(gui.colorbar,'off');
  end
  return
end
if verbose>1,disppercent(inf);,end

% for debugging, clears caches when user holds down shift key
if ~isempty(fig) && any(strcmp(get(fig,'CurrentModifier'),'shift'))
  disp(sprintf('(refreshMLRDisplay) Dumping caches'));
  view = viewSet(view,'roiCache','init');
  view = viewSet(view,'overlayCache','init');
  view = viewSet(view,'baseCache','init');
end

% Compute base coordinates and extract baseIm for the current slice
if verbose,disppercent(-inf,'extract base image');,end
base = viewGet(view,'baseCache');
if isempty(base)
  [base.im,base.coords,base.coordsHomogeneous] = ...
    getBaseSlice(view,slice,sliceIndex,rotate,baseNum,baseType);
  base.dims = size(base.im);

  % Rescale base volume
  if ~isempty(base.im)
    base.cmap = gray(256);
    base.clip = viewGet(view,'baseClip',baseNum);
    base.RGB = rescale2rgb(base.im,base.cmap,base.clip,baseGamma);
  else
    base.RGB = [];
  end
  % save extracted image
  view = viewSet(view,'baseCache',base);
  if verbose,disppercent(inf);disp('Recomputed base');end
else
  if verbose,disppercent(inf);end
end
view = viewSet(view,'cursliceBaseCoords',base.coords);

% Extract overlay images and overlay coords, and alphaMap
if verbose,disppercent(-inf,'extract overlay images');end
overlay = viewGet(view,'overlayCache');
if isempty(overlay)
  % get the transform from the base to the scan
  base2scan = viewGet(view,'base2scan',[],[],baseNum);
  % compute the overlay
  overlay = computeOverlay(view,base2scan,base.coordsHomogeneous,base.dims);
  % save in cache
  view = viewSet(view,'overlayCache',overlay);
  if verbose,disppercent(inf);disp('Recomputed overlay');end
else
  if verbose,disppercent(inf);end
end
view = viewSet(view,'cursliceOverlayCoords',overlay.coords);

% figure
% image(overlayRGB)
% colormap(overlayCmap)
% return

% Combine base and overlay
if verbose>1,disppercent(-inf,'combine base and overlay');,end
if ~isempty(base.RGB) & ~isempty(overlay.RGB)
  img = (1-overlay.alphaMap).*base.RGB + overlay.alphaMap.*overlay.RGB;
  cmap = overlay.cmap;
  cbarRange = overlay.range;
elseif ~isempty(base.RGB)
  img = base.RGB;
  cmap = base.cmap;
  cbarRange = base.clip;
else
  % If no image at this point then display blank
  img = 0;
  cmap = gray(1);
  cbarRange = [0 1];
end
if verbose>1,disppercent(inf);,end

% If no image at this point then return
if ieNotDefined('img')
  return
end

% if no figure, then just return, this is being called just to create the overlay
fig = viewGet(view,'figNum');
if isempty(fig),return,end

% Display the image
if verbose>1,disppercent(-inf,'Clearing figure');,end
gui = guidata(fig);
% note: This cla here is VERY important. Otherwise
% we keep drawing over old things on the axis and
% the rendering gets impossibly slow... -j.
cla(gui.axis);
if verbose>1,disppercent(inf);,end
if verbose>1,disppercent(-inf,'Displaying image');,end
if baseType <= 1
  % set the renderer to painters (this seems
  % to avoid some weird gliches in the OpenGL
  % renderer. It also appears about 20ms or so
  % faster for displaying images as opposed
  % to the 3D surfaces.
  set(fig,'Renderer','painters')
  image(img,'Parent',gui.axis);
else
  % set the renderer to OpenGL, this makes rendering
  % *much* faster -- from about 30 seconds to 30ms
  % it is also marginally faster than the zbuffer
  % renderer (order of 5-10 ms)
  set(fig,'Renderer','OpenGL')
  % get the base surface
  baseSurface = viewGet(view,'baseSurface');
  % display the surface
  patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', squeeze(img),'facecolor','interp','edgecolor','none','Parent',gui.axis);
  % make sure x direction is normal to make right/right
  set(gui.axis,'XDir','reverse');
  set(gui.axis,'YDir','normal');
  set(gui.axis,'ZDir','normal');
  % set the camera taret to center of surface
  camtarget(gui.axis,mean(baseSurface.vtcs))
  % set the size of the field of view in degrees
  % i.e. 90 would be very wide and 1 would be ver
  % narrow. 9 seems to fit the whole brain nicely
  camva(gui.axis,9);
  % set the view angle
  setMLRViewAngle(view);
end
if verbose>1,disppercent(inf);,end
if verbose>1,disppercent(-inf,'Setting axis');,end
axis(gui.axis,'off');
axis(gui.axis,'image');
if verbose>1,disppercent(inf);,end


% Display colorbar
if verbose>1,disppercent(-inf,'colorbar');,end
cbar = rescale2rgb([1:256],cmap,[1,256],baseGamma);
image(cbar,'Parent',gui.colorbar);
set(gui.colorbar,'YTick',[]);
set(gui.colorbar,'XTick',[1 64 128 192 256]);
set(gui.colorbar,'XTicklabel',num2str(linspace(cbarRange(1),cbarRange(2),5)',3));
if verbose>1,disppercent(inf);,end

% Display the ROIs
if baseType <= 1
  drawnow;
end
nROIs = viewGet(view,'numberOfROIs');
if nROIs
%  if baseType <= 1
    roi = displayROIs(view,slice,sliceIndex,rotate,baseNum,base.coordsHomogeneous,base.dims,verbose);
%  end
else
  roi = [];
end

if verbose>1,disppercent(-inf,'rendering');,end
drawnow
if verbose>1,disppercent(inf);,end
if verbose,toc,end

return

%%%%%%%%%%%%%%%%%%%%%%
%%   getBaseSlice   %%
%%%%%%%%%%%%%%%%%%%%%%
function [baseIm,baseCoords,baseCoordsHomogeneous] = ...
  getBaseSlice(view,sliceNum,sliceIndex,rotate,baseNum,baseType)
%
% getBaseSlice: extracts base image and corresponding coordinates

baseCoords = [];
baseCoordsHomogeneous = [];
baseIm = [];

% viewGet
volSize = viewGet(view,'baseDims',baseNum);
baseData = viewGet(view,'baseData',baseNum);

% get the crop type
% for regular images we want loose, for
% flat maps, we want crop
if baseType == 0,cropType = 'loose';else cropType = 'crop';end

if ~isempty(volSize)

  % Generate coordinates with meshgrid
  switch sliceIndex
    case 1
      x = sliceNum * ones(volSize(2),volSize(3));
      [z,y] = meshgrid(1:volSize(3),1:volSize(2));
    case 2
      y = sliceNum * ones(volSize(1),volSize(3));
      [z,x] = meshgrid(1:volSize(3),1:volSize(1));
    case 3
      z = sliceNum * ones(volSize(1),volSize(2));
      [y,x] = meshgrid(1:volSize(2),1:volSize(1));
  end

  % if there are baseCoords associated with this base
  % then use those instead of image coordinates (this
  % is the case for flat files)
  baseCoordMap = viewGet(view,'baseCoordMap');
  if ~isempty(baseCoordMap)
    % only use the baseCoordMap for when the slice
    % in the third dimension (no other view of a 
    % flat map is really valid).
    if sliceIndex == 3
      x = baseCoordMap.coords(:,:,sliceNum,1);
      y = baseCoordMap.coords(:,:,sliceNum,2);
      z = baseCoordMap.coords(:,:,sliceNum,3);
    else
      oneTimeWarning('badSliceIndex',sprintf('(refreshMLRDisplay:getBaseSlice) Trying to display a flat/surface with the sliceIndex set to %i instead of 3. This is probably because there is something wrong with the Nifti Qforms at your site -- specifically you should check in mrAlign whether your volume displays correctly for when you have click the saggital, coronal and axial buttons. If not, you will need to swap dimensions until they do and then make sure all of your qforms have their dimensions in the same order. Your overlays will not display correctly on this volume.',sliceIndex));
    end
  end

  % Rotate coordinates
  if (rotate ~= 0)
    x = imrotate(x,rotate,'bilinear',cropType);
    y = imrotate(y,rotate,'bilinear',cropType);
    z = imrotate(z,rotate,'bilinear',cropType);
  end

  % Reformat base coordinates
  imageDims = size(x);
  numPixels = prod(imageDims);
  xvec = reshape(x,1,numPixels);
  yvec = reshape(y,1,numPixels);
  zvec = reshape(z,1,numPixels);
  baseCoordsHomogeneous = [xvec; yvec; zvec; ones(1,numPixels)];
  baseCoords = reshape(baseCoordsHomogeneous(1:3,:)',[imageDims 3]);
end

% Extract base image
if ~isempty(baseData)
  switch sliceIndex
    case 1
      baseIm = squeeze(baseData(sliceNum,:,:));
    case 2
      baseIm = squeeze(baseData(:,sliceNum,:));
    case 3
      baseIm = squeeze(baseData(:,:,sliceNum));
  end
  if (rotate ~= 0)
    baseIm = imrotate(baseIm,rotate,'bilinear',cropType);
  end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getROIBaseCoords: extracts ROI coords transformed to the base image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiBaseCoords= getROIBaseCoords(view,baseNum,roiNum);

roiBaseCoords = [];

% viewGet
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
roiCoords = viewGet(view,'roiCoords',roiNum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);
base2roi = viewGet(view,'base2roi',roiNum,baseNum);

if ~isempty(roiCoords) & ~isempty(base2roi)
  % Use xformROI to supersample the coordinates
  roiBaseCoords = round(xformROIcoords(roiCoords,inv(base2roi),roiVoxelSize,baseVoxelSize));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getROIImageCoords: get the roiBaseCoords as x,y and s positions for this
% particular view of the volume (i.e these are coorinates that can
% be plotted on the image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s]= getROIImageCoords(view,roiBaseCoords,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims);

x=[];y=[];s=[];

if ~isempty(roiBaseCoords)
  % make sure that baseCoords are rounded (they may not be
  % if we are working with a baseCoordMap'd flat map
  baseCoordsHomogeneous = round(baseCoordsHomogeneous);

  % base dims
  baseDims = viewGet(view,'baseDims');
  % get the mapping between the image plane and the
  % actual voxel numbers. This will be non-empty for flat surfaces
  baseCoordMap = viewGet(view,'baseCoordMap',baseNum);
  if isfield(baseCoordMap,'dims')
    baseDims = baseCoordMap.dims;
  end

  % get base and roi coordinates linearly (this is so that
  % we don't have to intersect "by rows" which is slower
  % than working with linear indexes by about a factor of 3 -j.
  if isempty(baseCoordMap)
    inplaneIndexes = setdiff(1:3,sliceIndex);
    baseCoordsLinear = mrSub2ind(baseDims(inplaneIndexes),baseCoordsHomogeneous(inplaneIndexes(1),:),baseCoordsHomogeneous(inplaneIndexes(2),:));
    roiBaseCoordsLinear = mrSub2ind(baseDims(inplaneIndexes),roiBaseCoords(inplaneIndexes(1),:),roiBaseCoords(inplaneIndexes(2),:));
    % we use ismember here since it will keep duplicates.
    % in this case we have duplicate roi coordinates that will
    % be found (since we are ignoring which slice we are on-so
    % that we can calculate all the base coordinates irregardless
    % of which slice it is on.
    [roiIndices baseIndices] = ismember(roiBaseCoordsLinear,baseCoordsLinear);
    roiIndices = find(roiIndices);
    baseIndices = baseIndices(roiIndices);
  else
    % for flat maps, we have only one slice and all coordinates
    % are important to match
    baseCoordsLinear = mrSub2ind(baseDims,baseCoordsHomogeneous(1,:),baseCoordsHomogeneous(2,:),baseCoordsHomogeneous(3,:));
    roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
    % we use ismember here since it will keep duplicates.
    % in this case the base may have multiple coordinates, but
    % the roi will have unique coordinates, so we switch
    % the order of arguments from above
    [baseIndices roiIndices] = ismember(baseCoordsLinear,roiBaseCoordsLinear);
    baseIndices = find(baseIndices);
    roiIndices = roiIndices(baseIndices);
  end
  % now transform the coordinates that exist in this
  % base into image coordinates.
  [x,y] = ind2sub(imageDims,baseIndices);
  s = roiBaseCoords(sliceIndex,roiIndices);
end


%%%%%%%%%%%%%%%%%%%%%
%%   displayROIs   %%
%%%%%%%%%%%%%%%%%%%%%
function roi = displayROIs(view,sliceNum,sliceIndex,rotate,baseNum,baseCoordsHomogeneous,imageDims,verbose);
%
% displayROIs: draws the ROIs in the current slice.
%
% viewGet(view,'showROIs') can be:
% 'all'
% 'selected'
% 'all perimeter'
% 'selected perimeter'
% 'hide'
roi = {};

selectedROI = viewGet(view,'currentroi');
labelROIs = viewGet(view,'labelROIs');
n = viewGet(view,'numberOfROIs');

% Order in which to draw the ROIs
option = viewGet(view,'showROIs');
switch option
  case{'all','all perimeter'}
    if selectedROI
      order = [1:selectedROI-1,selectedROI+1:n,selectedROI];
    else
      order = 1:n;
    end
  case{'selected','selected perimeter'}
    order = selectedROI;
  otherwise
    return
end

% Loop through ROIs in order
for r = order
  % look in cache for roi
  roiCache = viewGet(view,'ROICache',r);
  % if not found
  if isempty(roiCache)
    disppercent(-inf,sprintf('Computing ROI base coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    % Get ROI coords transformed to the base dimensions
    roi{r}.roiBaseCoords = getROIBaseCoords(view,baseNum,r);
    % save to cache
    view = viewSet(view,'ROICache',roi{r},r);
    disppercent(inf);
  else
    roi{r} = roiCache;
  end
end

% get figure
fig = viewGet(view,'figNum');
gui = guidata(fig);

% see if this is a flat
baseType = viewGet(view,'baseType',baseNum);

% Draw it
if verbose>1,disppercent(-inf,'Drawing ROI');,end
lineWidth = 0.5;
w = 0.5;

% get which color to draw the selected ROI in
selectedROIColor = mrGetPref('selectedROIColor');
if isempty(selectedROIColor),selectedROIColor = [1 1 1];end

for r = order
  if (r == selectedROI)
    % Selected ROI: set color
    if strcmp(selectedROIColor,'none')
      color = viewGet(view,'roicolor',r);
    else
      color = selectedROIColor;
    end
  else
    % Non-selected ROI, get roi color
    color = viewGet(view,'roicolor',r);
  end
  % then transform them to the image coordinates
  % see if they are in the cached version
  % of the roi. We store here by baseName and
  % sliceIndex. Many different base anatomies
  % may use the same cached roi (since they
  % might have the same xform--i.e. for flat images).
  % but then each base and each sliceIndex of the
  % base needs to have its imageCoords calculated
  baseName = fixBadChars(viewGet(view,'baseName',baseNum));
  % for bases with a baseCoord, need to recalculate roi
  % image coords if the cortical depth has changed
  if viewGet(view,'baseType',baseNum)
    % note we only keep cortical depth to a precision of 1 decimal place
    corticalDepth = sprintf('%0.1g',viewGet(view,'corticalDepth'));
    baseName = sprintf('%s%s',baseName,fixBadChars(corticalDepth,{'.','_'}));
  end
  if ((~isfield(roi{r},baseName)) || ...
      (length(roi{r}.(baseName)) < sliceIndex) || ...
      (isempty(roi{r}.(baseName){sliceIndex})))
    disppercent(-inf,sprintf('Computing ROI image coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    [x y s] = getROIImageCoords(view,roi{r}.roiBaseCoords,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims);
    % keep the coordinates
    roi{r}.(baseName){sliceIndex}.x = x;
    roi{r}.(baseName){sliceIndex}.y = y;
    roi{r}.(baseName){sliceIndex}.s = s;
    disppercent(inf);
    % save in cache
    view = viewSet(view,'ROICache',roi{r},r);
  else
    % just get the coordinates out of the cached ones
    x = roi{r}.(baseName){sliceIndex}.x;
    y = roi{r}.(baseName){sliceIndex}.y;
    s = roi{r}.(baseName){sliceIndex}.s;
  end
  % If it's a 'text' color, translate it...
  roi{r}.color = color2RGB(color);
  % decide whether we are drawing perimeters or not
  doPerimeter = ismember(option,{'all perimeter','selected perimeter'});
  if baseType == 2
    baseSurface = viewGet(view,'baseSurface');
    if 0 %%doPerimeter
      disppercent(-inf,'(refreshMLRDisplay) Computing perimeter');
      baseCoordMap = viewGet(view,'baseCoordMap');
      newy = [];
      for i = 1:length(y)
	% find all the triangles that this vertex belongs to
	[row col] = find(ismember(baseCoordMap.tris,y(i)));
	% get all the neighboring vertices
	neighboringVertices = baseCoordMap.tris(row,:);
	neighboringVertices = setdiff(neighboringVertices(:),y(i));
	% if there are any neighboring vertices that are 
	% not in he roi then this vertex is an edge
	numNeighbors(i) = length(neighboringVertices);
	numROINeighbors(i) = sum(ismember(neighboringVertices,y));
	if numNeighbors(i) ~= numROINeighbors(i)
	  newy = union(newy,baseCoordMap.tris(row(1),:));
	end
	disppercent(i/length(y));
      end
      disppercent(inf);
      disp(sprintf('%i/%i edges',length(newy),length(y)));
      y = newy;
    end
    % display the surface
    roiColors = zeros(size(baseSurface.vtcs));
    roiColors(:) = nan;
    roiColors(y,1) = roi{r}.color(1);
    roiColors(y,2) = roi{r}.color(2);
    roiColors(y,3) = roi{r}.color(3);
    patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', roiColors,'facecolor','interp','edgecolor','none','FaceAlpha',0.4,'Parent',gui.axis);
    continue
  end
  % get image coords for this slice
  if baseType == 0
    % for regular volumes, only get coordinates that match slice
    x = x(s==sliceNum);y = y(s==sliceNum);
  else
    % for flat maps, all coordinates are possible
    % so we just keep all the slices. Note that if
    % we ever make flat maps with multiple slices, then
    % we will have to fix this up here to take that
    % into account.
    sliceNum = 1;
  end
  if ~isempty(x) & ~isempty(y)
    % removing the loop from this algorithm speeds things
    % up dramatically. Calculating the perimeter before
    % for a medium size ROI on the flat could easily take 30
    % seconds on my intel mac. With this verison it takes
    % less than 30ms. yeah :-)...-jg.
    % first get positions of all vertical lines
    % this includes one on either side of each voxel
    % and then make that into a linear coordinate
    % and sort them. Note the +1 on imageDims
    % is to allow for voxels that are at the edge of the image
    % also notice the flip of the imageDims and switching
    % of y and x. This is so that we can get consecutive line
    % segments (see below).
    vlines = sort(sub2ind(fliplr(imageDims)+1,[y y],[x x+1]));
    % now we look for any lines that are duplicates
    % that means they belong to two voxels, so that
    % they should not be drawn. we look for duplicates
    % as ones in which the voxel number is the same
    % as the next one in the list.
    duplicates = diff(vlines)==0;
    % and add the last one in.
    duplicates = [duplicates 0];
    if doPerimeter
      % make sure to score both copies as duplicates
      % only necessary for perimeter drawing, for
      % drawing all boundaries, we need to keep
      % at least one
      duplicates(find(duplicates)+1) = 1;
    end
    % now get everything that is not a duplicate
    vlines = vlines(~duplicates);
    % now we look for consecutive line segments
    % so that we can just draw a single line. This
    % speeds things up because each line segment
    % that matlab draws is an independent child
    % and reducing the number of things it has
    % to keep track of helps things out a lot.
    % so, first we look for the start of a line.
    % this is any index which the difference from
    % its neighbor is greater than one (i.e.
    % the next vertical line being drawn does not
    % happen in the next voxel).
    vlinesBottom = vlines(diff([vlines inf])~=1);
    % now we flip everything and look for the other
    % end of the line in the same way.
    vlinesflip = fliplr(vlines);
    vlinesTop = fliplr(vlinesflip(diff([vlinesflip inf])~=-1));
    % now convert the top and bottom coordinates back
    % to image coordinates
    [vty vtx] = ind2sub(fliplr(imageDims)+1,vlinesTop);
    [vby vbx] = ind2sub(fliplr(imageDims)+1,vlinesBottom);
    % now do the same for the horizontal lines
    hlines = sort(sub2ind(imageDims+1,[x x],[y y+1]));
    duplicates = diff(hlines)==0;
    duplicates = [duplicates 0];
    if doPerimeter
      duplicates(find(duplicates)+1) = 1;
    end
    hlines = hlines(~duplicates);
    hlinesRight = hlines(diff([hlines inf])~=1);
    hlinesflip = fliplr(hlines);
    hlinesLeft = fliplr(hlinesflip(diff([hlinesflip inf])~=-1));
    [hlx hly] = ind2sub(imageDims+1,hlinesLeft);
    [hrx hry] = ind2sub(imageDims+1,hlinesRight);
    % and make them into lines (draw -0.5 and +0.5 so
    % that we draw around the pixel not through the center
    % and note that x/y are flipped
    roi{r}.lines.x = [vty-0.5 hly-0.5;vby+0.5 hry-0.5];
    roi{r}.lines.y = [vtx-0.5 hlx-0.5;vbx-0.5 hrx+0.5];
    % save to cache (since other functions like mrPrint need this
    view = viewSet(view,'ROICache',roi{r},r);
    % now render those lines
    line(roi{r}.lines.x,roi{r}.lines.y,'Color',roi{r}.color,'LineWidth',lineWidth,'Parent',gui.axis);
  else
    roi{r}.lines.x = [];
    roi{r}.lines.y = [];
  end
end
if baseType == 2,return,end
% label them, if labelROIs is set
if labelROIs
  for r = order
    % get x, y lines for this roi on this slice (calculated above)
    x = roi{r}.lines.x(~isnan(roi{r}.lines.x));
    y = roi{r}.lines.y(~isnan(roi{r}.lines.y));
    if ~isempty(x) & ~isempty(y)
      % draw roi label text
      h = text(median(x),median(y),viewGet(view,'roiName',r),'Parent',gui.axis);
      % and set properties
      set(h,'Color','w');
      set(h,'Interpreter','None');
      set(h,'EdgeColor',roi{r}.color);
      set(h,'BackgroundColor','k');
      set(h,'FontSize',10);
      set(h,'HorizontalAlignment','center');
    end
  end
end


if verbose>1,disppercent(inf);,end
return;


