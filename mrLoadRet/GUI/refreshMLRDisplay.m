function [img base roi] = refreshMLRDisplay(viewNum)
%	$Id$

mrGlobals
img = [];base = [];roi = [];
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
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
  interpMethod = 'linear';
end
interpExtrapVal = NaN;
view = viewGet(viewNum,'view');
scan = viewGet(view,'curscan');
slice = viewGet(view,'curslice');
alpha = viewGet(view,'alpha');
rotate = viewGet(view,'rotate');
baseNum = viewGet(view,'currentBase');
sliceIndex = viewGet(view,'baseSliceIndex',baseNum);
% if no base then clear axis and return
if isempty(baseNum)
  fig = viewGet(view,'figNum');
  gui = guidata(fig);
  mlrGuiSet(view,'basemin',0);
  mlrGuiSet(view,'basemax',0);
  set(fig,'CurrentAxes',gui.axis);
  cla
  set(fig,'CurrentAxes',gui.colorbar);
  cla
  axis off
  return
end
if verbose>1,disppercent(inf);,end

% for debugging, clears caches when user holds down shift key
if (exist('mglGetKeys')==3) &&  mglGetKeys(57)
  view = viewSet(view,'roiCache','init');
  view = viewSet(view,'overlayCache','init');
  view = viewSet(view,'baseCache','init');
end

% Compute base coordinates and extract baseIm for the current slice
if verbose,disppercent(-inf,'extract base image');,end
base = viewGet(view,'baseCache');
if isempty(base)
  [base.im,base.coords,base.coordsHomogeneous] = ...
    getBaseSlice(view,slice,sliceIndex,rotate,baseNum);
  base.dims = size(base.im);

  % Rescale base volume
  if ~isempty(base.im)
    base.cmap = gray(256);
    base.clip = viewGet(view,'baseClip',baseNum);
    base.RGB = rescale2rgb(base.im,base.cmap,base.clip);
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
  curOverlay = viewGet(view,'currentOverlay');
  analysisNum = viewGet(view,'currentAnalysis');
  [overlayImages,overlay.coords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(view,scan,slice,sliceIndex,rotate,...
    baseNum,base.coordsHomogeneous,base.dims,...
    analysisNum,interpMethod,interpExtrapVal);
  if ~isempty(overlayImages)
    overlayIm = overlayImages(:,:,curOverlay);
  else
    overlayIm = [];
  end

  % get the number of the overlay that shoudl be
  % used as an alpha map. This allows the alpha
  % to be set according to another overlay
  alphaOverlay = viewGet(view,'overlayNum',viewGet(view,'alphaOverlay'));

  numOverlays = viewGet(view,'numberofOverlays');
  mask = ones(size(base.im));
  % Loop through overlays, filling in NaNs according to clip values.
  if ~isempty(overlayImages)
    for ov = 1:numOverlays
      im = overlayImages(:,:,ov);
      clip = viewGet(view,'overlayClip',ov);
      % Find pixels that are within clip
      if diff(clip) > 0
        pts = (im >= clip(1) & im <= clip(2));
      else
        pts = (im >= clip(1) | im <= clip(2));
      end
      mask = mask & pts;
      % if this is the alpha overlay then keep it.
      if ov == alphaOverlay
	alphaOverlayImage = im;
      end
    end
  end
  % Finally, make the alphaMap. Normally this is just set
  % to 0 or 1 so that voxels are hard thresholded. If the
  % overlay has an alphaOverlay field set to the name
  % of another overlay, then we will use the values in
  % that overlay to set the alpha
  if isempty(alphaOverlay)
    overlay.alphaMap = repmat(alpha*mask,[1 1 3]);
  else
    % get the range of the alpha overlay
    range = viewGet(view,'overlayRange',alphaOverlay);
    % handle setRangeToMax
    if strcmp(viewGet(view,'overlayCtype',alphaOverlay),'setRangeToMax')
      maxRange = max(clip(1),min(alphaOverlayImage(mask)));
      if ~isempty(maxRange),range(1) = maxRange;end
      minRange = min(max(alphaOverlayImage(mask)),clip(2));
      if ~isempty(minRange),range(2) = minRange;end
    end
    % now compute the alphaOverlay as a number from
    % 0 to 1 of the range
    alphaOverlayImage = alpha*((alphaOverlayImage-range(1))./diff(range));
    alphaOverlayImage(alphaOverlayImage>alpha) = alpha;
    alphaOverlayImage(alphaOverlayImage<0) = 0;
    alphaOverlayImage = alphaOverlayImage.^viewGet(view,'alphaOverlayExponent');
    alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
    overlay.alphaMap = repmat(alphaOverlayImage.*mask,[1 1 3]); 
  end   
  
  % Rescale current overlay.
  if ~isempty(overlayIm)
    overlay.cmap = viewGet(view,'overlayCmap',curOverlay);
    overlay.range = viewGet(view,'overlayRange',curOverlay);
    if strcmp(viewGet(view,'overlayCtype',curOverlay),'setRangeToMax')
      clip = viewGet(view,'overlayClip',curOverlay);
      if ~isempty(overlayIm(mask))
        overlay.range(1) = max(clip(1),min(overlayIm(mask)));
        overlay.range(2) = min(max(overlayIm(mask)),clip(2));
      else
        overlay.range = clip;
      end
    end
    overlay.RGB = rescale2rgb(overlayIm,overlay.cmap,overlay.range);
  else
    overlay.RGB = [];
  end

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

% Display the image
if verbose>1,disppercent(-inf,'Clearing figure');,end
fig = viewGet(view,'figNum');
gui = guidata(fig);
% note: This cla here is VERY important. Otherwise
% we keep drawing over old things on the axis and
% the rendering gets impossibly slow... -j.
cla(gui.axis);
if verbose>1,disppercent(inf);,end
if verbose>1,disppercent(-inf,'Displaying image');,end
image(img,'Parent',gui.axis);
if verbose>1,disppercent(inf);,end
if verbose>1,disppercent(-inf,'Setting axis');,end
axis(gui.axis,'off');
axis(gui.axis,'image');
if verbose>1,disppercent(inf);,end


% Display colorbar
if verbose>1,disppercent(-inf,'colorbar');,end
cbar = rescale2rgb([1:256],cmap,[1,256]);
image(cbar,'Parent',gui.colorbar);
set(gui.colorbar,'YTick',[]);
set(gui.colorbar,'XTick',[1 64 128 192 256]);
set(gui.colorbar,'XTicklabel',num2str(linspace(cbarRange(1),cbarRange(2),5)',3));
if verbose>1,disppercent(inf);,end

% Display the ROIs
drawnow;
nROIs = viewGet(view,'numberOfROIs');
if nROIs
  roi = displayROIs(view,slice,sliceIndex,rotate,baseNum,base.coordsHomogeneous,base.dims,verbose);
else
  roi = [];
end

if verbose>1,disppercent(-inf,'rendering');,end
drawnow
if verbose>1,disppercent(inf);,end
if verbose,toc,end

return


%-------------------------------------------------------------------------
function rgb = rescale2rgb(image,cmap,clip)
%
% function rgb = rescale2rgb(image,cmap,[clipMin,clipMax])
%
% Clips top and bottom tails of image histogram.
% Sets NaNs to the lowest index in the colormap

if ~exist('clip','var')
  % Choose clipping based on histogram
  histThresh = length(image(:))/1000;
  [cnt, val] = hist(image(:),100);
  goodVals = find(cnt>histThresh);
  clipMin = val(min(goodVals));
  clipMax = val(max(goodVals));
  clip = [clipMin,clipMax];
else
  clipMin = clip(1);
  clipMax = clip(2);
end

% Clip
result = image;
result(find(image < clipMin)) = clipMin;
result(find(image > clipMax)) = clipMax;

% Scale
indices = round(255 * (result-clipMin)/(clipMax-clipMin)) + 1;
indices = max(1,min(indices,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(image));
g = zeros(size(image));
b = zeros(size(image));
r(:) = cmap(indices,1);
g(:) = cmap(indices,2);
b(:) = cmap(indices,3);

% Stuff them into rgb
dims = [size(image),3];
rgb = cat(length(dims),r,g,b);


%-------------------------------------------------------------------------
function [baseIm,baseCoords,baseCoordsHomogeneous] = ...
  getBaseSlice(view,sliceNum,sliceIndex,rotate,baseNum)
%
% getBaseSlice: extracts base image and corresponding coordinates

baseCoords = [];
baseCoordsHomogeneous = [];
baseIm = [];

% viewGet
volSize = viewGet(view,'baseDims',baseNum);
baseData = viewGet(view,'baseData',baseNum);

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
    end
  end

  % Rotate coordinates
  if (rotate ~= 0)
    x = imrotate(x,rotate,'nearest','loose');
    y = imrotate(y,rotate,'nearest','loose');
    z = imrotate(z,rotate,'nearest','loose');
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
    baseIm = imrotate(baseIm,rotate,'nearest','loose');
  end
end

return


%-------------------------------------------------------------------------
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
  getOverlaySlice(view,...
  scanNum,sliceNum,sliceIndex,rotate,...
  baseNum,baseCoordsHomogeneous,imageDims,...
  analysisNum,interpMethod,interpExtrapVal);
%
% getOverlaySlice: extracts overlay image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
baseXform = viewGet(view,'baseXform',baseNum);
overlayXform = viewGet(view,'overlayXform',scanNum);
numOverlays = viewGet(view,'numberofoverlays',analysisNum);
interpFnctn = viewGet(view,'overlayInterpFunction',analysisNum);

% Transform base coords corresponding to this slice/image to overlay
% coordinate frame.
% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin.
shiftXform = shiftOriginXform;
if ~isempty(overlayXform) & ~isempty(baseXform) & ~isempty(baseCoordsHomogeneous)
  xform = inv(shiftXform) * inv(overlayXform) * baseXform * shiftXform;
  % Transform coordinates
  overlayCoordsHomogeneous = xform * baseCoordsHomogeneous;
  overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[imageDims 3]);
end

% Extract overlayImages
if ~isempty(interpFnctn)
  overlayImages = feval(interpFnctn,view,scanNum,imageDims,...
    analysisNum,overlayCoords,interpMethod,interpExtrapVal);
else
  overlayImages = zeros([imageDims,numOverlays]);
  for ov = 1:numOverlays
    overlayData = viewGet(view,'overlayData',scanNum,ov,analysisNum);
    if ~isempty(overlayData) & ~isempty(overlayCoords)
      % Extract the slice
      overlayImages(:,:,ov) = interp3(overlayData,...
        overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
        interpMethod,interpExtrapVal);
    end
  end
end

%-------------------------------------------------------------------------
function overlayImages = corAnalInterp(view,scanNum,...
  imageDims,analysisNum,overlayCoords,interpMethod,interpExtrapVal);
%
% corAnalInterp: special case for corAnal. Need to treat amp/ph as complex
% valued.

% Initialize
numOverlays = viewGet(view,'numberofoverlays',analysisNum);
overlayImages = zeros([imageDims,numOverlays]);

% Interpolate complex values for amp and ph
ampNum = viewGet(view,'overlayNum','amp',analysisNum);
phNum = viewGet(view,'overlayNum','ph',analysisNum);
ampData = viewGet(view,'overlaydata',scanNum,ampNum,analysisNum);
phData = viewGet(view,'overlaydata',scanNum,phNum,analysisNum);
zData = ampData .* exp(i*phData);
if ~isempty(zData) & ~isempty(overlayCoords)
  zInterp = interp3(zData,...
    overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
    interpMethod,interpExtrapVal);
  overlayImages(:,:,ampNum) = abs(zInterp);
  overlayImages(:,:,phNum) = angle(zInterp)+pi;
end

% Loop through other overlays and extract images using normal interpolation
% methods. 
for ov = setdiff(1:numOverlays,[phNum ampNum])
  overlayData = viewGet(view,'overlayData',scanNum,ov,analysisNum);
  if ~isempty(overlayData) & ~isempty(overlayCoords)
    % Extract the slice
    overlayImages(:,:,ov) = interp3(overlayData,...
      overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
      interpMethod,interpExtrapVal);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getROIBaseCoords: extracts ROI coords transformed to the base image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiBaseCoords= getROIBaseCoords(view,baseNum,roiNum);

roiBaseCoords = [];

% viewGet
baseXform = viewGet(view,'baseXform',baseNum);
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
roiCoords = viewGet(view,'roiCoords',roiNum);
roiXform = viewGet(view,'roiXform',roiNum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);

if ~isempty(roiCoords) & ~isempty(roiXform) & ~isempty(baseXform)
  % Use xformROI to supersample the coordinates
  roiBaseCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));
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
  % base into image coordinates
  [x,y] = ind2sub(imageDims,baseIndices);
  s = roiBaseCoords(sliceIndex,roiIndices);
end


%-------------------------------------------------------------------------
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
    % Selected ROI: set color=white
    color = selectedROIColor;
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
  if isstr(color)
    switch (color)
     case {'yellow','y'}, color = [1 1 0];
     case {'magenta','m'}, color = [1 0 1];
     case {'cyan','c'}, color = [0 1 1];
     case {'red','r'}, color = [1 0 0];
     case {'green','g'}, color = [0 1 0];
     case {'blue','b'}, color = [0 0 1];
     case {'white','w'}, color = [1 1 1];
     case {'black','k'}, color = [0 0 0];
     otherwise, color = [1 1 1];
    end % end switch statement
  end
  roi{r}.color = color;
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
    % decide whether we are drawing perimeters or not
    doPerimeter = ismember(option,{'all perimeter','selected perimeter'});
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
    vlines = sort(sub2ind(imageDims+1,[x x+1],[y y]));
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
    % now make back into x,y coordinates
    [vx vy] = ind2sub(imageDims+1,vlines);
    % now do the same for the horizontal lines
    hlines = sort(sub2ind(imageDims+1,[x x],[y y+1]));
    duplicates = diff(hlines)==0;
    duplicates = [duplicates 0];
    if doPerimeter
      duplicates(find(duplicates)+1) = 1;
    end
    hlines = hlines(~duplicates);
    [hx hy] = ind2sub(imageDims+1,hlines);
    % and make them into lines (draw -0.5 and +0.5 so
    % that we draw around the pixel not through the center
    % and note that x/y are flipped
    roi{r}.lines.x = [vy-0.5 hy-0.5;vy+0.5 hy-0.5];
    roi{r}.lines.y = [vx-0.5 hx-0.5;vx-0.5 hx+0.5];
    % save to cache (since other functions like mrPrint need this
    view = viewSet(view,'ROICache',roi{r},r);
    % now render those lines
    line(roi{r}.lines.x,roi{r}.lines.y,'Color',color,'LineWidth',lineWidth,'Parent',gui.axis);
  else
    roi{r}.lines.x = [];
    roi{r}.lines.y = [];
  end
end

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


