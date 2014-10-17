function [img base roi overlays] = refreshMLRDisplay(viewNum)
%	$Id: refreshMLRDisplay.m 2838 2013-08-12 12:52:20Z julien $

mrGlobals
img = [];base = [];roi = [];overlays=[];
% for debugging/performance tests
% set to 0 for no info
% set to 1 for info on caching
% set to 2 for all info
% if you want to comment out all of this, replace the string
% "if verbose" to "%if verbose"
verbose = 0;
if verbose,tic,end

% get current base
view = viewGet(viewNum,'view');
fig = viewGet(view,'figNum');
gui = guidata(fig);
baseNum = viewGet(view,'currentBase');
baseType = viewGet(view,'baseType');

% if no base then clear axis and return
if isempty(baseNum)
  if ~isempty(fig)
    cla(gui.axis,'reset');
    set(fig,'CurrentAxes',gui.colorbar);
    cla(gui.colorbar,'reset');
    axis(gui.axis,'off');
    axis(gui.colorbar,'off');
  end
  return
end

% set pointer to watch
set(fig,'Pointer','watch');drawnow;

% for debugging, clears caches when user holds down alt key
if ~isempty(fig) && any(strcmp(get(fig,'CurrentModifier'),'alt'))
  disp(sprintf('(refreshMLRDisplay) Dumping caches'));
  view = viewSet(view,'roiCache','init');
  view = viewSet(view,'overlayCache','init');
  view = viewSet(view,'baseCache','init');
end

if verbose>1,disppercent(-inf,'Clearing figure');,end
% note: This cla here is VERY important. Otherwise
% we keep drawing over old things on the axis and
% the rendering gets impossibly slow... -j.
cla(gui.axis);
if verbose>1,disppercent(inf);,end

% check if these are inplanes (not flats or surfaces)
% and see if we should draw all three possible views
if (baseType == 0) && isequal(true,mrGetPref('dispAllPlanesOfAnatomy'))
  mlrGuiSet(view,'multiAxis','on');
  curCoords = viewGet(view,'curCoords');
  for iSliceIndex = 1:3
    [view img{iSliceIndex} base roi overlays curSliceBaseCoords] = dispBase(gui.sliceAxis(iSliceIndex),view,baseNum,gui,true,verbose,iSliceIndex,curCoords(iSliceIndex));
    % get slice to display
%    baseRawSlice = viewGet(view,'baseRawSlice',iSliceOrientation,curCoords(iSliceOrientation));
%    imagesc(flipud(baseRawSlice'),'Parent',gui.sliceAxis(iSliceOrientation));
%    colormap(gui.sliceAxis(iSliceOrientation),gray);
%    axis(gui.sliceAxis(iSliceOrientation),'off');
  end
  % now draw the 3D intersection
  baseDims = viewGet(view,'baseDims',baseNum);

  [x y] = meshgrid(1:baseDims(1),1:baseDims(2));
  imgSurface = surf(gui.axis,x,y,ones(size(x,1),size(x,2))*curCoords(3));
  set(imgSurface,'CData',permute(img{3},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);
  hold(gui.axis,'on');
  
  [x z] = meshgrid(1:baseDims(1),1:baseDims(3));
  imgSurface = surf(gui.axis,x,ones(size(x,1),size(x,2))*curCoords(2),z);
  set(imgSurface,'CData',permute(img{2},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);

  [y z] = meshgrid(1:baseDims(2),1:baseDims(3));
  imgSurface = surf(gui.axis,ones(size(y,1),size(y,2))*curCoords(1),y,z);
  set(imgSurface,'CData',permute(img{1},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);

  axis(gui.axis,'off');
else
  % set the gui to display only a single axis
  mlrGuiSet(view,'multiAxis','off');
  % just draw a single base
  [view img base roi overlays curSliceBaseCoords] = dispBase(gui.axis,view,baseNum,gui,true,verbose);
  % keep cursliceBaseCoords
  view = viewSet(view,'cursliceBaseCoords',curSliceBaseCoords);
end

axes(gui.axis);

% draw any other base that has multiDisplay set
% only do this for surfaces for now
if (baseType >= 2) || ((baseType == 0) && isequal(true,mrGetPref('dispAllPlanesOfAnatomy')))
  for iBase = setdiff(1:viewGet(view,'numBase'),baseNum)
    if viewGet(view,'baseType',iBase)>=2
      if viewGet(view,'baseMultiDisplay',iBase)
	dispBase(gui.axis,view,iBase,gui,false,verbose);
      end
    end
  end
end

if verbose>1,disppercent(-inf,'rendering');end
%this is really stupid: the ListboxTop property of listbox controls seems to be updated only when the control is drawn
%In cases where it is more than the number of overlay names in the box
%it has to be changed, but it is necessary to wait until drawnow before changing it  
%otherwise the change is not taken into account
%Even in this case, It still outputs a warning, that has to be disabled
if strcmp(get(gui.overlayPopup,'style'),'listbox')
  warning('off','MATLAB:hg:uicontrol:ListboxTopMustBeWithinStringRange');
  set(gui.overlayPopup,'ListboxTop',min(get(gui.overlayPopup,'ListboxTop'),length(get(gui.overlayPopup,'string'))));
end
%draw the figure
drawnow
if strcmp(get(gui.overlayPopup,'style'),'listbox')
  warning('on','MATLAB:hg:uicontrol:ListboxTopMustBeWithinStringRange');
end
if verbose>1,disppercent(inf);end
if verbose,toc,end

% set pointer back
set(fig,'Pointer','arrow');

%%%%%%%%%%%%%%%%%%
%    dispBase    %
%%%%%%%%%%%%%%%%%%
function [v img base roi overlays curSliceBaseCoords] = dispBase(hAxis,v,baseNum,gui,dispColorbar,verbose,sliceIndex,slice)

% slight hack here - we set the current base to draw
% each base since that is the way refreshMLRDisplay used
% to work - it defaulted to display the curBase - now it
% needs to display any base that is asked for, so we 
% do that by temporarily making whatever base is is that
% we wish to draw the curBase.
curBase = viewGet(v,'curBase');
v = viewSet(v,'curBase',baseNum);

% Get current view and baseNum.
% Get interp preferences.
% Get slice, scan, alpha, rotate, and sliceIndex from the gui.
if verbose>1,disppercent(-inf,'viewGet');,end
% get variables for current base, but only if they are not set in input 
% if the arguments sliceIndex and slice are set it means we are being
% called to do a "mutliAxis" plot -one in which we are plotting each
% different slice plan in a different axis. For these rotate is going
% to always be 0 since we don't allow that option (it's slow and unnecessary
% and also makes getting the mouse coordinates harder). Also, now that
% volumes with a proper nifti get loaded up in the proper orientation
% it should no longer be necessary to actually rotate images
if nargin < 7
  slice = viewGet(v,'curslice');
  rotate = viewGet(v,'rotate');
  sliceIndex = viewGet(v,'baseSliceIndex',baseNum);
else
  rotate = 0;
end
baseType = viewGet(v,'baseType',baseNum);
baseGamma = viewGet(v,'baseGamma',baseNum);
if verbose>1,disppercent(inf);,end

% Compute base coordinates and extract baseIm for the current slice
if verbose,disppercent(-inf,'extract base image');,end
base = viewGet(v,'baseCache',baseNum,slice,sliceIndex,rotate);
base = [];
if isempty(base)
  [base.im,base.coords,base.coordsHomogeneous] = ...
    getBaseSlice(v,slice,sliceIndex,rotate,baseNum,baseType);
  base.dims = size(base.im);

  % Rescale base volume
  if ~isempty(base.im)
    base.cmap = gray(256);
    base.clip = viewGet(v,'baseClip',baseNum);
    base.RGB = rescale2rgb(base.im,base.cmap,base.clip,baseGamma);
  else
    base.RGB = [];
  end
  % save extracted image
  v = viewSet(v,'baseCache',base,baseNum,slice,sliceIndex,rotate);
  if verbose,disppercent(inf);disp('Recomputed base');end
else
  disp(sprintf('(refershMLRDisplay) Found base'));
  if verbose,disppercent(inf);end
end

if size(base.coords,4)>1
  corticalDepth = viewGet(v,'corticalDepth');
  corticalDepthBins = viewGet(v,'corticalDepthBins');
  corticalDepths = 0:1/(corticalDepthBins-1):1;
  slices = corticalDepths>=corticalDepth(1)-eps & corticalDepths<=corticalDepth(end)+eps; %here I added eps to account for round-off erros
  curSliceBaseCoords = mean(base.coords(:,:,:,slices),4);
else
  curSliceBaseCoords = base.coords;
end

% Extract overlays images and overlays coords, and alphaMap
% right now combinations of overlays are cached as is
% but it would make more sense to cache them separately
% because actual blending occurs after they're separately computed
if verbose,disppercent(-inf,'extract overlays images');end
overlays = viewGet(v,'overlayCache',baseNum,slice,sliceIndex,rotate);
if isempty(overlays)
  % get the transform from the base to the scan
  base2scan = viewGet(v,'base2scan',[],[],baseNum);
  % compute the overlays
  overlays = computeOverlay(v,base2scan,base.coordsHomogeneous,base.dims);
  % save in cache
  v = viewSet(v,'overlayCache',overlays,baseNum,slice,sliceIndex,rotate);
  if verbose,disppercent(inf);disp('Recomputed overlays');end
else
  if verbose,disppercent(inf);end
end
v = viewSet(v,'cursliceOverlayCoords',overlays.coords);

% Combine base and overlays
if verbose>1,disppercent(-inf,'combine base and overlays');,end
if ~isempty(base.RGB) & ~isempty(overlays.RGB)
  switch(mrGetPref('colorBlending'))
    case 'Alpha blend'
      % alpha blending (non-commutative 'over' operator, each added overlay is another layer on top)
      % since commutative, depends on order
      img = base.RGB;
      for iOverlay = 1:size(overlays.RGB,4)
        img = overlays.alphaMaps(:,:,:,iOverlay).*overlays.RGB(:,:,:,iOverlay)+(1-overlays.alphaMaps(:,:,:,iOverlay)).*img;
      end
  
    case 'Additive'
      %additive method (commutative: colors are blended in additive manner and then added as a layer on top of the base)
      % 1) additively multiply colors weighted by their alpha
      RGB = overlays.alphaMaps.*overlays.RGB;
      % 2) add pre-multipliedcolormaps (and alpha channels)
      img = RGB(:,:,:,1);
      alpha = overlays.alphaMaps(:,:,:,1);
      for iOverlay = 2:size(overlays.RGB,4)
        img = img.*(1-RGB(:,:,:,iOverlay))+ RGB(:,:,:,iOverlay);
        alpha = alpha .* (1-overlays.alphaMaps(:,:,:,iOverlay))+overlays.alphaMaps(:,:,:,iOverlay);
      end
      % 2) overlay result on base  using the additively computed alpha (but not for the blended overlays because pre-multiplied)
       img = img+(1-alpha).*base.RGB;
  end
  cmap = overlays.cmap;
  cbarRange = overlays.range;
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
  v = viewSet(v,'curBase',curBase);
  return
end

% if no figure, then just return, this is being called just to create the overlays
fig = viewGet(v,'figNum');
if isempty(fig)
  nROIs = viewGet(v,'numberOfROIs');
  if nROIs
    %  if baseType <= 1
    roi = displayROIs(v,slice,sliceIndex,rotate,baseNum,base.coordsHomogeneous,base.dims,verbose);
    %  end
  else
    roi = [];
  end
  v = viewSet(v,'curBase',curBase);
  return;
end 

% Display the image
if verbose>1,disppercent(-inf,'Displaying image');,end
if baseType <= 1
  % set the renderer to painters (this seems
  % to avoid some weird gliches in the OpenGL
  % renderer. It also appears about 20ms or so
  % faster for displaying images as opposed
  % to the 3D surfaces.
  % 
  % Getting rid of this - to try to rotate the images
  % using view - which apparently cannot be done
  % with 2D images. This way we can display
  % all the images in the correct orientation
  % without having to manually rotate 270 degrees
  %  set(fig,'Renderer','painters')
  %  image(img,'Parent',hAxis);
  set(fig,'Renderer','OpenGL');
  imgSurface = surf(hAxis,zeros(size(img,1),size(img,2)));
  set(imgSurface,'CData',img,'FaceColor','texturemap','EdgeAlpha',0);
  view(hAxis,-90,90);
  set(hAxis,'yDir','reverse');
else
  % set the renderer to OpenGL, this makes rendering
  % *much* faster -- from about 30 seconds to 30ms
  % it is also marginally faster than the zbuffer
  % renderer (order of 5-10 ms)
  set(fig,'Renderer','OpenGL')
  % get the base surface
  baseSurface = viewGet(v,'baseSurface');
  % display the surface
  patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', squeeze(img),'facecolor','interp','edgecolor','none','Parent',hAxis);
  % make sure x direction is normal to make right/right
  set(hAxis,'XDir','reverse');
  set(hAxis,'YDir','normal');
  set(hAxis,'ZDir','normal');
  % set the camera taret to center of surface
  camtarget(hAxis,mean(baseSurface.vtcs))
  % set the size of the field of view in degrees
  % i.e. 90 would be very wide and 1 would be ver
  % narrow. 9 seems to fit the whole brain nicely
  camva(hAxis,9);
  % set the view angle
  setMLRViewAngle(v);
end
if verbose>1,disppercent(inf);,end
if verbose>1,disppercent(-inf,'Setting axis');,end
axis(hAxis,'off');
axis(hAxis,'image');
if verbose>1,disppercent(inf);,end

% Display colorbar
if dispColorbar
  if verbose>1,disppercent(-inf,'colorbar');,end
  cbar = permute(NaN(size(cmap)),[3 1 2]);
  for iOverlay = 1:size(cmap,3)
    cbar(iOverlay,:,:) = rescale2rgb(1:256,cmap(:,:,iOverlay),[1,256],1); 
  end
  image(cbar,'Parent',gui.colorbar);
  if size(cbar,1)==1
    set(gui.colorbar,'YTick',[]);
    set(gui.colorbar,'XTick',[1 64 128 192 256]);
    set(gui.colorbar,'XTicklabel',num2str(linspace(cbarRange(1),cbarRange(2),5)',3));
    if isfield(gui,'colorbarRightBorder')
      set(gui.colorbarRightBorder,'YTick',[]);
    end
  else 
    set(gui.colorbar,'XTick',[]);
    set(gui.colorbar,'YTick',(1:size(cbar,1)));
    set(gui.colorbar,'YTickLabel',cbarRange(:,1));
    if isfield(gui,'colorbarRightBorder')
      set(gui.colorbarRightBorder,'Ylim',[.5 size(cbar,1)+.5],'YTick',(1:size(cbar,1)));
      set(gui.colorbarRightBorder,'YTickLabel',flipud(cbarRange(:,2)));
    end
  end
  if verbose>1,disppercent(inf);,end
end

% Display the ROIs
nROIs = viewGet(v,'numberOfROIs');
if nROIs
  %  if baseType <= 1
  roi = displayROIs(v,hAxis,slice,sliceIndex,rotate,baseNum,base.coordsHomogeneous,base.dims,verbose);
  %  end
else
  roi = [];
end

% reset curBase
v = viewSet(v,'curBase',curBase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getROIBaseCoords: extracts ROI coords transformed to the base image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiBaseCoords= getROIBaseCoords(view,baseNum,roiNum)

roiBaseCoords = [];

% viewGet
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
roiCoords = viewGet(view,'roiCoords',roiNum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);
base2roi = viewGet(view,'base2roi',roiNum,baseNum);

if ~isempty(roiCoords) && ~isempty(base2roi)
  % Use xformROI to supersample the coordinates
  roiBaseCoords = round(xformROIcoords(roiCoords,inv(base2roi),roiVoxelSize,baseVoxelSize));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getROIImageCoords: get the roiBaseCoords as x,y and s positions for this
% particular view of the volume (i.e these are coordinates that can
% be plotted on the image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s]= getROIImageCoords(view,roiBaseCoords,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims)

x=[];y=[];s=[];

if ~isempty(roiBaseCoords)

  % base dims
  % get the mapping between the image plane and the
  % actual voxel numbers. This will be non-empty for flat surfaces
  baseCoordMap = viewGet(view,'baseCoordMap',baseNum);
  if isfield(baseCoordMap,'dims')
    baseDims = baseCoordMap.dims;
  else
    baseDims = viewGet(view,'baseDims');
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
    s = roiBaseCoords(sliceIndex,roiIndices); %s will only be used for volume, not flat maps/surfaces
  else
    % for flat maps, we have only one slice and all coordinates
    % are important to match
    
    if size(baseCoordsHomogeneous,3)>1%if it is a flat map with more than one depth
      corticalDepth = viewGet(view,'corticalDepth');
      corticalDepthBins = viewGet(view,'corticalDepthBins');
      corticalDepths = 0:1/(corticalDepthBins-1):1;
      slices = corticalDepths>=corticalDepth(1)-eps & corticalDepths<=corticalDepth(end)+eps; %here I added eps to account for round-off erros
      nDepths = nnz(slices);
      baseCoordsHomogeneous = reshape(baseCoordsHomogeneous(:,:,slices),[4 prod(imageDims)*nDepths]);
      
    else
      nDepths=1;
    end
    
    % make sure that baseCoords are rounded (they may not be
    % if we are working with a baseCoordMap's flat map
    baseCoordsHomogeneous = round(baseCoordsHomogeneous);

    baseCoordsLinear = mrSub2ind(baseDims,baseCoordsHomogeneous(1,:),baseCoordsHomogeneous(2,:),baseCoordsHomogeneous(3,:));
    roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
    % we use ismember here since it will keep duplicates.
    % in this case the base may have multiple coordinates, but
    % the roi will have unique coordinates, so we switch
    % the order of arguments from above
    [baseIndices roiIndices] = ismember(baseCoordsLinear,roiBaseCoordsLinear);
    
    if nDepths>1%if it is a flat map with more than one depth
      %only keep base voxels that are in a given propotion of the depth levels (given by mrPref roiCorticalDepthDisplayRatio)
      baseIndices=(sum(reshape(baseIndices,[prod(imageDims) nDepths]),2)>=nDepths*mrGetPref('roiCorticalDepthDisplayRatio'))';
      %Since there are several depths, there can be several base roi voxels at each image voxel
      % roiIndices=reshape(roiIndices,[prod(imageDims) nDepths]);
      %but which voxel is displayed is actually not used later because we can only display one depth at a time
      %(roiIndices used to go into s, which was not used for flat maps)
      %so let's just forget about roiIndices 
      %(this also avoid having to think about at what deptht o take the base coords)
    end
    baseIndices = find(baseIndices);
%     roiIndices = roiIndices(baseIndices);
    s=[];%s will only be used for volume, not flat maps/surfaces
  end
  % now transform the coordinates that exist in this
  % base into volume coordinates.
  [x,y] = ind2sub(imageDims,baseIndices);
end


%%%%%%%%%%%%%%%%%%%%%
%%   displayROIs   %%
%%%%%%%%%%%%%%%%%%%%%
function roi = displayROIs(view,hAxis,sliceNum,sliceIndex,rotate,baseNum,baseCoordsHomogeneous,imageDims,verbose);
%
% displayROIs: draws the ROIs in the current slice.
%
% viewGet(view,'showROIs') can be:
% 'all'
% 'selected'
% 'all perimeter'
% 'selected perimeter'
% 'group'
% 'group perimeter'
% 'hide'
roi = {};

selectedROI = viewGet(view,'currentroi');
labelROIs = viewGet(view,'labelROIs');

% Order in which to draw the ROIs
order = viewGet(view,'visibleROIs');
option = viewGet(view,'showROIs');

% Loop through ROIs in order
for r = order
  % look in cache for roi
  roiCache = viewGet(view,'ROICache',r);
  % if not found
  if isempty(roiCache)
    if verbose
      disppercent(-inf,sprintf('Computing ROI base coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    end
    % Get ROI coords transformed to the base dimensions
    roi{r}.roiBaseCoords = getROIBaseCoords(view,baseNum,r);
    % save to cache
    view = viewSet(view,'ROICache',roi{r},r);
    if verbose,disppercent(inf);end
  else
    roi{r} = roiCache;
  end
end

% get figure
fig = viewGet(view,'figNum');
if ~isempty(fig)
  gui = guidata(fig);
end
% see if this is a flat
baseType = viewGet(view,'baseType',baseNum);

% Draw it
if verbose>1,disppercent(-inf,'Drawing ROI');,end

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
    if verbose
      disppercent(-inf,sprintf('Computing ROI image coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    end
    [x y s] = getROIImageCoords(view,roi{r}.roiBaseCoords,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims);
    % keep the coordinates
    roi{r}.(baseName){sliceIndex}.x = x;
    roi{r}.(baseName){sliceIndex}.y = y;
    roi{r}.(baseName){sliceIndex}.s = s;
    if verbose, disppercent(inf); end
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
  doPerimeter = ismember(option,{'all perimeter','selected perimeter','group perimeter'});
  if baseType == 2
    baseSurface = viewGet(view,'baseSurface');
    if 0 %%doPerimeter
      if verbose, disppercent(-inf,'(refreshMLRDisplay) Computing perimeter'); end
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
	if verbose, disppercent(i/length(y)); end;
      end
      if verbose, disppercent(-inf); end;
      disp(sprintf('%i/%i edges',length(newy),length(y)));
      y = newy;
    end
    % display the surface
    roiColors = zeros(size(baseSurface.vtcs));
    roiColors(:) = nan;
    roiColors(y,1) = roi{r}.color(1);
    roiColors(y,2) = roi{r}.color(2);
    roiColors(y,3) = roi{r}.color(3);
 		roi{r}.vertices=y;
		roi{r}.overlayImage = roiColors;
    if ~isempty(fig)
      patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', roiColors,'facecolor','interp','edgecolor','none','FaceAlpha',0.4,'Parent',hAxis);
    end
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
    % JB edit: I've now implemented flat maps/surfaces that average overlay values across
    % different depths. But for ROIs, it is not clear how that would translate
    % in terms of showing the 3D ROI across depths onto the 2D display...
    % so I chose to only display ROI voxels that appear in at least half of the averaged depths
    % (see subfunction getROIImageCoords)
    % Therefore, here we still ignore in which slice (at which depth) ROI voxels actually are
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
    line(roi{r}.lines.x,roi{r}.lines.y,'Color',roi{r}.color,'LineWidth',mrGetPref('roiContourWidth'),'Parent',hAxis);
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
      h = text(median(x),median(y),viewGet(view,'roiName',r),'Parent',hAxis);
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


