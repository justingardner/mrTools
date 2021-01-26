function [img, base, roi, overlays, altBase, v] = refreshMLRDisplay(viewNum)
%	$Id: refreshMLRDisplay.m 2838 2013-08-12 12:52:20Z julien $

mrGlobals
img = [];base = [];roi = [];overlays=[];altBase = [];
% for debugging/performance tests
% set to 0 for no info
% set to 1 for info on caching
% set to 2 for all info
% if you want to comment out all of this, replace the string
% "if verbose" to "%if verbose"
verbose = 0;
if verbose,tic,end

% get current base
v = viewGet(viewNum,'view');
viewNum = viewGet(v,'viewNum');
fig = viewGet(v,'figNum');
if ~isempty(fig)
  gui = guidata(fig);
elseif nargout>0
  gui.axis = [];
else
  return  % if there is no GUI and no output arguments, then there is nothing to do
end
baseNum = viewGet(v,'currentBase');
baseType = viewGet(v,'baseType');

% debug code
%curCoords = viewGet(v,'curCoords');
%disp(sprintf('(refreshMLRDisplay) DEBUG: [%i %i %i]',curCoords(1),curCoords(2),curCoords(3)));

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
  v = viewSet(v,'roiCache','init');
  v = viewSet(v,'overlayCache','init');
  v = viewSet(v,'baseCache','init');
end

if verbose>1,mlrDispPercent(-inf,'Clearing figure');,end
% note: This cla here is VERY important. Otherwise
% we keep drawing over old things on the axis and
% the rendering gets impossibly slow... -j.
if ~isempty(fig), cla(gui.axis);end
if verbose>1,mlrDispPercent(inf);,end

% check if these are inplanes (not flats or surfaces)
% and see if we should draw all three possible views
baseMultiAxis = viewGet(v,'baseMultiAxis');
if (baseType == 0) && (baseMultiAxis>0)
  % set the axis, either to dispaly all 3 views + a 3D
  % or just a 3d
  if baseMultiAxis == 1
    mlrGuiSet(v,'multiAxis',1);
  else
    mlrGuiSet(v,'multiAxis',2);
  end
  % get what coords display
  curCoords = viewGet(v,'curCoords');

  % display the images, but if we are set to 3D just compute
  % the images, so that we can display the 3D.
  for iSliceIndex = 1:3
    % if we are displaying all slice dimensions and the #d
    if baseMultiAxis == 1
      % display each slice index
      [v img{iSliceIndex} base roi{iSliceIndex} overlays curSliceBaseCoords] = dispBase(gui.sliceAxis(iSliceIndex),v,baseNum,gui,iSliceIndex==3,verbose,iSliceIndex,curCoords(iSliceIndex));
    else
      % display each slice index
      [v img{iSliceIndex} base roi{iSliceIndex}] = dispBase([],v,baseNum,gui,iSliceIndex==3,verbose,iSliceIndex,curCoords(iSliceIndex));
    end
  end

  % now draw the 3D intersection
  baseDims = viewGet(v,'baseDims',baseNum);

  roiContourWidth = mrGetPref('roiContourWidth');

  % draw each dimension in turn using images that were returned
  % by each call to dispBase from above
  % 1st dimension (axial)
  cla(gui.axis);
  set(fig,'Renderer','OpenGL');
  [x y] = meshgrid(1:baseDims(1),1:baseDims(2));
  imgSurface = surf(gui.axis,y,x,ones(size(x,1),size(x,2))*curCoords(3));
  set(imgSurface,'CData',permute(img{3},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);
  hold(gui.axis,'on');
  % draw rois
  for iROI = 1:length(roi{3})
    if ~isempty(roi{3}{iROI})
      z = ones(size(roi{3}{iROI}.lines.x))*curCoords(3);
      line(roi{3}{iROI}.lines.x,roi{3}{iROI}.lines.y,z,'Color',roi{3}{iROI}.color,'LineWidth',roiContourWidth,'Parent',gui.axis);
    end
  end

  % 2nd dimension (coronal)
  [x z] = meshgrid(1:baseDims(1),1:baseDims(3));
  imgSurface = surf(gui.axis,ones(size(x,1),size(x,2))*curCoords(2),x,z);
  set(imgSurface,'CData',permute(img{2},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);
  % draw rois
  for iROI = 1:length(roi{2})
    if ~isempty(roi{2}{iROI})
      z = ones(size(roi{2}{iROI}.lines.x))*curCoords(2);
      line(z,roi{2}{iROI}.lines.y,roi{2}{iROI}.lines.x,'Color',roi{2}{iROI}.color,'LineWidth',roiContourWidth,'Parent',gui.axis);
    end
  end

  % 3rd dimension (saggital)
  [y z] = meshgrid(1:baseDims(2),1:baseDims(3));
  imgSurface = surf(gui.axis,y,ones(size(y,1),size(y,2))*curCoords(1),z);
  set(imgSurface,'CData',permute(img{1},[2 1 3]),'FaceColor','texturemap','EdgeAlpha',0);
  % draw rois
  for iROI = 1:length(roi{1})
    if ~isempty(roi{1}{iROI})
      z = ones(size(roi{1}{iROI}.lines.x))*curCoords(1);
      line(roi{1}{iROI}.lines.y,z,roi{1}{iROI}.lines.x,'Color',roi{1}{iROI}.color,'LineWidth',roiContourWidth,'Parent',gui.axis);
    end
  end

  % set the axis
  axis(gui.axis,'off');
  axis(gui.axis,'tight');
  % make sure x direction is normal to make right/right
  set(gui.axis,'XDir','reverse');
  set(gui.axis,'YDir','normal');
  set(gui.axis,'ZDir','normal');
  % set the size of the field of view in degrees
  % i.e. 90 would be very wide and 1 would be ver
  % narrow. 9 seems to fit the whole brain nicely
  % this hard-coded value also appears in dispBase
  camva(gui.axis,9);
  axis(gui.axis,'equal');
  daspect(gui.axis,1./viewGet(v,'baseVoxelSize'))
  setMLRViewAngle(v,gui.axis);
  
else
  % set the gui to display only a single axis
  mlrGuiSet(v,'multiAxis',0);
  % just draw a single base
  [v img base roi overlays curSliceBaseCoords] = dispBase(gui.axis,v,baseNum,gui,true,verbose);
  % keep cursliceBaseCoords
  v = viewSet(v,'cursliceBaseCoords',curSliceBaseCoords);
end

% turn on 3D free rotate if we are just displaying the one 3D axis
if ~isempty(fig) && ~mrInterrogator('isactive',viewNum)
  if (baseType == 2) || (baseMultiAxis == 2)
    mlrSetRotate3d(v,'on');
  else
    mlrSetRotate3d(v,'off');
  end
end

if ~isempty(gui.axis)
  axes(gui.axis);
end

% draw any other base that has multiDisplay set
% do this for surfaces for now or for 3D anatomies
% when we have multiAxis on (since that draws a 3D
% view with all the relevant planes)
if (baseType >= 2) || ((baseType == 0) && (baseMultiAxis>0))
  for iBase = setdiff(1:viewGet(v,'numBase'),baseNum)
    if viewGet(v,'baseType',iBase)>=2
      if viewGet(v,'baseMultiDisplay',iBase)
	% display the base and keep the information about what is drawn for functions like mrPrint
	[v altBase(iBase).img altBase(iBase).base altBase(iBase).roi altBase(iBase).overlays] = dispBase(gui.axis,v,iBase,gui,true,verbose);
      end
    end
  end
end

if ~isempty(gui.axis)
  if (baseType == 0) && (baseMultiAxis>0)
    % set the camera target to center of the volume
    camtarget(gui.axis,baseDims([2 1 3])/2)
  end
end

if verbose>1,mlrDispPercent(-inf,'rendering');end

if ~isempty(fig)
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
  drawnow('update');
  if strcmp(get(gui.overlayPopup,'style'),'listbox')
    warning('on','MATLAB:hg:uicontrol:ListboxTopMustBeWithinStringRange');
  end
end

if verbose>1,mlrDispPercent(inf);end
if verbose,toc,end

% set pointer back
set(fig,'Pointer','arrow');

%%%%%%%%%%%%%%%%%%
%    dispBase    %
%%%%%%%%%%%%%%%%%%
function [v img base roi overlays curSliceBaseCoords] = dispBase(hAxis,v,baseNum,gui,dispColorbar,verbose,sliceIndex,slice)

baseType = viewGet(v,'baseType',baseNum);
baseGamma = viewGet(v,'baseGamma',baseNum);

% Get current view and baseNum.
% Get interp preferences.
% Get slice, scan, alpha, rotate, and sliceIndex from the gui.
if verbose>1,mlrDispPercent(-inf,'viewGet');,end
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
  if baseType < 2
    rotate = viewGet(v,'rotate');
  else
    rotate = 0;
  end
  sliceIndex = viewGet(v,'baseSliceIndex',baseNum);
else
  rotate = 0;
end
if verbose>1,mlrDispPercent(inf);end

fig = viewGet(v,'figNum');
%disp(sprintf('(refreshMLRDIsplay:dispBase) DEBUG: sliceIndex: %i slice: %i',sliceIndex,slice));

% Compute base coordinates and extract baseIm for the current slice
if verbose,mlrDispPercent(-inf,'extract base image');end
base = viewGet(v,'baseCache',baseNum,slice,sliceIndex,rotate);
if isempty(base)
  [base.im,base.coords,base.coordsHomogeneous] = ...
    getBaseSlice(v,slice,sliceIndex,rotate,baseNum,baseType);
  base.dims = size(base.im);

  % Rescale base volume
  if ~isempty(base.im)
    base.cmap = gray(256);
    base.clip = viewGet(v,'baseClip',baseNum);
    base.RGB = rescale2rgb(base.im,base.cmap,base.clip,baseGamma);
    switch lower(mrGetPref('baseNaNsColor'))
      case 'transparent' %display NaNs the same color as figure background
        if ~isempty(fig) 
          backgroundColor = get(fig,'color');
        else
          backgroundColor = [0 0 0];
        end
      case 'black'
        backgroundColor = [0 0 0];
      case 'white'
        backgroundColor = [1 1 1];
    end
    if baseType==1 %make smooth transition beetween figure background and flat map
      alpha = zeros(base.dims(1),base.dims(2));
      alpha(isnan(base.im))=1;
      kernel = gaussianKernel2D(3)   ;     
      alpha = conv2(alpha,ones(4,4),'same');
      alpha = conv2(1-(alpha>0),kernel,'same');
      alpha = repmat(alpha,[1 1 3]);
      mask = repmat(permute(backgroundColor,[1 3 2]),[base.dims(1) base.dims(2) 1]);
      base.RGB = base.RGB.*alpha+(1-alpha).*mask;
      base.RGB=min(1,max(0,base.RGB));
      base.gyrusSulcusBoundary = edge(base.im>0.5+0)&edge(base.im<0.5+0); %we assume that 0.5 represents 
      % the curvature boundary, which should be the case for flat maps made from freesurfer-imported surfaces 
    else
      base.RGB = reshape(base.RGB,base.dims(1)*base.dims(2),3);
      base.RGB(isnan(base.im),:)=repmat(backgroundColor,[nnz(isnan(base.im)) 1]);
      base.RGB = reshape(base.RGB,[base.dims 3]);
    end
  else
    base.RGB = [];
  end
  % save extracted image
  v = viewSet(v,'baseCache',base,baseNum,slice,sliceIndex,rotate);
  if verbose,mlrDispPercent(inf);disp('Recomputed base');end
else
  if verbose,mlrDispPercent(inf);end
end

% for surfaces and flats calculate things based on cortical depth
if size(base.coords,4)>1
  % code allows for averaging across cortical depth
  corticalDepth = viewGet(v,'corticalDepth',baseNum);
  corticalDepthBins = mrGetPref('corticalDepthBins');
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
if verbose,mlrDispPercent(-inf,'extract overlays images');end
overlays = viewGet(v,'overlayCache',baseNum,slice,sliceIndex,rotate);
if isempty(overlays)
  % get the transform from the base to the scan
  base2scan = viewGet(v,'base2scan',[],[],baseNum);
  % compute the overlays
  overlays = computeOverlay(v,base2scan,base.coordsHomogeneous,base.dims,[],baseNum);
  overlays = addBaseOverlays(v,baseNum,overlays);
  % save in cache
  v = viewSet(v,'overlayCache',overlays,baseNum,slice,sliceIndex,rotate);
  if verbose,mlrDispPercent(inf);disp('Recomputed overlays');end
else
  if verbose,mlrDispPercent(inf);end
end
v = viewSet(v,'cursliceOverlayCoords',overlays.coords);

% Combine base and overlays
if verbose>1,mlrDispPercent(-inf,'combine base and overlays');,end
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
       
    case 'Contours'
      if size(overlays.RGB,4)==1
        img=base.RGB; %only display base in the background
        contours=overlays.overlayIm;  %display first overlay as contours 
        contours(~overlays.alphaMaps(:,:,1,1))=NaN;   % use non-zero alpha map as mask
        contourCscale=overlays.colorRange(1,:); 
        contourCmap = overlays.cmap(:,:,1);
      else
        %combine first overlay and base
        img=overlays.alphaMaps(:,:,:,1).*overlays.RGB(:,:,:,1)+(1-overlays.alphaMaps(:,:,:,1)).*base.RGB;
        contours=overlays.overlayIm(:,:,2); %display second overlay as contours
        contours(~overlays.alphaMaps(:,:,1,2))=NaN;   % use non-zero alpha map as mask
        contourCscale=overlays.colorRange(2,:);
        contourCmap = overlays.cmap(:,:,2);
      end
      if size(overlays.RGB,4)>2
        mrWarnDlg('(refreshMLRDisplay) Number of overlays limited to 2 for ''Contours'' display option');
      end
  end
  displayGyrusSulcusBoundary = viewGet(v,'displayGyrusSulcusBoundary');
  if baseType==1 && ~isempty(displayGyrusSulcusBoundary) && displayGyrusSulcusBoundary
      img = reshape(img,[prod(base.dims) 3]);
      img(base.gyrusSulcusBoundary>0,:) = repmat([0 0 0],nnz(base.gyrusSulcusBoundary>0),1);
      img = reshape(img,[base.dims 3]);
  end
  cmap = overlays.cmap;
  cbarRange = overlays.colorRange;
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
if verbose>1,mlrDispPercent(inf);,end

% If no image at this point then return
if ieNotDefined('img')
  return
end
img=double(img); %in case data are stored as singles (later instances of surf/set don't work with singles)


% if no figure, then just return, this is being called just to create the overlays
if isempty(fig)
  nROIs = viewGet(v,'numberOfROIs');
  if nROIs
    %  if baseType <= 1
    roi = displayROIs(v,hAxis,slice,sliceIndex,baseNum,base.coordsHomogeneous,base.dims,rotate,verbose);
    %  end
  else
    roi = [];
  end
  return;
end 

% Display colorbar (this should be done before returning when there is a
% figure, but the slice is not drawn (3D view), otherwise the colorbar
% won't get drwan/updated
if dispColorbar
  displayColorbar(gui,cmap,cbarRange,verbose)
end

% if we are not displaying then, just return    %JB: this seems redundant with the previous, identical, block of code
% after computing rois                          %     Does it ever happen that there is a figure, but the axis is empty? (a figure with just a color bar?)
if isempty(hAxis)
  % just compute the axis (displayROIs will not draw
  % if hAxis is set to empty. We need the ROI x,y for
  % drawing the coordinates into the multiAxis 3D
  nROIs = viewGet(v,'numberOfROIs');
  if nROIs
    roi = displayROIs(v,hAxis,slice,sliceIndex,baseNum,base.coordsHomogeneous,base.dims,rotate,verbose);
  else
    roi = [];
  end
  return
end

% Display the image
if verbose>1,mlrDispPercent(-inf,'Displaying image');,end
if baseType <= 1
  % set the renderer to painters (this seems
  % to avoid some weird gliches in the OpenGL
  % renderer. It also appears about 20ms or so
  % faster for displaying images as opposed
  % to the 3D surfaces.
  % 
  % Just draw with img for regular image,
  if isequal(0,viewGet(v,'baseMultiAxis',baseNum))
    set(fig,'Renderer','painters')
    h = image(img,'Parent',hAxis);
    v = viewSet(v,'baseHandle',h,baseNum);
    if ~isempty(overlays.RGB) && strcmp(mrGetPref('colorBlending'),'Contours') 
      hold(hAxis,'on');
      nContours = 9;
      contourStepSize =diff(contourCscale)/(nContours-1);
      minContourStep = ceil(min(min(contours))/contourStepSize)*contourStepSize;
      maxContourStep = floor(max(max(contours))/contourStepSize)*contourStepSize;
      contour(hAxis,contours,minContourStep:contourStepSize:maxContourStep,'cdatamapping','scaled')
      caxis(hAxis,contourCscale);
      colormap(contourCmap);
    end
    % get rotation matrix for image so that we can fix the data
    % aspect ratio properlybelow
    theta = pi*rotate/180;
    m = [cos(theta) sin(theta);-sin(theta) cos(theta)];
  else  
    % for multi axis, we want to have them rotated by 270
    % which apparently cannot be done using view with 2D images. 
    % We want 270 rotation since This way we can display
    % all the images in the correct orientation
    % without having to manually rotate 270 degrees
    % so we display them as textures on planer surfaces and
    % display those rotated correctly using view
    set(fig,'Renderer','OpenGL');
    imgSurface = surf(hAxis,zeros(size(img,1),size(img,2)));
    set(imgSurface,'CData',img,'FaceColor','texturemap','EdgeAlpha',0);
    v = viewSet(v,'baseHandle',imgSurface,baseNum);
    view(hAxis,-90,90);
    set(hAxis,'yDir','reverse');
    % set the rotation matrix for the data aspect ratio to identity
    % since we do not do any rotation here.
    m = eye(2);
  end
  % for anatomies to respect their rotation aspect ratios
  % this code fixes things - but should not do this for flats
  if baseType == 0
    % get data aspect ratio from voxel sizes
    baseVoxelSize = viewGet(v,'baseVoxelSize',baseNum);
    baseVoxelSize = baseVoxelSize(setdiff([1 2 3],sliceIndex));
    % rotate axis according to the rotation of the image
    xAxisRotated = m*[baseVoxelSize(1) 0]';
    yAxisRotated = m*[0 baseVoxelSize(2)]';
    % get aspect ratio in rotated frame. Still not quite sure
    % this works but the idea is to take the unit vectors for the
    % x and y axis rotated (from above) and see which is larger the
    % spread between them, or the distance from 0. And take that
    % as how much you have to scale each axis by.
    aspectRatio = max(abs(max(xAxisRotated,yAxisRotated)),abs(xAxisRotated-yAxisRotated));
    % now set the data aspect ratio so that images showup
    % with the aspect ratio appropriately for the voxel size
    daspect(hAxis,[aspectRatio' 1]);
  else
    axis(hAxis,'image'); %set aspect ratio to 1:1 for flat maps
  end
else
  % set the renderer to OpenGL, this makes rendering
  % *much* faster -- from about 30 seconds to 30ms
  % it is also marginally faster than the zbuffer
  % renderer (order of 5-10 ms)
  set(fig,'Renderer','OpenGL')
  % get the base surface
  baseSurface = getBaseSurface(v,baseNum); %get baseSurface coordinates, convert to different base space if necessary
  % get alpha
  baseAlpha = viewGet(v,'baseAlpha',baseNum);
  % display the surface
  hSurface = patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', squeeze(img),'facecolor','interp','edgecolor','none','Parent',hAxis,'FaceAlpha',baseAlpha);
  % keep the surface handle
  v = viewSet(v,'baseHandle',hSurface,baseNum);
  % make sure x direction is normal to make right/right
  set(hAxis,'XDir','reverse');
  set(hAxis,'YDir','normal');
  set(hAxis,'ZDir','normal');
  % set the camera target to center of surface (but only if curBase)
  if isequal(viewGet(v,'curBase'),baseNum)
    camtarget(hAxis,mean(baseSurface.vtcs))
  end
  % set the size of the field of view in degrees
  % i.e. 90 would be very wide and 1 would be ver
  % narrow. 9 seems to fit the whole brain nicely
  camva(hAxis,9);
  % set the view angle
  setMLRViewAngle(v);
  if baseNum == viewGet(v,'currentBase') % set aspect ratio to 1:1:1, but only if this is the current base
    axis(hAxis,'image');                 % (if not, the aspect has already been set for the current base)
  end
end
if verbose>1,mlrDispPercent(inf);,end
if verbose>1,mlrDispPercent(-inf,'Setting axis');,end
axis(hAxis,'off');
if verbose>1,mlrDispPercent(inf);,end

% Display ROIs
nROIs = viewGet(v,'numberOfROIs');
if nROIs
  [roi,v] = displayROIs(v,hAxis,slice,sliceIndex,baseNum,base.coordsHomogeneous,base.dims,rotate,verbose);
else
  roi = [];
end

%Display Surface contours on volume
if ~baseType
  v = displaySurfaceOnVolume(v,hAxis,slice,sliceIndex,rotate,base.dims,slice,baseNum);
end

% if the interrogator is on, then we need to reinitialize it
if baseType==2 && ~verLessThan('matlab','8.4') && mrInterrogator('isactive',viewGet(v,'viewNum'))
  mrInterrogator('init',viewGet(v,'viewNum'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displayColorbar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayColorbar(gui,cmap,cbarRange,verbose)
  if isempty(cmap)
    set(gui.colorbar,'Visible','off');
  else
    set(gui.colorbar,'Visible','on');
    
    if verbose>1,mlrDispPercent(-inf,'colorbar');,end
    cbar = permute(NaN(size(cmap)),[3 1 2]);
    for iOverlay = 1:size(cmap,3)
      cbar(iOverlay,:,:) = rescale2rgb(1:size(cmap,1),cmap(:,:,iOverlay),[1,size(cmap,1)],1); 
    end
    image(cbar,'Parent',gui.colorbar);
    if size(cbar,1)==1
      set(gui.colorbar,'YTick',[]);
      set(gui.colorbar,'XTick',linspace(0.5,size(cmap,1)+0.5,5));
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
    if verbose>1,mlrDispPercent(inf);,end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getBaseSurface: gets baseSurface using viewGet, but converts vertices coordinates to
% the current base space (in case surface is overlaid onto different base in 3D view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function baseSurface = getBaseSurface(v,baseNum)
  
if baseNum == viewGet(v,'currentBase')
  baseSurface = viewGet(v,'baseSurface',baseNum);
else
  % if this is not the current base, this means we're using mutlibase
  % we need to get the actual surface coordinates
  % and xform them so that it will display in the same space
  % as the current one
  baseSurface = viewGet(v,'baseSurfaceCoords',baseNum);
  % get xform that xforms coordinates from this base to
  % the current base
  base2base = inv(viewGet(v,'base2base',baseNum));
  if isempty(base2base)
    disp(sprintf('(refreshMLRDisplay:getBaseSurface) Could not compute base2base xform'));
    return
  end
  % remove some sigfigs so that the isequal check works
  if ~isequal(round(base2base*100000)/100000,eye(4))
    % homogenous coordinates
    baseSurface.vtcs(:,4) = 1;
    swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];
    % xform to current base coordinates
    baseSurface.vtcs = (swapXY * base2base * swapXY * baseSurface.vtcs')';
    % and remove the homogenous 1
    baseSurface.vtcs = baseSurface.vtcs(:,1:3);
  end
end

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
      corticalDepth = viewGet(view,'corticalDepth',baseNum);
      corticalDepthBins = mrGetPref('corticalDepthBins');
      corticalDepths = 0:1/(corticalDepthBins-1):1;
      slices = corticalDepths>=corticalDepth(1)-eps & corticalDepths<=corticalDepth(end)+eps; %here I added eps to account for round-off errors
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
      %(this also avoids having to think about at what depth to take the base coords)
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
%    displayROIs    %
%%%%%%%%%%%%%%%%%%%%%
function [roi,view] = displayROIs(view,hAxis,sliceNum,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims,rotate,verbose)
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
  roiCache = viewGet(view,'ROICache',r,baseNum,rotate);
  % if not found
  if isempty(roiCache)
    if verbose
      mlrDispPercent(-inf,sprintf('Computing ROI base coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    end
    % Get ROI coords transformed to the base dimensions
    roi{r}.roiBaseCoords = getROIBaseCoords(view,baseNum,r);
    % save to cache
    view = viewSet(view,'ROICache',roi{r},r,baseNum,rotate);
    if verbose,mlrDispPercent(inf);end
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
if verbose>1,mlrDispPercent(-inf,'Drawing ROI');,end

% get which color to draw the selected ROI in
selectedROIColor = mrGetPref('selectedROIColor');
if isempty(selectedROIColor),selectedROIColor = [1 1 1];end

c=0;
for r = order
  c = c+1;
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
      mlrDispPercent(-inf,sprintf('Computing ROI image coordinates for %i:%s',r,viewGet(view,'roiName',r)));
    end
    [x y s] = getROIImageCoords(view,roi{r}.roiBaseCoords,sliceIndex,baseNum,baseCoordsHomogeneous,imageDims);
    % keep the coordinates
    roi{r}.(baseName){sliceIndex}.x = x;
    roi{r}.(baseName){sliceIndex}.y = y;
    roi{r}.(baseName){sliceIndex}.s = s;
    if verbose, mlrDispPercent(inf); end
    % save in cache
    view = viewSet(view,'ROICache',roi{r},r,baseNum,rotate);
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
    baseSurface = getBaseSurface(view,baseNum); %get baseSurface coordinates, converted to different base space if necessary
    if 0 %%doPerimeter
      if verbose, mlrDispPercent(-inf,'(refreshMLRDisplay) Computing perimeter'); end
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
	if verbose, mlrDispPercent(i/length(y)); end;
      end
      if verbose, mlrDispPercent(-inf); end;
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
    if ~isempty(fig) && ~isempty(hAxis)
      hSurface(c) = patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', roiColors,'facecolor','interp','edgecolor','none','FaceAlpha',0.4,'Parent',hAxis);
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
    view = viewSet(view,'ROICache',roi{r},r,baseNum,rotate);
    % now render those lines
    if ~isempty(hAxis)
      line(roi{r}.lines.x,roi{r}.lines.y,'Color',roi{r}.color,'LineWidth',mrGetPref('roiContourWidth'),'Parent',hAxis);
    end
  else
    roi{r}.lines.x = [];
    roi{r}.lines.y = [];
  end
end
if baseType == 2
  if ~isempty(order) % store surface ROI handles for use in mrInterogator
    view = viewSet(view,'surfaceROIHandle',hSurface,baseNum);
  end
  return
end

% label them, if labelROIs is set
if labelROIs && ~isempty(hAxis)
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


if verbose>1,mlrDispPercent(inf);,end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    displaySurfaceOnVolume    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view = displaySurfaceOnVolume(view,axis,currentSlice,sliceIndex,rotate,sliceDims,slice,baseNum)

colors = [.9 .4 .9;... % color of first inner surface
          .4 .7 .9;... % color of first outer surface
          .9 .4 .4;... % color of second inner surface
          .4 .9 .7;... % color of second outer surface
          .9 .9 .4;... % color of third inner surface
          .6 .4 .9];   % color of third outer surface
% colors = [.4 .4 .4;... % color of first inner surface GREY
%           .2 .2 .2]; % color of first outer surface GREY
numberDifferentColors=1; %maximum number of inner+outer with different colors (up to 3)

cSurf=0;        
for iBase = 1:viewGet(view,'numbase')
  if viewGet(view,'basetype',iBase)==2 && viewGet(view,'baseMultiDisplayContours',iBase)
    if  mod(rotate,90)
      mrWarnDlg('(refreshMLRDisplay) Surface contours not displayed for rotations that are not 0 or multiples of 90');
      return
    elseif isempty(viewGet(view,'base2mag'))
      mrWarnDlg('(refreshMLRDisplay) Cannot display surface contours because the base coordinate system is not defined');
      return
    end
    surfSlice = viewGet(view,'baseCache',baseNum,slice,sliceIndex,rotate,iBase);
    if isempty(surfSlice)

      cSurf=cSurf+1;
      baseCoordMap = viewGet(view,'baseCoordmap',iBase,0);

      innerCoords = permute(baseCoordMap.innerCoords,[2 4 1 3]);
      outerCoords = permute(baseCoordMap.outerCoords,[2 4 1 3]);
      %convert surface coordinates to base coordinates
      base2surf=viewGet(view,'base2base',iBase);
      innerCoords = (base2surf\[innerCoords ones(size(innerCoords,1),1)]')';
      outerCoords = (base2surf\[outerCoords ones(size(outerCoords,1),1)]')';

      % %switch x and y because that's the way everything is plotted in the GUI
      % innerCoords = innerCoords(:,[2 1 3]);
      % outerCoords = outerCoords(:,[2 1 3]);

      %compute the intersection of the surface mesh with the slice
      %each triangle intersects with the current slice plane along a segment defined by two intersections
      %first the WM/GM boundary
%       [surfSlice.firstIntersectionWM,surfSlice.secondIntersectionWM]=computeIntersection(baseCoordMap.tris,innerCoords,currentSlice,sliceIndex,rotate, sliceDims);
%       %and the GM/CSF boundary
%       [surfSlice.firstIntersectionGM,surfSlice.secondIntersectionGM]=computeIntersection(baseCoordMap.tris,outerCoords,currentSlice,sliceIndex,rotate, sliceDims);
      surfSlice.contoursWM=computeIntersection(baseCoordMap.tris,innerCoords,currentSlice,sliceIndex,rotate, sliceDims);
      %and the GM/CSF boundary
      surfSlice.contoursGM=computeIntersection(baseCoordMap.tris,outerCoords,currentSlice,sliceIndex,rotate, sliceDims);

      view = viewSet(view,'baseCache',surfSlice,baseNum,slice,sliceIndex,rotate,iBase);
    end

    %plot intersection contours for WM/GM
    if ~isempty(surfSlice.contoursWM)
      for i=1:length(surfSlice.contoursWM)
        line(surfSlice.contoursWM{i}(:,1), surfSlice.contoursWM{i}(:,2),...
            'color',colors(2*rem(cSurf-1,numberDifferentColors)+1,:),'Parent',axis);
      end
    end
    %and for GM/CSF
    if ~isempty(surfSlice.contoursWM)
      for i=1:length(surfSlice.contoursGM)
        line(surfSlice.contoursGM{i}(:,1), surfSlice.contoursGM{i}(:,2),...
            'color',colors(2*rem(cSurf-1,numberDifferentColors)+2,:),'Parent',axis);
      end
    end
%     %plot all the segments between the first and second intersections for WM/GM 
%     %(this is much slower than plotting connected contours as above)
%     h = line([surfSlice.firstIntersectionWM(:,1) surfSlice.secondIntersectionWM(:,1)]',...
%           [surfSlice.firstIntersectionWM(:,2) surfSlice.secondIntersectionWM(:,2)]',...
%           'color',colors(2*rem(cSurf-1,numberDifferentColors)+1,:),'Parent',axis);
%     h = line([surfSlice.firstIntersectionGM(:,1) surfSlice.secondIntersectionGM(:,1)]',...
%           [surfSlice.firstIntersectionGM(:,2) surfSlice.secondIntersectionGM(:,2)]',...
%           'color',colors(2*rem(cSurf-1,numberDifferentColors)+2,:),'Parent',axis);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     computeIntersection     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function contours=computeIntersection(triangles,vertices,currentSlice,sliceIndex,rotate, sliceDims)
axes2D=fliplr(setdiff([1 2 3],sliceIndex)); %axes of the display depending on the base slice index (slice orientation)
                                             %(this is set in getBaseSlice and axes are switched for unknown reasons)
switch(mod(rotate,360))   % the order of the slice dimensions depends on the set rotation
  case {-180,0,180}       % (not sure why, but didn't look into it too much)
    sliceDims = fliplr(sliceDims);
  case {-270,-90,90,270}
end

%look at x and y coordinates of each triangles vertices and keep those
%that have at least one vertex within sliceDims (irrespective of the z coordinates)
isWithinSlice = reshape(vertices(triangles,axes2D),[size(triangles) 2]);
isWithinSlice = any(isWithinSlice(:,:,1)>0.5 & isWithinSlice(:,:,1)<sliceDims(1)+0.5 & ...
                isWithinSlice(:,:,2)>0.5 & isWithinSlice(:,:,2)<sliceDims(2)+0.5,2);
triangles = triangles(isWithinSlice,:);
% axes2D=fliplr(setdiff([1 2 3],sliceIndex));
%look at the z coordinates of each triangle's vertices and see if it is above (1) or below (-1) the centre of the current slice
isAboveCurrentSlice = reshape(vertices(triangles,sliceIndex),size(triangles))>currentSlice;
% we can use this to find
%1) which triangles to keep (those whose vertices are both above and below the current slice)
trianglesToKeep = find(any(isAboveCurrentSlice,2) & ~all(isAboveCurrentSlice,2));
triangles = triangles(trianglesToKeep,:);
isAboveCurrentSlice = isAboveCurrentSlice(trianglesToKeep,:);
%2) which sides of the triangles intersect with the current slice (those between vertices that are on both sides of the current slice)
%find the vertices that are alone on one side of the slice
loneVertices = [isAboveCurrentSlice(:,1)~=isAboveCurrentSlice(:,2) & isAboveCurrentSlice(:,1)~=isAboveCurrentSlice(:,3) ...
                isAboveCurrentSlice(:,2)~=isAboveCurrentSlice(:,1) & isAboveCurrentSlice(:,2)~=isAboveCurrentSlice(:,3) ...
                isAboveCurrentSlice(:,3)~=isAboveCurrentSlice(:,1) & isAboveCurrentSlice(:,3)~=isAboveCurrentSlice(:,2)]';
triangles = triangles'; %transpose so that we can index more easily
pairVertices = reshape(triangles(~loneVertices),2,size(triangles,2))';
loneVertices = triangles(loneVertices);
%find the intersection between each of these segments and the centre of the current slice
% firstIntersection = currentSlice*ones(size(triangles,2),3);
a = (currentSlice-vertices(loneVertices,sliceIndex))./(vertices(pairVertices(:,1),sliceIndex) - vertices(loneVertices,sliceIndex));
firstIntersection(:,1) = vertices(loneVertices,axes2D(1)) + a .* (vertices(pairVertices(:,1),axes2D(1)) - vertices(loneVertices,axes2D(1)));
firstIntersection(:,2) = vertices(loneVertices,axes2D(2)) + a .* (vertices(pairVertices(:,1),axes2D(2)) - vertices(loneVertices,axes2D(2)));
% secondIntersection = currentSlice*ones(size(triangles,2),3);
b = (currentSlice-vertices(loneVertices,sliceIndex))./(vertices(pairVertices(:,2),sliceIndex) - vertices(loneVertices,sliceIndex));
secondIntersection(:,1) = vertices(loneVertices,axes2D(1)) + b .* (vertices(pairVertices(:,2),axes2D(1)) - vertices(loneVertices,axes2D(1)));
secondIntersection(:,2) = vertices(loneVertices,axes2D(2)) + b .* (vertices(pairVertices(:,2),axes2D(2)) - vertices(loneVertices,axes2D(2)));

%apply rotation to 2D coordinates (only working for rotations that are multiples of 90 deg)
switch mod(rotate,360)
  case 0
    %do nothing
  case {90,-270}
    %flip x coordinates
    firstIntersection(:,1) = sliceDims(1)+1 - firstIntersection(:,1);
    secondIntersection(:,1) = sliceDims(1)+1 - secondIntersection(:,1);
    %swap x and y
    firstIntersection = firstIntersection * [0 1 ;1 0];
    secondIntersection = secondIntersection * [0 1 ;1 0];
  case {180,-180}
    %flip x and y coordinates
    firstIntersection(:,1) = sliceDims(1)+1 - firstIntersection(:,1);
    firstIntersection(:,2) = sliceDims(2)+1 - firstIntersection(:,2);
    secondIntersection(:,1) = sliceDims(1)+1 - secondIntersection(:,1);
    secondIntersection(:,2) = sliceDims(2)+1 - secondIntersection(:,2);
  case {270,-90}
    %flip y coordinates
    firstIntersection(:,2) = sliceDims(2)+1 - firstIntersection(:,2);
    secondIntersection(:,2) = sliceDims(2)+1 - secondIntersection(:,2);
    %swap x and y
    firstIntersection = firstIntersection * [0 1 ;1 0];
    secondIntersection = secondIntersection * [0 1 ;1 0];
  otherwise
    %the more general case would require to figure out how exactly imrotate pads the rotated image
    %(see getBaseSlice.m). No time to do that + not very useful
    %However, this is a start, but has not been tested
    %first need to center the coordinates at the centre of the original (non-rotated) slice 
    %(not available in this function at the moment)
    
    %then rotate the coordinates
    rotateRad = rotate/180*pi;
    xform = eye(3);
    xform(1:2,1:2) = [cos(rotateRad) -sin(rotateRad); sin(rotateRad) cos(rotateRad)];
    firstIntersection = (xform * [firstIntersection ones(size(firstIntersection,1),1)]')';
    secondIntersection = (xform * [secondIntersection ones(size(secondIntersection,1),1)]')';
    
    %then would need to de-center the coordinates using the dimensions of the rotated image output by imrotate
    %this depends on how exactly imrotate pads the rotated image and is probably wrong
    firstIntersection(:,1) = sliceDims(1)/2 + firstIntersection(:,1);
    firstIntersection(:,2) = sliceDims(2)/2 + firstIntersection(:,2);
    secondIntersection(:,1) = sliceDims(1)/2 + secondIntersection(:,1);
    secondIntersection(:,2) = sliceDims(2)/2 + secondIntersection(:,2);
    
    % update 22/07/2016: now that imrotate has been replaced by mrImRotate, should be easy to
    % look inside mrImrotate and take the padding code (probably not worth bothering though)
end 

%now connect contiguous segments to speed-up plotting
c=0;
contours = [];
while ~isempty(firstIntersection)
  c=c+1;
  contours{c} = [firstIntersection(1,:);secondIntersection(1,:)]; % start from first segment coords
  firstIntersection(1,:)=[];  %remove this segment from  both the first and second intersection segment coords
  secondIntersection(1,:)=[];
  % find segment sharing one intersection with first segment
  % this intersection could be found either in the remaing first or second intersections
  [~,nextSegment] =  ismember(contours{c}(end,:),firstIntersection,'rows'); 
  if ~nextSegment
    [~,nextSegment] =  ismember(contours{c}(end,:),secondIntersection,'rows');
    foundInFirst = false;
  else
    foundInFirst = true;
  end
  while nextSegment %repeat while a shared intersection is found
    if foundInFirst %if the shared intersection is in the first intersections,  
      %then the next intersection in the contour is in the second
      contours{c}(end+1,:) = secondIntersection(nextSegment,:); %add it to the list of contour coordinates
    else %and vice versa
      contours{c}(end+1,:) = firstIntersection(nextSegment,:);
    end      
    firstIntersection(nextSegment,:)=[]; %remove current segment from both
    secondIntersection(nextSegment,:)=[];  %the first and second intersection coords
    [~,nextSegment] =  ismember(contours{c}(end,:),firstIntersection,'rows');  % ..etc..
    if ~nextSegment
      [~,nextSegment] =  ismember(contours{c}(end,:),secondIntersection,'rows');
      foundInFirst = false;
    else
      foundInFirst = true;
    end
  end
end

function kernel = gaussianKernel2D(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
kernelDims = 2*[w w]+1;
kernelCenter = ceil(kernelDims/2);
[X,Y] = meshgrid(1:kernelDims(1),1:kernelDims(2));
kernel = exp(-((X-kernelCenter(1)).^2+(Y-kernelCenter(2)).^2)/(2*sigma_d^2)); %Gaussian function
kernel = kernel./sum(kernel(:));


