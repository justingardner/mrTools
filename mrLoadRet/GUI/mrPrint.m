% mrPrint.m
%
%        $Id$
%      usage: mrPrint(v,<'params', params>,<'useDefault=1'>,<'justGetParams=1'>)
%         by: justin gardner
%       date: 10/04/07
%    purpose: puts a printable version of the data into the graph win
%
%    To change any default parameter within a script:
%       [~,printParams] = mrPrint(thisView,'justGetParams=1','useDefault=1');
%       % change parameters...
%       mrPrint(thisView,'params',printParams)
%
function [f,params] = mrPrint(v,varargin)

f=gobjects(0);
% check arguments
if nargin < 1
  help mrPrint
  return
end

if viewGet(getMLRView,'baseMultiAxis')>0
  mrWarnDlg('(mrPrint) Not implemented for multi-axes display');
  return;
end

mrGlobals;

% get input arguments
getArgs(varargin,{'params=[]','useDefault=0','justGetParams=0','roiSmooth=0','roiLabels=0'});

% get base type
baseType = viewGet(v,'baseType');
visibleROIs = viewGet(v,'visibleROIs');
if baseType<2
  imageDimensions = viewGet(v,'baseDims');
  switch(baseType)
    case 0
      sliceIndex = viewGet(v,'baseSliceIndex');
      imageDimensions = imageDimensions(setdiff(1:3,sliceIndex));
      imageDimensions = imageDimensions([2 1]);
    case 1
      imageDimensions  = imageDimensions(1:2);
  end
end
overlayList = viewGet(v,'curOverlay');
nOverlays = numel(overlayList);

if ieNotDefined('params')
  % default parameters
  defaultTitle = sprintf('%s: %s',getLastDir(MLR.homeDir),viewGet(v,'description'));
  defaultFontSize = 16;
  defaultMosaicNrows = 0;
  defaultMosaicMargins = [.1 .1];
  defaultColorbarScaleFunction = '@(x)x';
  colobarLocs = {'South','North','East','West','None'};
  if nOverlays>1
    colobarLocs = putOnTopOfList('None',colobarLocs);
  end
  defaultColorbarTickNumber = 4;
  
  % first get parameters that the user wants to display
  paramsInfo = {};
  paramsInfo{end+1} = {'title',defaultTitle,'Title of figure'};
  paramsInfo{end+1} = {'fontSize',defaultFontSize,'Font size of the title and colorbar title'};
  paramsInfo{end+1} = {'backgroundColor',{'white','black'},'type=popupmenu','Background color, either white or black'};
  if baseType == 2 % options for surfaces
    % if this surface has inf in the name, then guess that it is inflated and default to thresholding
    baseName = viewGet(v,'baseName');
    if ~isempty(strfind(lower(baseName),'inf')) thresholdCurvature = 1;else thresholdCurvature = 0;end
    paramsInfo{end+1} = {'thresholdCurvature',thresholdCurvature,'type=checkbox','Thresholds curvature so that the surface is two tones rather than has smooth tones'};

    % compute a good threshold value
    [base.im,base.coords,base.coordsHomogeneous] = getBaseSlice(v); % get the base surface
    baseImg = rescale2rgb(base.im,gray(256),viewGet(v,'baseClip'),viewGet(v,'baseGamma')); % and compute mesh (gray) RGB value like in refreshMLRDisplay
    thresholdValue = mean(baseImg(1,(baseImg(1,:,1)==baseImg(1,:,2))&(baseImg(1,:,3)==baseImg(1,:,2)),1));
    thresholdValue = round(thresholdValue*100)/100;

    paramsInfo{end+1} = {'thresholdValue',thresholdValue,'minmax=[0 1]','incdec=[-0.01 0.01]','contingent=thresholdCurvature','Threshold point - all values below this will turn to the thresholdMin value and all values above this will turn to thresholdMax if thresholdCurvature is turned on.'};
    paramsInfo{end+1} = {'thresholdMin',0.2,'minmax=[0 1]','incdec=[-0.1 0.1]','contingent=thresholdCurvature','The color that all values less than thresholdValue will turn to if thresholdCurvature is set.'};
    paramsInfo{end+1} = {'thresholdMax',0.5,'minmax=[0 1]','incdec=[-0.1 0.1]','contingent=thresholdCurvature','The color that all values greater than thresholdValue will turn to if thresholdCurvature is set.'};
  else % options for flat maps and volume slices
    paramsInfo{end+1} = {'cropX',[1 imageDimensions(2)],sprintf('minmax=[1 %d]',imageDimensions(2)),'incdec=[-10 10]','type=array','X coordinates of a rectangle in pixels to crop the image ([xOrigin width]), before upsampling. X origin is on the left of the image. Not implemented for surfaces'};
    paramsInfo{end+1} = {'cropY',[1 imageDimensions(1)],sprintf('minmax=[1 %d]',imageDimensions(1)),'incdec=[-10 10]','type=array','Y coordinates of a rectangle in pixels to crop the image ([yOrigin height]), before upsampling. Y origin is at the top of the image. Not implemented for surfaces'};
    paramsInfo{end+1} = {'upSampleFactor',0,'type=numeric','round=1','incdec=[-1 1]','minmax=[0 inf]','How many to upsample image by. Each time the image is upsampled it increases in dimension by a factor of 2. So, for example, setting this to 2 will increase the image size by 4'};
    if baseType == 1
      paramsInfo{end+1} = {'maskType',{'Circular','Remove black','None'},'type=popupmenu','Masks out anatomy image for flat maps. Circular finds the largest circular aperture to view the anatomy through. Remove black keeps the patch the same shape, but removes pixels at the edge that are black.'};
    end
  end
  if nOverlays>1
    paramsInfo{end+1} = {'mosaic',false,'type=checkbox','Displays each overlay in a separate panel'};
    paramsInfo{end+1} = {'mosaicNrows',defaultMosaicNrows,'incdec=[-1 1]',sprintf('minmax=[1 %d]',nOverlays),'contingent=mosaic','Number of rows in the mosaic. 0 = set the number of rows automatically'};
    paramsInfo{end+1} = {'mosaicMargins',defaultMosaicMargins,'incdec=[-.01 .01]','minmax=[0 1]','contingent=mosaic','X and Y margins between images, expressed as a proportion of each image''s width and height'};
    contingentString = 'contingent=mosaic';
    colorbarTitle = '';
  else
    contingentString = '';
    colorbarTitle = viewGet(v,'overlayName');
  end
  paramsInfo{end+1} = {'colorbarLoc',colobarLocs,'type=popupmenu',contingentString,'Location of colorbar, select ''None'' if you do not want a colorbar'};
  paramsInfo{end+1} = {'colorbarTitle',colorbarTitle,contingentString,'Title of the colorbar'};
  if nOverlays==1
    paramsInfo{end+1} = {'colorbarScale',viewGet(v,'overlayColorRange'),'type=array',contingentString,'Lower and upper limits of the color scale to display on the color bar'};
  end
  paramsInfo{end+1} = {'colorbarScaleFunction',defaultColorbarScaleFunction,'type=string',contingentString,'Anonymous function to apply to the colorbar scale values [e.g. @(x)exp(x) for data on a logarithmic scale]. This will be applied after applying the colorbarScale parameter. The function must accept and return a one-dimensional array of color scale values.'};
  paramsInfo{end+1} = {'colorbarTickNumber',defaultColorbarTickNumber,'type=numeric','round=1','incdec=[-1 1]','minmax=[2 inf]',contingentString,'Number of ticks on the colorbar'};
  if ~isempty(visibleROIs)
    if baseType == 2 % ROI options for surfaces
      paramsInfo{end+1} = {'roiAlpha',0.4,'minmax=[0 1]','incdec=[-0.1 0.1]','Sets the alpha of the ROIs'};
    else % ROI options for flatmaps and images
      paramsInfo{end+1} = {'roiLineWidth',mrGetPref('roiContourWidth'),'incdec=[-1 1]','minmax=[0 inf]','Line width for drawing ROIs. Set to 0 if you don''t want to display ROIs.'};
      paramsInfo{end+1} = {'roiColor',putOnTopOfList('default',color2RGB),'type=popupmenu','Color to use for drawing ROIs. Select default to use the color currently being displayed.'};
      paramsInfo{end+1} = {'roiOutOfBoundsMethod',{'Remove','Max radius'},'type=popupmenu','If there is an ROI that extends beyond the circular aperture, you can either not draw the lines (Remove) or draw them at the edge of the circular aperture (Max radius). This is only important if you are using a circular aperture.'};
      paramsInfo{end+1} = {'roiLabels',roiLabels,'type=checkbox','Print ROI name at center coordinate of ROI'};
      if baseType == 1
        paramsInfo{end+1} = {'roiSmooth',roiSmooth,'type=checkbox','Smooth the ROI boundaries'};
        paramsInfo{end+1} = {'whichROIisMask',0,'incdec=[-1 1]', sprintf('minmax=%s',mat2str([0 length(visibleROIs)])) 'Which ROI to use as a mask. 0 does no masking'};
        paramsInfo{end+1} = {'filledPerimeter',1,'type=numeric','round=1','minmax=[0 1]','incdec=[-1 1]','Fills the perimeter of the ROI when drawing','contingent=roiSmooth'};
      end
    end
  end

  if useDefault
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo,'Print figure options');
  end
  
  if ~isempty(params)  % if some fields are undefined (because of parameter dependencies), use the defaults
    if fieldIsNotDefined(params,'mosaic')
      params.mosaic = false;
    end
    if fieldIsNotDefined(params,'mosaicNrows')
      params.mosaicNrows = defaultMosaicNrows;
    end
    if fieldIsNotDefined(params,'mosaicMargins')
      params.mosaicMargins = defaultMosaicMargins;
    end
    if fieldIsNotDefined(params,'colorbarLoc')
      params.colorbarLoc = 'None';
    end
    if fieldIsNotDefined(params,'colorbarTitle')
      params.colorbarTitle = colorbarTitle;
    end
    if fieldIsNotDefined(params,'colorbarScaleFunction')
      params.colorbarScaleFunction = defaultColorbarScaleFunction;
    end
    if fieldIsNotDefined(params,'colorbarTickNumber')
      params.colorbarTickNumber = defaultColorbarTickNumber;
    end
  end
end

if isempty(params) || justGetParams,return,end

if baseType<2
  cropX = params.cropX;
  cropY = params.cropY;
end
% grab the image(s)
if params.mosaic
  nImages = nOverlays;
else
  nImages = 1;
end

if nOverlays>1 && ~params.mosaic
  mrWarnDlg('(mrPrint) Printing colorbar for multiple overlays is not implemented');
end

mlrDispPercent(-inf,'(mrPrint) Rerendering image');
for iImage = 1:nImages
  if nOverlays>1 && params.mosaic
    v = viewSet(v,'curOverlay',overlayList(iImage)); % set each overlay in the view one by one
  end
  [img{iImage}, base, roi, ~, altBase{iImage}] = refreshMLRDisplay(viewGet(v,'viewNum'));
  % get the gui, so that we can extract colorbar
  fig = viewGet(v,'figNum');
  if nOverlays>1 && ~params.mosaic
      cmap{iImage} = [];
  elseif ~isempty(fig) % this won't work with the view doesn't have a GUI figure associated with it
    gui = guidata(fig);
    % grab the colorbar data
    H = get(gui.colorbar,'children');
    cmap{iImage} = get(H(end),'CData');
    cmap{iImage}=squeeze(cmap{iImage}(1,:,:));
    yTicks{iImage} = get(gui.colorbar,'YTick');
    xTicks{iImage} = (get(gui.colorbar,'XTick')-0.5)/length(colormap);
    xTickLabels{iImage} = str2num(get(gui.colorbar,'XTicklabel'));
  else
    cmap{iImage} = [];
  end
end
roi = roi(visibleROIs);
if nOverlays>1 && params.mosaic
  v = viewSet(v,'curOverlay',overlayList); % set the overlays back in the view
  refreshMLRDisplay(viewGet(v,'viewNum'));
end
mlrDispPercent(inf);

% just so the code won't break. roiSmooth is only fro baseType = 1
if ~isfield(params,'roiSmooth') params.roiSmooth = 0;end


% display in graph window
f = selectGraphWin;
clf(f);drawnow;
set(f,'Pointer','watch');drawnow;

set(f,'Name','Print figure');
set(f,'NumberTitle','off');
set(f,'color',params.backgroundColor)

% value to consider to be "black" in image
blackValue = 0;

if baseType~=1,params.maskType = 'None';end

% get the mask(s)
if strcmp(params.maskType,'None')
  mask = zeros(size(img{1}));
elseif strcmp(params.maskType,'Remove black')
  mask(:,:,1) = (base.im<=blackValue);
  mask(:,:,2) = (base.im<=blackValue);
  mask(:,:,3) = (base.im<=blackValue);
elseif strcmp(params.maskType,'Circular')
  % find the largest circular aperture the fits the image
  xCenter = (size(base.im,2)/2);
  yCenter = (size(base.im,1)/2);
  x = (1:size(base.im,2))-xCenter;
  y = (1:size(base.im,1))-yCenter;
  [x, y] = meshgrid(x,y);
  % now compute the distance from the center for
  % every point
  d = sqrt(x.^2+y.^2);
  % now go an find the largest distance for which there
  % are no black voxels
  maxd = max(d(:));
  mind = min(d(:));
  for circd = maxd:-1:mind
    if ~any(base.im(d<=circd)<=blackValue)
      break;
    end
  end
  mask(:,:,1) = (d>circd);
  mask(:,:,2) = (d>circd);
  mask(:,:,3) = (d>circd);
end

% get foregroundColor
if strcmp(params.backgroundColor,'white')
  foregroundColor = [0 0 0];
else
  foregroundColor = [1 1 1];
end

if isfield(params,'upSampleFactor')
  % convert upSampleFactor into power of 2
  upSampleFactor = 2^params.upSampleFactor;
  % up sample if called for
  if upSampleFactor > 1
    upSampMask(:,:,1) = upBlur(double(mask(:,:,1)),params.upSampleFactor);
    upSampMask(:,:,2) = upBlur(double(mask(:,:,2)),params.upSampleFactor);
    upSampMask(:,:,3) = upBlur(double(mask(:,:,3)),params.upSampleFactor);
    for iImage = 1:nImages
      upSampImage(:,:,1) = upSample(img{iImage}(:,:,1),params.upSampleFactor);
      upSampImage(:,:,2) = upSample(img{iImage}(:,:,2),params.upSampleFactor);
      upSampImage(:,:,3) = upSample(img{iImage}(:,:,3),params.upSampleFactor);
      img{iImage} = upSampImage;
    end
    upSampMask(upSampMask>0) = upSampMask(upSampMask>0)/max(upSampMask(:));
    mask = upSampMask;
    % make sure we clip to 0 and 1
    mask(mask<0) = 0;mask(mask>1) = 1;
    for iImage = 1:nImages
      img{iImage}(img{iImage}<0) = 0;img{iImage}(img{iImage}>1) = 1;
    end
    % fix the parameters that are used for clipping to a circular aperture
    if exist('circd','var')
      circd = circd*upSampleFactor;
      xCenter = xCenter*upSampleFactor;
      yCenter = yCenter*upSampleFactor;
    end
    cropX(1) = (cropX(1)-1)*upSampleFactor+1;
    cropY(1) = (cropY(1)-1)*upSampleFactor+1;
    cropX(2) = cropX(2)*upSampleFactor;
    cropY(2) = cropY(2)*upSampleFactor;
  end
end

if (baseType == 1) && ~isempty(roi) && params.roiSmooth
  % get the roiImage and mask
  [roiImage, roiMask, dataMask] = getROIPerimeterRGB(v,roi,size(img{1}),params);
  % now set img correctly
  [roiY, roiX] = find(roiMask);
  for i = 1:length(roiX)
    for j = 1:3
      if (roiX(i) <= size(img{1},1)) && (roiY(i) <= size(img{1},2))
        for iImage = 1:nImages
          img{iImage}(roiX(i),roiY(i),j) = roiImage(roiY(i),roiX(i),j);
        end
      end
    end
  end
end

% ROI-based masking
if isfield(params,'whichROIisMask') && params.whichROIisMask
  dataMask = permute(dataMask, [2 1]);
  dataMask = 1-repmat(dataMask, [1 1 3]);

  baseMask = base.RGB;
  baseMask = dataMask .* baseMask;
  for iImage = 1:nImages
    img{iImage} = (1-dataMask) .* img{iImage};
    img{iImage} = img{iImage} + baseMask;
  end
end

% mask out the image
if ~strcmp(params.maskType,'None')
  for iImage = 1:nImages
    if strcmp(params.backgroundColor,'white')
      img{iImage} = (1-mask).*img{iImage} + mask;
    else
      img{iImage} = (1-mask).*img{iImage};
    end
  end
end
for iImage = 1:nImages
  img{iImage}(img{iImage}<0) = 0;img{iImage}(img{iImage}>1) = 1;
end

if nImages>1
  if params.mosaicNrows==0
    figPosition = get(f,'position');
    [nRows,nCols] = getArrayDimensions(nImages,figPosition(4)/figPosition(3));
  else
    nRows = params.mosaicNrows;
    nCols = ceil(nImages/nRows);
  end
  xOuterMargin = .2;
  yOuterMargin = .2;
end

for iImage = 1:nImages
  if nImages>1
    subplotPosition = getSubplotPosition(1+ceil(iImage/nRows),1+iImage-floor((iImage-1)/nRows)*nRows,...
                      [xOuterMargin ones(1,nCols) xOuterMargin],[yOuterMargin ones(1,nRows) yOuterMargin],...
                      params.mosaicMargins(1),params.mosaicMargins(2));
    axisHandle{iImage} = axes('parent',f,'position',subplotPosition);
  else
    axisHandle{iImage} = gca(f);
  end

  % set the colormap
  if ~isempty(cmap{iImage})
    colormap(cmap{iImage});
  end

  % now display the images
  if baseType == 2
    % this is the surface display

    curBase = viewGet(v,'curBase');
    for iBase = 1:viewGet(v,'numBase')
      if viewGet(v,'baseType',iBase)>=2
        if viewGet(v,'baseMultiDisplay',iBase) || isequal(iBase,curBase)
    % get the img (returned by refreshMLRDisplay. This is different
    % for each base when we are displaying more than one. Note
    % that this code hasn't been fully tested with all options yet (jg 2/18/2015)
    if iBase ~= curBase
      thisimg = altBase{iImage}(iBase).img;
    else
      thisimg = img{iImage};
    end
    % taken from refreshMLRDisplay
    baseSurface = viewGet(v,'baseSurface',iBase);
    % threshold curvature if asked for
    if params.thresholdCurvature
      % get all grayscale points (assuming these are the ones that are from the surface)
      grayscalePoints = find((thisimg(1,:,1)==thisimg(1,:,2))&(thisimg(1,:,3)==thisimg(1,:,2)));
      % get points less than 0.5
      lowThresholdPoints = grayscalePoints(thisimg(1,grayscalePoints,1) < params.thresholdValue);
      hiThresholdPoints = grayscalePoints(thisimg(1,grayscalePoints,1) >= params.thresholdValue);
      % set the values to the threshold values
      thisimg(1,lowThresholdPoints,:) = params.thresholdMin;
      thisimg(1,hiThresholdPoints,:) = params.thresholdMax;
    end
    % display the surface
    patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', squeeze(thisimg),'facecolor','interp','edgecolor','none','Parent',axisHandle{iImage});
    hold on
    % make sure x direction is normal to make right/right
    set(axisHandle{iImage},'XDir','reverse');
    set(axisHandle{iImage},'YDir','normal');
    set(axisHandle{iImage},'ZDir','normal');
    % set the camera taret to center of surface
    camtarget(axisHandle{iImage},mean(baseSurface.vtcs))
    % set the size of the field of view in degrees
    % i.e. 90 would be very wide and 1 would be ver
    % narrow. 9 seems to fit the whole brain nicely
    camva(axisHandle{iImage},9);
    setMLRViewAngle(v,axisHandle{iImage});
    % draw the rois
    for roiNum = 1:length(roi)
      patch('vertices', baseSurface.vtcs, 'faces', baseSurface.tris,'FaceVertexCData', roi{roiNum}.overlayImage,'facecolor','interp','edgecolor','none','FaceAlpha',params.roiAlpha,'Parent',axisHandle{iImage});
    end
        end
      end
    end
  else

    %crop image
    cropX(1) = min(cropX(1),size(img{iImage},2));
    cropY(1) = min(cropY(1),size(img{iImage},1));
    cropX(2) = min(cropX(2), size(img{iImage},2) - cropX(1))+1;
    cropY(2) = min(cropY(2), size(img{iImage},1) - cropY(1))+1;
    if cropX(1)>1 || cropY(1)>1 || cropX(2)<size(img{iImage},2) || cropY(2)<size(img{iImage},1)
      img{iImage} = img{iImage}( cropY(1)+(0:cropY(2)-1) , cropX(1)+(0:cropX(2)-1) ,:);
    end

    % display the image (this is for flat maps and images)
    image(img{iImage});
  end

  % set axis
  axis(axisHandle{iImage},'equal');
  axis(axisHandle{iImage},'off');
  axis(axisHandle{iImage},'tight');
  hold(axisHandle{iImage},'on');

  % calculate directions
  params.plotDirections = 0;
  if params.plotDirections
    % calculate gradient on baseCoords
    baseCoords = viewGet(v,'cursliceBaseCoords');
    baseCoords(baseCoords==0) = nan;

    fxx = baseCoords(round(end/2),:,1);
    fxx = fxx(~isnan(fxx));
    fxx = fxx(end)-fxx(1);
    fxy = baseCoords(:,round(end/2),1);
    fxy = fxy(~isnan(fxy));
    fxy = fxy(end)-fxy(1);

    fyx = baseCoords(round(end/2),:,2);
    fyx = fyx(~isnan(fyx));
    fyx = fyx(end)-fyx(1);
    fyy = baseCoords(:,round(end/2),2);
    fyy = fyy(~isnan(fyy));
    fyy = fyy(end)-fyy(1);

    fzx = baseCoords(round(end/2),:,3);
    fzx = fzx(~isnan(fzx));
    fzx = fzx(end)-fzx(1);
    fzy = baseCoords(:,round(end/2),3);
    fzy = fzy(~isnan(fzy));
    fzy = fzy(end)-fzy(1);

  %samplingSize = 4;
  %[fxx fxy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,1));
  %[fyx fyy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,2));
  %[fzx fzy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,3));

  % get mean direction
  %fxx = mean(fxx(~isnan(fxx(:))));fxy = mean(fxy(~isnan(fxy)));
  %fyx = mean(fyx(~isnan(fyx(:))));fyy = mean(fyy(~isnan(fyy)));
  %fzx = mean(fzx(~isnan(fzx(:))));fzy = mean(fzy(~isnan(fzy)));

    startx = 0.15;starty = 0.8;maxlength = 0.075;
    scale = maxlength/max(abs([fxx fxy fyx fyy fzx fzy]));

    annotation('textarrow',startx+[scale*fzx 0],starty+[scale*fzy 0],'String','Left','HeadStyle','none');
    annotation('arrow',startx+[0 scale*fzx],starty+[0 scale*fzy]);
    annotation('textarrow',startx+[-scale*fxx 0],starty+[-scale*fxy 0],'String','Dorsal','HeadStyle','none');
    annotation('arrow',startx+[0 -scale*fxx],starty+[0 -scale*fxy]);
    annotation('textarrow',startx+[scale*fyx 0],starty+[scale*fyy 0],'String','Anterior','HeadStyle','none');
    annotation('arrow',startx+[0 scale*fyx],starty+[0 scale*fyy]);
  end

  % display the colormap
  if ~isempty(cmap{1}) && ~strcmpi(params.colorbarLoc,'None')
    if baseType == 2 && ismember(lower(params.colorbarLoc),{'south','north','east'})
      colorbarLoc = params.colorbarLoc; % put the color bar inside the axes for a tighter figure
    else
      colorbarLoc = [params.colorbarLoc 'Outside'];
    end
    H = colorbar(colorbarLoc);
    set(H,'axisLocation','out'); %make sure ticks and labels are pointing away from the figure (not the case by default when colorbar is inside the axes)

    % set the colorbar ticks, making sure to switch
    % them if we have a vertical as opposed to horizontal bar
    if ~fieldIsNotDefined(params,'colorbarTickNumber')
      xTicks{iImage} = linspace(xTicks{iImage}(1),xTicks{iImage}(end),params.colorbarTickNumber);
      xTickLabels{iImage} = linspace(xTickLabels{iImage}(1),xTickLabels{iImage}(end),params.colorbarTickNumber)';
    end
    if ~fieldIsNotDefined(params,'colorbarScale')
      colorScale = viewGet(v,'overlayColorRange');
      xTickLabels{iImage} =(xTickLabels{iImage}-colorScale(1))/diff(colorScale)*diff(params.colorbarScale)+params.colorbarScale(1);
    end
    if ~fieldIsNotDefined(params,'colorbarScaleFunction')
      colorbarScaleFunction = str2func(params.colorbarScaleFunction);
      xTickLabels{iImage} = colorbarScaleFunction(xTickLabels{iImage});
    end
    xTickLabels{iImage} = num2str( xTickLabels{iImage}, '%.3f');
    % remove trailing zeros (can't use %g or %f to directly get both a fixed number of decimal points and no trailing zeros)
    totalColNum = size(xTickLabels{iImage},2);
    for iTick = 1:size(xTickLabels{iImage},1)
      colNum = totalColNum;
      while xTickLabels{iImage}(iTick,colNum)=='0'
        colNum=colNum-1;
      end
      if xTickLabels{iImage}(iTick,colNum)=='.'
        colNum=colNum-1;
      end
      xTickLabels{iImage}(iTick,:) = circshift(xTickLabels{iImage}(iTick,:),-colNum);
      xTickLabels{iImage}(iTick,1:totalColNum-colNum) = repmat(' ',1,totalColNum-colNum);
    end
    if ismember(lower(params.colorbarLoc),{'east','west'})
      set(H,'XTick',yTicks{iImage});
      set(H,'Ytick',xTicks{iImage});
      set(H,'YTickLabel',xTickLabels{iImage});
    else
      set(H,'YTick',yTicks{iImage});
      set(H,'Xtick',xTicks{iImage});
      set(H,'XTickLabel',xTickLabels{iImage});
    end
    set(H,'XColor',foregroundColor);
    set(H,'YColor',foregroundColor);
    if nImages>1 && params.mosaic
      set(H,'tickDirection','out'); % for mosaic display, colorbars are smaller, so orient the ticks outwards
    end

    %color bar title (label)
    set(get(H,'Label'),'String',params.colorbarTitle);
    set(get(H,'Label'),'Interpreter','none');
    set(get(H,'Label'),'Color',foregroundColor);
    set(get(H,'Label'),'FontSize',params.fontSize-4);
    switch(lower(params.colorbarLoc)) % change default position of label depending on location of colorbar
      case 'south'
        set(get(H,'Label'),'position',[0.5 2])
      case 'north'
  %       do nothing
      case 'east'
        set(get(H,'Label'),'position',[-1 0.5]);
        imagePosition = get(axisHandle{iImage},'position');
        imagePosition(1) = imagePosition(1)-0.01;
        set(axisHandle{iImage},'position',imagePosition);
        labelPosition = get(H,'position');
        labelPosition(1) = labelPosition(1)+0.01;
        set(H,'position',labelPosition);
      case 'west'
        set(get(H,'Label'),'position',[1 0.5])
        set(get(H,'Label'),'rotation',270)
    end
  end

  % create a title
  if ~isempty(params.title) && iImage==1
    H = title(params.title);
    set(H,'Interpreter','none');
    set(H,'Color',foregroundColor);
    set(H,'FontSize',params.fontSize);
  end
  
  % draw the roi
  if baseType ~= 2
    if iImage==1
      sliceNum = viewGet(v,'currentSlice');
      label = {};
      mlrDispPercent(-inf,'(mrPrint) Rendering ROIs');
      visibleROIs = viewGet(v,'visibleROIs');
      for rnum = 1:length(roi)
        % check for lines
        if params.roiLineWidth > 0 && ~isempty(roi{rnum}) && isfield(roi{rnum},'lines') && ~isempty(roi{rnum}.lines.x)
          % get color
          if strcmp(params.roiColor,'default')
            color = roi{rnum}.color;
          else
            color = color2RGB(params.roiColor);
          end
          % deal with upSample factor
          roi{rnum}.lines.x = roi{rnum}.lines.x*upSampleFactor;
          roi{rnum}.lines.y = roi{rnum}.lines.y*upSampleFactor;
          % labels for rois, just create here
          % and draw later so they are always on top
          if params.roiLabels
            x = roi{rnum}.lines.x;
            y = roi{rnum}.lines.y;
            label{end+1}.x = median(x(~isnan(x)));
            label{end}.y = median(y(~isnan(y)));
            label{end}.str = viewGet(v,'roiName',visibleROIs(rnum));
            label{end}.color = color;
          end
          % if we have a circular apertuer then we need to
          % fix all the x and y points so they don't go off the end
          if strcmp(params.maskType,'Circular')
            % get the distance from center
            x = roi{rnum}.lines.x-xCenter;
            y = roi{rnum}.lines.y-yCenter;
            d = sqrt(x.^2+y.^2);
            if strcmp(params.roiOutOfBoundsMethod,'Max radius')
        % find the angle of all points
        ang = atan(y./x);
        ysign = (y > 0)*2-1;
        xsign = (x > 0)*2-1;
        newx = circd*cos(ang);
        newy = circd*sin(ang);
        % now reset all points past the maximum radius
        % with values at the outermost edge of the aperture
        x(d>circd) = newx(d>circd);
        y(d>circd) = newy(d>circd);
        % set them back in the structure
        roi{rnum}.lines.x = xsign.*abs(x)+xCenter;
        roi{rnum}.lines.y = ysign.*abs(y)+yCenter;
            else
        % set all values greater than the radius to nan
        x(d>circd) = nan;
        y(d>circd) = nan;

        % set them back in the strucutre
        roi{rnum}.lines.x = x+xCenter;
        roi{rnum}.lines.y = y+yCenter;
            end
          end

          % if cropping, correct x and y coordinates and remove lines falling outside the crop box
          roi{rnum}.lines.x = roi{rnum}.lines.x - cropX(1) + 1;
          roi{rnum}.lines.y(:,all(roi{rnum}.lines.x > cropX(2)+0.5)) = [];
          roi{rnum}.lines.x(:,all(roi{rnum}.lines.x > cropX(2)+0.5)) = [];
          roi{rnum}.lines.x(roi{rnum}.lines.x > cropX(2)+0.5) = cropX(2)+0.5;
          roi{rnum}.lines.y(:,all(roi{rnum}.lines.x < 0.5)) = [];
          roi{rnum}.lines.x(:,all(roi{rnum}.lines.x < 0.5)) = [];
          roi{rnum}.lines.x(roi{rnum}.lines.x < 0.5) = 0.5;
          roi{rnum}.lines.y = roi{rnum}.lines.y - cropY(1) + 1;
          roi{rnum}.lines.x(:,all(roi{rnum}.lines.y > cropY(2)+0.5)) = [];
          roi{rnum}.lines.y(:,all(roi{rnum}.lines.y > cropY(2)+0.5)) = [];
          roi{rnum}.lines.y(roi{rnum}.lines.y > cropY(2)+0.5) = cropY(2)+0.5;
          roi{rnum}.lines.x(:,all(roi{rnum}.lines.y < 0.5)) = [];
          roi{rnum}.lines.y(:,all(roi{rnum}.lines.y < 0.5)) = [];
          roi{rnum}.lines.y(roi{rnum}.lines.y < 0.5) = 0.5;

        end
      end
    end
    mlrDispPercent(inf);

    % draw the lines
    for rnum = 1:length(roi)
      if params.roiLineWidth > 0 && ~isempty(roi{rnum}) && isfield(roi{rnum},'lines') && ~isempty(roi{rnum}.lines.x) && ~params.roiSmooth
        line(roi{rnum}.lines.x,roi{rnum}.lines.y,'Color',color,'LineWidth',params.roiLineWidth);
      end
    end
    
    % draw the labels
    for i = 1:length(label)
      h = text(label{i}.x,label{i}.y,label{i}.str);
      set(h,'Color',foregroundColor);
      set(h,'Interpreter','None');
      set(h,'EdgeColor',label{i}.color);
      set(h,'BackgroundColor',params.backgroundColor);
      set(h,'FontSize',10);
      set(h,'HorizontalAlignment','center');
    end
    
  end
end
set(f,'Pointer','arrow');

% bring up print dialog
global mrPrintWarning
if isempty(mrPrintWarning)
  mrWarnDlg('(mrPrint) If you are having trouble getting colors to print correctly, try exporting the figure to an eps (use its File/Export Setup menu) and printing that');
  mrPrintWarning = 1;
end
%printdlg(f);
%printpreview(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function modified from makeROIPerimeterRGB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [roiRGB,roiMask,dataMask] = getROIPerimeterRGB(v, roi, imageSize, params)
%
% [roiRGB,roiMask,dataMask] = getROIPerimeterRGB(view, imageSize, params)
% 
% Draws a line around the ROIs.  
%
% roiRGB is the RGB color image with the ROI outlines drawn in. 
% roiMask is a mask image (binary) showing where the all the ROIs are. You
% can use this to draw in the roiRGB data without disturbing the rest of
% the image.
% dataMask is created if an ROI that is called 'mask' exists- it is an
% empty matrix if it doesn't exist. It is a binary mask image with ones 
% inside the ROI called 'mask' and zeros outside it. It can be used to 
% show only a selected portion of the data.
%
% 2002.12.18 RFD: added dataMask output and some comments.

if ~exist('lineWidth', 'var')
  lineWidth = 0.5;
end

% make sure the image size is 2D
upsampleFactor = 2^params.upsampleFactor;
upSampImSize = imageSize(1:2)*upSampleFactor;

% Initialize the output RGB image
roiRGB = zeros([upSampImSize 3]);
roiMask = zeros(upSampImSize);
dataMask = [];

% get baseName, used for retrieving image coordinates
% of ROI calculated by refreshMLRDisplay
baseName = fixBadChars(viewGet(v,'baseName'));
sliceIndex = viewGet(v,'baseSliceIndex');
if viewGet(v,'baseType')
  % note we only keep cortical depth to a precision of 1 decimal place
  corticalDepth = sprintf('%0.1g',viewGet(v,'corticalDepth'));
  baseName = sprintf('%s%s',baseName,fixBadChars(corticalDepth,{'.','_'}));
end

mlrDispPercent(-inf,'(mrPrint) Calculating smoothed ROIs');
for r=1:length(roi)
  if ~isempty(roi{r})
    % get the x and y image coordinates
    x = roi{r}.(baseName){sliceIndex}.x;
    y = roi{r}.(baseName){sliceIndex}.y;
    s = roi{r}.(baseName){sliceIndex}.s;
    
    % check for regular images
    baseType = viewGet(v,'baseType');
    if baseType == 0
      %FIX FIX, must select correct voxels for 2D inplanes
      %    voxelsOnThisSlice = find(s==viewGet(v,'curSlice'));
      %    x = x(voxelsOnThisSlice);
      %    y = y(voxelsOnThisSlice);
    end
    
    % upSample the coords
    x = x.*upSampleFactor;
    y = y.*upSampleFactor;
    
    upSampSq = upSampleFactor^2;
    n = length(x)*upSampSq;
    hiResX = zeros(1,n);
    hiResY = zeros(1,n);
    
    for ii=1:upSampleFactor
      offsetX = (ii-upSampleFactor-.5)+upSampleFactor/2;
      for jj=1:upSampleFactor
        offsetY = (jj-upSampleFactor-.5)+upSampleFactor/2;
        hiResX((ii-1)*upSampleFactor+jj:upSampSq:end) = x+offsetX;
        hiResY((ii-1)*upSampleFactor+jj:upSampSq:end) = y+offsetY;
      end
    end
    
    hiResX = round(hiResX);
    hiResY = round(hiResY);
    
    goodVals=find((hiResX>0) & (hiResY>0) & (hiResX<upSampImSize(2)) & (hiResY<upSampImSize(1)));
    
    hiResX=hiResX(goodVals);
    hiResY=hiResY(goodVals);
    
    % Draw the whole ROI into an image
    roiBits = zeros(upSampImSize);
    roiBits(sub2ind(upSampImSize,hiResY,hiResX)) = 1;

    % blur it some, but only need to do this
    % if we haven't already upsampled
    if upSampleFactor <= 2
      roiBits = blur(roiBits,2);
      roiBits(roiBits<median(roiBits(:)))=0;
      roiBits(roiBits~=0) = 1;
    end
    % is it a mask roi?
    if r==params.whichROIisMask
      dataMask = roiBits;
    else
      if (params.filledPerimeter)
        % Do the ROI plotting using the image processing
        % toolbox's morphological ops
        
        % Dilate, erode, fill
        se=strel('disk',32);
        tmat=imdilate(logical(roiBits),se);
        se=strel('disk',32);
        tmat=imerode(logical(tmat),se);
        tmat=imfill(tmat,'holes');
        tmat=tmat-min(tmat(:));
        roiBits=bwperim(tmat);
        se=strel('square',params.roiLineWidth*2);
        
        roiBits=double(imdilate(logical(roiBits),se));
        
      else
        % grow the region and subtract
        % create a circular filter of ones
        e = exp(-([-params.roiLineWidth:params.roiLineWidth]./params.roiLineWidth).^2);
        filt = (e'*e)>.367;
        roiBits = (conv2(roiBits,double(filt),'same')>0.1)-roiBits;
      end
      
      % get color
      if strcmp(params.roiColor,'default')
        color = roi{r}.color;
      else
        color = color2RGB(params.roiColor);
      end
      
      % convert to rgb
      roiRGB(:,:,1) = roiRGB(:,:,1) + roiBits.*color(1);
      roiRGB(:,:,2) = roiRGB(:,:,2) + roiBits.*color(2);
      roiRGB(:,:,3) = roiRGB(:,:,3) + roiBits.*color(3);
      roiMask = roiMask | roiBits;
      mlrDispPercent(r/length(roi));
    end 
  end
end    
mlrDispPercent(inf);

