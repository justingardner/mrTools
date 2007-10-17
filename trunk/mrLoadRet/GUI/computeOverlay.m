% computeOverlay.m
%
%        $Id$
%      usage: overlay = computeOverlay(view,baseXform,baseCoordsHomogeneous,baseDims)
%         by: David Heeger
%       date: 10/16/07
%    purpose: this used to live within the refreshMLRDisplay
%             function, but has been pulled out so that it
%             can be called by other functions
%
function overlay = computeOverlay(view,baseXform,baseCoordsHomogeneous,baseDims)

% check arguments
if ~any(nargin == [4])
  help computeOverlay
  return
end

% view get some stuff
curOverlay = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
alpha = viewGet(view,'alpha');
scan = viewGet(view,'curscan');
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
  interpMethod = 'linear';
end
interpExtrapVal = NaN;

% pull out overlay for this slice
[overlayImages,overlay.coords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(view,scan,baseXform,baseCoordsHomogeneous,baseDims,...
			 analysisNum,interpMethod,interpExtrapVal);
if ~isempty(overlayImages)
  overlayIm = overlayImages(:,:,curOverlay);
else
  overlayIm = [];
end

% get the number of the overlay that should be
% used as an alpha map. This allows the alpha
% to be set according to another overlay
alphaOverlay = viewGet(view,'overlayNum',viewGet(view,'alphaOverlay'));

numOverlays = viewGet(view,'numberofOverlays');
mask = ones(baseDims);
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

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getOverlaySlice   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
  getOverlaySlice(view,scanNum,baseXform,baseCoordsHomogeneous,imageDims,...
		       analysisNum,interpMethod,interpExtrapVal);
%
% getOverlaySlice: extracts overlay image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
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

%%%%%%%%%%%%%%%%%%%%%%%
%%   corAnalInterp   %%
%%%%%%%%%%%%%%%%%%%%%%%
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
  ang = angle(zInterp);
  ang(ang < 0) = ang(ang < 0)+2*pi;
  overlayImages(:,:,phNum) = ang;
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


