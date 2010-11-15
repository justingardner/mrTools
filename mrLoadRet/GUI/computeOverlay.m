% computeOverlay.m
%
%        $Id$
%      usage: overlay = computeOverlay(view,base2overlay,baseCoordsHomogeneous,baseDims,<alpha>)
%         by: David Heeger
%       date: 10/16/07
%    purpose: this used to live within the refreshMLRDisplay
%             function, but has been pulled out so that it
%             can be called by other functions (e.g. mrPrint)-jg
%
function overlay = computeOverlay(view,base2overlay,baseCoordsHomogeneous,baseDims,alpha)

% check arguments
if ~any(nargin == [4 5])
  help computeOverlay
  return
end

% view get some stuff
curOverlay = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
if nargin == 4
  alpha = viewGet(view,'alpha');
end
scan = viewGet(view,'curscan');
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
  interpMethod = 'linear';
end
interpExtrapVal = NaN;

% pull out overlay for this slice
[overlayImages,overlay.coords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(view,scan,base2overlay,baseCoordsHomogeneous,baseDims,...
			 analysisNum,interpMethod,interpExtrapVal);
if ~isempty(overlayImages)
  overlay.overlayIm = overlayImages(:,:,curOverlay);
else
  overlay.overlayIm = [];
end

% get the number of the overlay that should be
% used as an alpha map. This allows the alpha
% to be set according to another overlay
alphaOverlay = viewGet(view,'overlayNum',viewGet(view,'alphaOverlay'));


%JB: should replace the following by a call to maskOverlay.m, but make sure it does the exact same thing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUT HERE
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
    elseif diff(clip) < 0
      pts = (im >= clip(1) | im <= clip(2));
    else
      pts = false(size(im));
    end
    % do not clip out for any points that are set to nan
    % this can happen if the current overlay does not
    % exist for this scan
    if ov ~= curOverlay
      pts = pts | isnan(im);
    end
    % now make the mask
    mask = mask & pts;
    % if this is the alpha overlay then keep it.
    if ov == alphaOverlay
      alphaOverlayImage = im;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUT HERE


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
  alphaOverlayExponent = viewGet(view,'alphaOverlayExponent');
  if alphaOverlayExponent<0   % if the alpha overlay exponent is negative, set it positive and take the 1-alpha for the alpha map
     alphaOverlayExponent = -alphaOverlayExponent;
     alphaOverlayImage = 1-alphaOverlayImage;
  end
  alphaOverlayImage = alphaOverlayImage.^alphaOverlayExponent;
  alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
  overlay.alphaMap = repmat(alphaOverlayImage.*mask,[1 1 3]); 
end   

% Rescale current overlay.
if ~isempty(overlay.overlayIm)
  overlay.cmap = viewGet(view,'overlayCmap',curOverlay);
  overlay.range = viewGet(view,'overlayRange',curOverlay);
  if strcmp(viewGet(view,'overlayCtype',curOverlay),'setRangeToMax')
    clip = viewGet(view,'overlayClip',curOverlay);
    if ~isempty(overlay.overlayIm(mask))
      overlay.range(1) = max(clip(1),min(overlay.overlayIm(mask)));
      overlay.range(2) = min(max(overlay.overlayIm(mask)),clip(2));
    else
      overlay.range = clip;
    end
  elseif strcmp(viewGet(view,'overlayCtype',curOverlay),'setRangeToMaxAroundZero')
    clip = viewGet(view,'overlayClip',curOverlay);
    if ~isempty(overlay.overlayIm(mask))
      maxval = max(abs(overlay.overlayIm(mask)));
      overlay.range(1) = -maxval;
      overlay.range(2) = maxval;
    else
      overlay.range = clip;
    end
  end
  overlay.RGB = rescale2rgb(overlay.overlayIm,overlay.cmap,overlay.range);
else
  overlay.RGB = [];
end

%%%%%%%%%%%%%%%%%%%%%%%
%   getOverlaySlice   %
%%%%%%%%%%%%%%%%%%%%%%%
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
  getOverlaySlice(view,scanNum,base2overlay,baseCoordsHomogeneous,imageDims,...
		       analysisNum,interpMethod,interpExtrapVal);
%
% getOverlaySlice: extracts overlay image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
numOverlays = viewGet(view,'numberofoverlays',analysisNum);
interpFnctn = viewGet(view,'overlayInterpFunction',analysisNum);

% Transform base coords corresponding to this slice/image to overlay
% coordinate frame.
if ~isempty(base2overlay) & ~isempty(baseCoordsHomogeneous)
  % Transform coordinates
  overlayCoordsHomogeneous = base2overlay * baseCoordsHomogeneous;
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
    else
      overlayImages(:,:,ov) = nan;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    corAnalInterp   %%
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


