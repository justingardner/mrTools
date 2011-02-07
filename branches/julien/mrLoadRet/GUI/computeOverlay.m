% computeOverlay.m
%
%        $Id$
%      usage: overlay = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,<alpha>)
%         by: David Heeger
%       date: 10/16/07
%    purpose: this used to live within the refreshMLRDisplay
%             function, but has been pulled out so that it
%             can be called by other functions (e.g. mrPrint)-jg
%
function overlay = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,alpha)

% check arguments
if ~any(nargin == [4 5])
  help computeOverlay
  return
end

% viewGet some stuff
curOverlay = viewGet(thisView,'currentOverlay');
analysisNum = viewGet(thisView,'currentAnalysis');
if nargin == 4
  alpha = viewGet(thisView,'alpha');
end
scan = viewGet(thisView,'curscan');
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
  interpMethod = 'linear';
end
interpExtrapVal = NaN;

% pull out overlay for this slice
[overlayImages,overlay.coords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(thisView,scan,base2overlay,baseCoordsHomogeneous,baseDims,...
			 analysisNum,interpMethod,interpExtrapVal);
%put slices on 4th dimensions
overlayImages = permute(overlayImages,[1 2 4 3]);

% get the number of the overlay that should be
% used as an alpha map. This allows the alpha
% to be set according to another overlay
alphaOverlayNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay'));

if ~isempty(overlayImages) && ~isempty(curOverlay)
  [mask,nonMasked]= maskOverlay(thisView,[curOverlay alphaOverlayNum],[],overlayImages);
  if ~isempty(alphaOverlayNum)
    alphaOverlayImage = nonMasked{2};
  end
  overlay.overlayIm = nonMasked{1};
  mask = mask{1};
else
  overlay.overlayIm = [];
  alphaOverlayImage = [];
end

% %JB: should replace the following by a call to maskOverlay.m, but make sure it does the exact same thing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUT HERE
% numOverlays = viewGet(thisView,'numberofOverlays');
% mask = ones(baseDims);
% % Loop through overlays, filling in NaNs according to clip values.
% if ~isempty(overlayImages)
%   for ov = 1:numOverlays
%     im = overlayImages(:,:,ov);
%     clip = viewGet(thisView,'overlayClip',ov);
%     % Find pixels that are within clip
%     if diff(clip) > 0
%       % show data between the clip values
%       % note: AND
%       pts = (im >= clip(1) & im <= clip(2));
%     elseif diff(clip) < 0
%       % show data more extreme than the clip values
%       % note: OR
%       pts = (im >= clip(1) | im <= clip(2));
%     else
%       % nothing in interval; don't show anything
%       pts = false(size(im));
%     end
%     % do not clip out for any points that are set to nan
%     % this can happen if the current overlay does not
%     % exist for this scan
%     if ov ~= curOverlay
%       pts = pts | isnan(im);
%     end
%     % now make the mask
%     mask = mask & pts;
%     % if this is the alpha overlay then keep it.
%     if ov == alphaOverlayNum
%       alphaOverlayImage = im;
%     end
%   end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUT HERE




if ~isempty(overlay.overlayIm)
  % Finally, make the alphaMap. Normally this is just set
  % to 0 or 1 so that voxels are hard thresholded. If the
  % overlay has an alphaOverlayNum field set to the name
  % of another overlay, then we will use the values in
  % that overlay to set the alpha
  if isempty(alphaOverlayNum)
    overlay.alphaMap = repmat(alpha*mask,[1 1 3]);
  else
    % get the range of the alpha overlay
    range = viewGet(thisView,'overlayRange',alphaOverlayNum);
    % handle setRangeToMax
    if strcmp(viewGet(thisView,'overlayCtype',alphaOverlayNum),'setRangeToMax')
      maxRange = max(clip(1),min(alphaOverlayImage(mask)));
      if ~isempty(maxRange),range(1) = maxRange;end
      minRange = min(max(alphaOverlayImage(mask)),clip(2));
      if ~isempty(minRange),range(2) = minRange;end
    end
    % now compute the alphaOverlayNum as a number from
    % 0 to 1 of the range
    alphaOverlayImage = alpha*((alphaOverlayImage-range(1))./diff(range));
    alphaOverlayImage(alphaOverlayImage>alpha) = alpha;
    alphaOverlayImage(alphaOverlayImage<0) = 0;
    alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
    if alphaOverlayExponent<0   % if the alpha overlay exponent is negative, set it positive and take the 1-alpha for the alpha map
       alphaOverlayExponent = -alphaOverlayExponent;
       alphaOverlayImage = 1-alphaOverlayImage;
    end
    alphaOverlayImage = alphaOverlayImage.^alphaOverlayExponent;
    alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
    overlay.alphaMap = repmat(alphaOverlayImage.*mask,[1 1 3]); 
  end   

  % Rescale current overlay.
  overlay.cmap = viewGet(thisView,'overlayCmap',curOverlay);
  switch(viewGet(thisView,'overlayCtype',curOverlay))
    case 'normal'
      overlay.colorRange = viewGet(thisView,'overlayColorRange',curOverlay);
    case 'setRangeToMax'
      overlay.colorRange = viewGet(thisView,'overlayClip',curOverlay);
      if ~isempty(overlay.overlayIm(mask))
        overlay.colorRange(1) = max(overlay.colorRange(1),min(overlay.overlayIm(mask)));
        overlay.colorRange(2) = min(max(overlay.overlayIm(mask)),overlay.colorRange(2));
      end
    case 'setRangeToMaxAroundZero'
      overlay.colorRange = viewGet(thisView,'overlayClip',curOverlay);
      if ~isempty(overlay.overlayIm(mask))
        maxval = max(abs(overlay.overlayIm(mask)));
        overlay.colorRange(1) = -maxval;
        overlay.colorRange(2) = maxval;
      end
    case 'setRangeToMaxAcrossSlices'
      overlay.colorRange = viewGet(thisView,'overlayClip',curOverlay);
      minOverlay = viewGet(thisView,'minOverlaydata',curOverlay,analysisNum,scan); 
      maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlay,analysisNum,scan);
      if ~isempty(minOverlay)
        overlay.colorRange(1) = max(minOverlay,overlay.colorRange(1));
      end
      if ~isempty(maxOverlay)
        overlay.colorRange(2) = min(maxOverlay,overlay.colorRange(2));
      end
    case 'setRangeToMaxAcrossSlicesAndScans'
      overlay.colorRange = viewGet(thisView,'overlayClip',curOverlay);
      minOverlay = viewGet(thisView,'minOverlaydata',curOverlay,analysisNum); 
      maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlay,analysisNum);
      if ~isempty(minOverlay)
        overlay.colorRange(1) = max(minOverlay,overlay.colorRange(1));
      end
      if ~isempty(maxOverlay)
        overlay.colorRange(2) = min(maxOverlay,overlay.colorRange(2));
      end
  end
  overlay.RGB = rescale2rgb(overlay.overlayIm,overlay.cmap,overlay.colorRange);
else
  overlay.RGB = [];
end

%%%%%%%%%%%%%%%%%%%%%%%
%   getOverlaySlice   %
%%%%%%%%%%%%%%%%%%%%%%%
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
  getOverlaySlice(thisView,scanNum,base2overlay,baseCoordsHomogeneous,imageDims,...
		       analysisNum,interpMethod,interpExtrapVal);
%
% getOverlaySlice: extracts overlay image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
numOverlays = viewGet(thisView,'numberofoverlays',analysisNum);
interpFnctn = viewGet(thisView,'overlayInterpFunction',analysisNum);

% Transform base coords corresponding to this slice/image to overlay
% coordinate frame.
if ~isempty(base2overlay) & ~isempty(baseCoordsHomogeneous)
  % Transform coordinates
  overlayCoordsHomogeneous = base2overlay * baseCoordsHomogeneous;
  overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[imageDims 3]);
end

% Extract overlayImages
if ~isempty(interpFnctn)
  overlayImages = feval(interpFnctn,thisView,scanNum,imageDims,...
    analysisNum,overlayCoords,interpMethod,interpExtrapVal);
else
  overlayImages = zeros([imageDims,numOverlays]);
  for ov = 1:numOverlays
    overlayData = viewGet(thisView,'overlayData',scanNum,ov,analysisNum);
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
function overlayImages = corAnalInterp(thisView,scanNum,...
  imageDims,analysisNum,overlayCoords,interpMethod,interpExtrapVal);
%
% corAnalInterp: special case for corAnal. Need to treat amp/ph as complex
% valued.

% Initialize
numOverlays = viewGet(thisView,'numberofoverlays',analysisNum);
overlayImages = zeros([imageDims,numOverlays]);

% Interpolate complex values for amp and ph
ampNum = viewGet(thisView,'overlayNum','amp',analysisNum);
phNum = viewGet(thisView,'overlayNum','ph',analysisNum);
ampData = viewGet(thisView,'overlaydata',scanNum,ampNum,analysisNum);
phData = viewGet(thisView,'overlaydata',scanNum,phNum,analysisNum);
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
  overlayData = viewGet(thisView,'overlayData',scanNum,ov,analysisNum);
  if ~isempty(overlayData) & ~isempty(overlayCoords)
    % Extract the slice
    overlayImages(:,:,ov) = interp3(overlayData,...
      overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
      interpMethod,interpExtrapVal);
  end
end


