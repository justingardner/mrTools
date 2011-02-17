% computeOverlay.m
%
%        $Id$
%      usage: overlays = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,<alpha>)
%         by: David Heeger
%       date: 10/16/07
%    purpose: this used to live within the refreshMLRDisplay
%             function, but has been pulled out so that it
%             can be called by other functions (e.g. mrPrint)-jg
%
function overlays = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,alpha)

% check arguments
if ~any(nargin == [4 5])
  help computeOverlay
  return
end

% viewGet some stuff
curOverlays = viewGet(thisView,'currentOverlay');
nCurOverlays = length(curOverlays);
analysisNum = viewGet(thisView,'currentAnalysis');
scan = viewGet(thisView,'curscan');
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
  interpMethod = 'linear';
end
interpExtrapVal = NaN;

% pull out overlays for this slice
[overlayImages,overlays.coords] = ...
    getOverlaySlice(thisView,scan,base2overlay,baseCoordsHomogeneous,baseDims,...
			 analysisNum,interpMethod,interpExtrapVal);
%put slices on 4th dimensions
overlayImages = permute(overlayImages,[1 2 4 3]);

if ~isempty(overlayImages) && ~isempty(curOverlays)
  %get the mask of clipped values across all overlays
  mask= maskOverlay(thisView,[],[],overlayImages);
  mask = mask{1};
  overlays.overlayIm = overlayImages(:,:,:,curOverlays);
  %alpha based on the alpha slider
  if nargin == 4
    alpha = viewGet(thisView,'alpha');
  end
  overlays.alphaMap = repmat(alpha*mask,[1 1 3]); 
  %alpha maps based on each overlay's alphaoverlay or alpha value
  overlays.alphaMaps=repmat(NaN(size(overlays.overlayIm)),[1 1 3 1]);
  overlays.RGB=repmat(NaN(size(overlays.overlayIm)),[1 1 3 1]);
  for iOverlay = 1:nCurOverlays
    % get the number of the overlays that should be
    % used as an alpha map. This allows the alpha
    % to be set according to another overlays
    alphaOverlayNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',curOverlays(iOverlay)));
    % Finally, make the alphaMap. Normally this is just set
    % to 0 or 1 so that voxels are hard thresholded. If the
    % overlays has an alphaOverlayNum field set to the name
    % of another overlays, then we will use the values in
    % that overlays to set the alpha
    if nargin == 4
      alpha = viewGet(thisView,'alpha',curOverlays(iOverlay));
    end
    if isempty(alphaOverlayNum)
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alpha*mask,[1 1 3]);
    else
      %get the alpha overlay data
      alphaOverlayImage = overlayImages(:,:,:,alphaOverlayNum);
      % get the range of the alpha overlays
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
      if alphaOverlayExponent<0   % if the alpha overlays exponent is negative, set it positive and take the 1-alpha for the alpha map
         alphaOverlayExponent = -alphaOverlayExponent;
         alphaOverlayImage = 1-alphaOverlayImage;
      end
      alphaOverlayImage = alphaOverlayImage.^alphaOverlayExponent;
      alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alphaOverlayImage.*mask,[1 1 3]); 
    end   

    % Rescale current overlays.
    overlays.cmap(:,:,iOverlay) = viewGet(thisView,'overlayCmap',curOverlays(iOverlay));
    switch(viewGet(thisView,'overlayCtype',curOverlays(iOverlay)))
      case 'normal'
        overlays.colorRange(iOverlay,:) = viewGet(thisView,'overlayColorRange',curOverlays(iOverlay));
      case 'setRangeToMax'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(mask))
          overlays.colorRange(1) = max(overlays.colorRange(1),min(overlays.overlayIm(mask)));
          overlays.colorRange(2) = min(max(overlays.overlayIm(mask)),overlays.colorRange(2));
        end
      case 'setRangeToMaxAroundZero'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(mask))
          maxval = max(abs(overlays.overlayIm(mask)));
          overlays.colorRange(1) = -maxval;
          overlays.colorRange(2) = maxval;
        end
      case 'setRangeToMaxAcrossSlices'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum,scan); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum,scan);
        if ~isempty(minOverlay)
          overlays.colorRange(1) = max(minOverlay,overlays.colorRange(1));
        end
        if ~isempty(maxOverlay)
          overlays.colorRange(2) = min(maxOverlay,overlays.colorRange(2));
        end
      case 'setRangeToMaxAcrossSlicesAndScans'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum);
        if ~isempty(minOverlay)
          overlays.colorRange(1) = max(minOverlay,overlays.colorRange(1));
        end
        if ~isempty(maxOverlay)
          overlays.colorRange(2) = min(maxOverlay,overlays.colorRange(2));
        end
    end
    overlays.RGB(:,:,:,iOverlay) = rescale2rgb(overlays.overlayIm(:,:,:,iOverlay),overlays.cmap(:,:,iOverlay),overlays.colorRange(iOverlay,:));
  end
else
  overlays.overlayIm = [];
  overlays.RGB = [];
end

%%%%%%%%%%%%%%%%%%%%%%%
%   getOverlaySlice   %
%%%%%%%%%%%%%%%%%%%%%%%
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
  getOverlaySlice(thisView,scanNum,base2overlay,baseCoordsHomogeneous,imageDims,...
		       analysisNum,interpMethod,interpExtrapVal)
%
% getOverlaySlice: extracts overlays image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
numOverlays = viewGet(thisView,'numberofoverlays',analysisNum);
interpFnctn = viewGet(thisView,'overlayInterpFunction',analysisNum);
tic
% Transform base coords corresponding to this slice/image to overlays
% coordinate frame.
if ~isempty(base2overlay) & ~isempty(baseCoordsHomogeneous) 
  % Transform coordinates
  if size(baseCoordsHomogeneous,3)>1%if it is a flat map with more than one depth
    corticalDepth = viewGet(thisView,'corticalDepth');
    corticalDepthBins = viewGet(thisView,'corticalDepthBins');
    corticalDepths = 0:1/corticalDepthBins:1;
    slices = corticalDepths>=corticalDepth(1) & corticalDepths<=corticalDepth(end);
    baseCoordsHomogeneous = baseCoordsHomogeneous(:,:,slices);
    nDepths = nnz(slices);
  else
    nDepths=1;
  end
  overlayCoordsHomogeneous = base2overlay * reshape(baseCoordsHomogeneous,4,prod(imageDims)*nDepths);
  %temporarily put depth dimensions with y dimension
  overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[imageDims(1) imageDims(2)*nDepths 3 ]);
  overlayCoordsHomogeneous = reshape(overlayCoordsHomogeneous,[4 prod(imageDims) nDepths]);
end


% Extract overlayImages

%if interpolation method is 'nearest', rounding the base coordinates 
%before interpolation should not change the result
%then, using unique to exclude duplicate coordinates,
%interpolation can be performed on far less voxels (at least for surfaces, 
%or if the base has a higher resolution than the scan)
%and this can substantially reduce the time needed for interpolation
if strcmp(interpMethod,'nearest') 
  overlayCoords = round(overlayCoords);   
  [overlayCoords,dump,coordsIndex] = unique(reshape(overlayCoords,[prod(imageDims)*nDepths 3]),'rows'); 
  overlayCoords = permute(overlayCoords,[1 3 2]); 
end
  
if ~isempty(interpFnctn)
  overlayImages = feval(interpFnctn,thisView,scanNum,imageDims,...
    analysisNum,overlayCoords,interpMethod,interpExtrapVal);
else
  overlayImages = zeros([size(overlayCoords,1) size(overlayCoords,2) numOverlays]);
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

if strcmp(interpMethod,'nearest')
  overlayImages = mean(permute(reshape(overlayImages(coordsIndex,:,:),[imageDims nDepths numOverlays]),[1 2 4 3]),4);
  overlayCoords = mean(permute(reshape(overlayCoords(coordsIndex,:,:),[imageDims nDepths 3]),[1 2 4 3]),4);
else
  overlayCoords = mean(permute(reshape(overlayCoords,[imageDims nDepths 3 ]),[1 2 4 3]),4);
  overlayImages = mean(permute(reshape(overlayImages,[imageDims nDepths numOverlays]),[1 2 4 3]),4);
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
overlayImages = zeros([size(overlayCoords,1),size(overlayCoords,2),numOverlays]);

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


