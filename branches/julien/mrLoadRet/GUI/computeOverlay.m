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

sliceInfo.base2overlay = base2overlay;
sliceInfo.baseCoordsHomogeneous = baseCoordsHomogeneous;
sliceInfo.baseDims = [baseDims 1];

if ~isempty(curOverlays)
  curAlphaOverlays = zeros(1,nCurOverlays);
  for iOverlay = 1:nCurOverlays
    thisNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',curOverlays(iOverlay)));
    if ~isempty(thisNum)
      curAlphaOverlays(iOverlay) = thisNum;
    end
  end
  curOverlays = [curOverlays curAlphaOverlays];

  %get overlays and masks of clipped values 
  [overlayMasks,overlayImages,overlays.coords]= maskOverlay(thisView,curOverlays,scan,sliceInfo);
else
  overlayImages = [];
end
if ~isempty(overlayImages)
  overlays.coords = overlays.coords{1};
  alphaOverlayImages = overlayImages{1}(:,:,:,nCurOverlays+1:end);
  alphaOverlayMasks = overlayMasks{1}(:,:,:,nCurOverlays+1:end);
  overlayImages = overlayImages{1}(:,:,:,1:nCurOverlays);
  overlayMasks = overlayMasks{1}(:,:,:,1:nCurOverlays);
  
  if size(overlayImages,3)>1
    overlayMasks = any(overlayMasks,3);
    alphaOverlayMasks = any(alphaOverlayMasks,3);
    switch(mrGetPref('multiSliceProjectionMethod'))
      case 'Average'
        overlayImages = nanmean(overlayImages,3);
        alphaOverlayImages = nanmean(alphaOverlayImages,3);
      case 'Maximum Intensity Projection'
        alphaOverlayImages= permute(alphaOverlayImages,[4 3 1 2]);
        newAlphaOverlayImages = NaN([size(overlayImages,4) size(overlayImages,1) size(overlayImages,2)]);
        newOverlayImages = NaN([size(overlayImages,1) size(overlayImages,2) 1 size(overlayImages,4)]);
        for iOverlay = 1:nCurOverlays
          [newOverlayImages(:,:,1,iOverlay),maxIndex] = max(overlayImages(:,:,:,iOverlay),[],3);
          %get the corresponding alpha values (?)
          for iSlice = 1:size(overlayImages,3)
            thisIndex = find(maxIndex==iSlice);
            newAlphaOverlayImages(iOverlay,thisIndex) = alphaOverlayImages(iOverlay,iSlice,thisIndex);
          end
        end
        alphaOverlayImages = permute(newAlphaOverlayImages,[2 3 4 1]);
        overlayImages = newOverlayImages;
    end
  end
  
  overlays.overlayIm = permute(overlayImages,[1 2 4 3]); 
  imageDims = [size(overlayImages,1) size(overlayImages,2)];
  
  
  %alpha maps based on each overlay's alphaoverlay or alpha value
  overlays.alphaMaps=repmat(NaN([imageDims 1 nCurOverlays],mrGetPref('defaultPrecision')),[1 1 3 1]);
  overlays.RGB=repmat(NaN([imageDims 1 nCurOverlays],mrGetPref('defaultPrecision')),[1 1 3 1]);
  
  % Finally, make the alphaMap. Normally this is just set
  % to 0 or 1 so that voxels are hard thresholded. If the
  % overlays has an alphaOverlay field set to the name
  % of another overlays, then we will use the values in
  % that overlays to set the alpha
  for iOverlay = 1:nCurOverlays
    if nargin == 4
      alpha = viewGet(thisView,'alpha',curOverlays(iOverlay));
    end
    if ~curAlphaOverlays(iOverlay)
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alpha*overlayMasks(:,:,:,iOverlay),[1 1 3]);
    else
      %get the alpha overlay data
      alphaOverlayImage = alphaOverlayImages(:,:,:,iOverlay);
      % get the range of the alpha overlays
      range = viewGet(thisView,'overlayRange',curAlphaOverlays(iOverlay));
      % handle setRangeToMax
      if strcmp(viewGet(thisView,'overlayCtype',curAlphaOverlays(iOverlay)),'setRangeToMax')
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
      alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
      if alphaOverlayExponent<0   % if the alpha overlays exponent is negative, set it positive and take the 1-alpha for the alpha map
         alphaOverlayExponent = -alphaOverlayExponent;
         alphaOverlayImage = 1-alphaOverlayImage;
      end
      alphaOverlayImage = alphaOverlayImage.^alphaOverlayExponent;
      alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
      %mask the alpha overlays with its own mask and the overlay mask
      alphaOverlayImage = alphaOverlayImage.*overlayMasks(:,:,:,iOverlay);
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alphaOverlayImage.*alphaOverlayMasks(:,:,:,iOverlay),[1 1 3]); 
    end   

    % Rescale current overlays.
    overlays.cmap(:,:,iOverlay) = viewGet(thisView,'overlayCmap',curOverlays(iOverlay));
    switch(viewGet(thisView,'overlayCtype',curOverlays(iOverlay)))
      case 'normal'
        overlays.colorRange(iOverlay,:) = viewGet(thisView,'overlayColorRange',curOverlays(iOverlay));
      case 'setRangeToMax'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(alphaOverlayMasks(:,:,:,iOverlay)))
          overlays.colorRange(1) = max(overlays.colorRange(1),min(overlays.overlayIm(alphaOverlayMasks(:,:,:,iOverlay))));
          overlays.colorRange(2) = min(max(overlays.overlayIm(alphaOverlayMasks(:,:,:,iOverlay))),overlays.colorRange(2));
        end
      case 'setRangeToMaxAroundZero'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(alphaOverlayMasks(:,:,:,iOverlay)))
          maxval = max(abs(overlays.overlayIm(alphaOverlayMasks(:,:,:,iOverlay))));
          overlays.colorRange(iOverlay,1) = -maxval;
          overlays.colorRange(iOverlay,2) = maxval;
        end
      case 'setRangeToMaxAcrossSlices'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum,scan); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum,scan);
        if ~isempty(minOverlay)
          overlays.colorRange(iOverlay,1) = max(minOverlay,overlays.colorRange(1));
        end
        if ~isempty(maxOverlay)
          overlays.colorRange(iOverlay,2) = min(maxOverlay,overlays.colorRange(2));
        end
      case 'setRangeToMaxAcrossSlicesAndScans'
        overlays.colorRange = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum);
        if ~isempty(minOverlay)
          overlays.colorRange(iOverlay,1) = max(minOverlay,overlays.colorRange(1));
        end
        if ~isempty(maxOverlay)
          overlays.colorRange(iOverlay,2) = min(maxOverlay,overlays.colorRange(2));
        end
    end
    overlays.RGB(:,:,:,iOverlay) = rescale2rgb(overlayImages(:,:,:,iOverlay),overlays.cmap(:,:,iOverlay),overlays.colorRange(iOverlay,:));
  end
  %alpha based on the alpha slider
  if nargin == 4
    alpha = viewGet(thisView,'alpha');
  end
  overlays.alphaMap = alpha*any(overlays.alphaMaps,4);
  %overlays.alphaMap = repmat(alpha*ones(imageDims,mrGetPref('defaultPrecision')),[1 1 3]); 

else
  overlays.coords = [];
  overlays.overlayIm = [];
  overlays.RGB = [];
end

