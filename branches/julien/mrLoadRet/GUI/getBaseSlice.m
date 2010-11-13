% getBaseSlice.m
%
%        $Id$ 
%      usage: [baseIm baseCoords baseCoordsHomogeneous] = getBaseSlice(view, sliceNum, sliceIndex, rotate, baseNum, baseType)
%         by: justin gardner
%       date: 03/05/10
%    purpose: DJH - pulled out from refreshMLRDisplay for use in other functions -jg
%             extracts base image and corresponding coordinates
function [baseIm,baseCoords,baseCoordsHomogeneous] = getBaseSlice(view,sliceNum,sliceIndex,rotate,baseNum,baseType)

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
  baseCoordMap = viewGet(view,'baseCoordMap',baseNum);
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
