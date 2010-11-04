% getNewSpaceOverlay.m
%
%      usage: newOverlayData = (overlayData, xform, newXCoords, newYCoords, newZCoords, interpMethod)
%
%         by: julien besle, loosely based on mrExport2SR by eli merriam
%       date: 11/05/2010
%    purpose: computes overlay values in a new space
%        $Id:$	
%

function newOverlayData = getNewSpaceOverlay(overlayData, xform, newXCoords, newYCoords, newZCoords, interpMethod)

% interpMethod and interpExtrapVal are used by calls to interp3 when
% extracting slices from the base and overlay arrays.
if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end
if isempty(interpMethod)
    interpMethod = 'linear';
end
interpExtrapVal = NaN;

sliceDims(1) = size(newXCoords,1);
sliceDims(2) = size(newXCoords,2);
numPixels = sliceDims(1)*sliceDims(2);
newOverlayData = zeros(size(newXCoords));

hWaitbar = mrWaitbar(0,'Resampling overlay to new space');
% Compute new overlay data by base slice
nSlices = size(newXCoords,3);
for i_slice = 1:nSlices
    
   baseCoordsHomogeneous = [reshape(newXCoords(:,:,i_slice),1,numPixels); reshape(newYCoords(:,:,i_slice),1,numPixels); reshape(newZCoords(:,:,i_slice),1,numPixels); ones(1,numPixels)];
   %baseCoords = reshape(baseCoordsHomogeneous(1:3,:)',[sliceDims 3]);

   % Transform coordinates
   overlayCoordsHomogeneous = xform * baseCoordsHomogeneous;
   overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[sliceDims 3]);

   % Extract slice from current overlay.
   if ~isempty(overlayCoords) 
     newOverlayData(:,:,i_slice) = interp3(overlayData,overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3), interpMethod,interpExtrapVal);
   end
   mrWaitbar(i_slice/nSlices,hWaitbar);
end
mrCloseDlg(hWaitbar);
