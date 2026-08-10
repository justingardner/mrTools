% function outputOverlay = keepKlargestConnectedRegions(overlay, <k>, <connectivity>)
%
%    Keeps only the k largest visible 3D connected clusters in (clipped) overlay (default k = 1)
%    If calling this function from combineTransformOverlays,
%    the "clip" or "alphaClip" checkbox must be checked.
%    For optional argument <connectivity> , see bwconncomp's help
%    options are 6, 18, 26 (default = 6)
%
%   author: julien besle (22/07/2020)


function outputOverlay = keepKlargestConnectedClusters(overlay,k, connectivity)

if ieNotDefined('k')
   k=1;
end
if ieNotDefined('connectivity')
   connectivity=6;
end

outputOverlay = nan(size(overlay));

%find connected clusters
cc = bwconncomp(~isnan(overlay),connectivity);
numPixels = cellfun(@numel,cc.PixelIdxList);
% sort them by descreasing size
[~,sizeIndex] = sort(numPixels,'descend');
for iCluster = 1:k
  outputOverlay(cc.PixelIdxList{sizeIndex(iCluster)}) = overlay(cc.PixelIdxList{sizeIndex(iCluster)});
end

