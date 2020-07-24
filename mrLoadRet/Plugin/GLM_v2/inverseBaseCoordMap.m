% function flat2volumeMap = inverseBaseCoordMap(baseCoordsMap,volumeDims,<xform>)
%
%   Computes a correspondence map between each voxel in a given volume and corresponding
%   voxels of a flat map and outputs it as a sparse matrix
%   By default this is computed for the base space on which the flat map is based.
%   If an xform (rotation matrix) is provided, it is computed for the corresponding volume.
%   In either case, volumeDims gives the dimensions of the destination volume (after rotation).
%
%   To transform data from flat space to volume space, call :
%       volumeData = applyInverseBaseCoordMap(flat2volumeMap,flatData)
%
%   Taken out of combineTransformOverlays.m (22/07/2020)


function flat2volumeMap = inverseBaseCoordMap(baseCoordsMap,volumeDims,xform)

if ieNotDefined('xform')
  xform = eye(4);
end
baseCoordsMap = reshape(baseCoordsMap,numel(baseCoordsMap)/3,3);
overlayCoordsMap = (xform*[baseCoordsMap';ones(1,size(baseCoordsMap,1))])';
overlayCoordsMap = overlayCoordsMap(:,1:3);
overlayCoordsMap(all(~overlayCoordsMap,2),:)=NaN;
overlayCoordsMap = round(overlayCoordsMap);
overlayCoordsMap(any(overlayCoordsMap>repmat(volumeDims,size(overlayCoordsMap,1),1)|overlayCoordsMap<1,2),:)=NaN;
%convert overlay coordinates to overlay indices for manipulation ease
overlayIndexMap = sub2ind(volumeDims, overlayCoordsMap(:,1), overlayCoordsMap(:,2), overlayCoordsMap(:,3));

%now make a coordinate map of which base map voxels each overlay index corresponds to
%(there will be several maps because each overlay voxels might correspond to several base voxels)

% %       %METHOD 1
% %       %sort base indices
% %       [sortedOverlayIndices,whichBaseIndices] = sort(overlayIndexMap);
% %       %remove NaNs (which should be at the end of the vector)
% %       whichBaseIndices(isnan(sortedOverlayIndices))=[];
% %       sortedOverlayIndices(isnan(sortedOverlayIndices))=[];
% %       %find the first instance of each unique index
% %       firstInstances = sortedIndices(1:end-1) ~= sortedIndices(2:end);
% %       firstInstances = [true;firstInstances];
% %       %get the unique overlay indices
% %       uniqueOverlayIndices = sortedOverlayIndices(firstInstances);
% %       %compute the number of instances for each  unique overlay index (= number
% %       %of base different indices for each unique overlay index)
% %       numberInstances = diff(find([firstInstances;true]));
% %       maxInstances = max(numberInstances);
% %       baseCoordsOverlay2{iScan} = sparse(prod(scanDims),maxInstances);
% %       hWaitBar = mrWaitBar(-inf,'(combineTransformOverlays) Creating base coordinates overlay map for scan');
% %       %for each unique overlay index, find all the corresponding base indices
% %       for i = 1:length(uniqueOverlayIndices)
% %         mrWaitBar( i/length(uniqueOverlayIndices), hWaitBar);
% %         theseBaseIndices = whichBaseIndices(sortedOverlayIndices==uniqueOverlayIndices(i));
% %         baseCoordsOverlay2{iScan}(uniqueOverlayIndices(i),1:length(theseBaseIndices))=theseBaseIndices';
% %       end
% %       mrCloseDlg(hWaitBar);

%METHOD 2 (faster)
%first find the maximum number of base voxels corresponding to a single overlay voxel (this is modified from function 'unique')
%sort base non-NaN indices
sortedIndices = sort(overlayIndexMap(~isnan(overlayIndexMap)));
%find the first instance of each unique index
firstInstances = sortedIndices(1:end-1) ~= sortedIndices(2:end);
firstInstances = [true;firstInstances];
%compute the number of instances for each unique overlay index
%(= number of base different indices for each unique overlay index)
numberInstances = diff(find([firstInstances;true]));
maxInstances = max(numberInstances);
flat2volumeMap = sparse(prod(volumeDims),maxInstances);
%Now for each set of unique overlay indices, find the corresponding base indices
hWaitBar = mrWaitBar(-inf,'(inverseBaseCoordMap) Computing inverse baseCoordMap');
for i=1:maxInstances
  mrWaitBar( i/maxInstances, hWaitBar);
  %find set of unique instances of overlay indices
  [uniqueOverlayIndices, whichBaseIndices]= unique(overlayIndexMap);
  %remove NaNs
  whichBaseIndices(isnan(uniqueOverlayIndices))=[];
  uniqueOverlayIndices(isnan(uniqueOverlayIndices))=[];
  %for each overlay voxel found, set the corresponding base index
  flat2volumeMap(uniqueOverlayIndices,i)=whichBaseIndices;
  %remove instances that were found from the overlay index map before going through the loop again
  overlayIndexMap(whichBaseIndices)=NaN;
end
mrCloseDlg(hWaitBar);
