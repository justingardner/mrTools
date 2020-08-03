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


function flat2volumeMap = inverseBaseCoordMap(coordsMap,volumeDims,xform)

if ieNotDefined('xform')
  xform = eye(4);
end
coordsMap = reshape(coordsMap,numel(coordsMap)/3,3);
coordsMap = (xform*[coordsMap';ones(1,size(coordsMap,1))])';
coordsMap = coordsMap(:,1:3);
coordsMap(all(~coordsMap,2),:)=NaN;
coordsMap = round(coordsMap);
coordsMap(any(coordsMap>repmat(volumeDims,size(coordsMap,1),1)|coordsMap<1,2),:)=NaN;
% convert overlay coordinates to overlay indices for manipulation ease
overlayIndexMap = sub2ind(volumeDims, coordsMap(:,1), coordsMap(:,2), coordsMap(:,3));
clearvars('coordsMap'); % save memory

% now make a coordinate map of which base map voxels each overlay index corresponds to
% (there will be several maps because each overlay voxels might correspond to several base voxels)

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

% METHOD 2 (faster)
% first find the maximum number of base voxels corresponding to a single overlay voxel (this is modified from function 'unique')
% sort base non-NaN indices
sortedIndices = sort(overlayIndexMap(~isnan(overlayIndexMap)));
nBaseVoxels = numel(sortedIndices);
% find the first instance of each unique index
firstInstances = sortedIndices(1:end-1) ~= sortedIndices(2:end);
firstInstances = [true;firstInstances];
% compute the number of instances for each unique overlay index
% (= number of different base indices for each unique overlay index)
numberInstances = diff(find([firstInstances;true]));
maxInstances = max(numberInstances);
% Now for each set of unique overlay indices, find the corresponding base indices
hWaitBar = mrWaitBar(-inf,'(inverseBaseCoordMap) Inverting baseCoordMap');
% flat2volumeMap = sparse(prod(volumeDims),maxInstances); % I used to create the sparse mapping matrix,
% % but supposedly it's faster to first gather the matrix's indices and data and create it later (doesn't
% % make much of a difference though, presumably because I minimized the number of iteratiosnin the loop)
allUniqueOverlayIndices = zeros(nBaseVoxels,1);
allWhichBaseIndices = zeros(nBaseVoxels,1);
allInstances = zeros(nBaseVoxels,1);
n = 0;
for i=1:maxInstances
  mrWaitBar( i/maxInstances, hWaitBar);
  % find set of unique instances of overlay indices
  [uniqueOverlayIndices, whichBaseIndices]= unique(overlayIndexMap);
  % remove NaNs
  whichBaseIndices(isnan(uniqueOverlayIndices))=[];
  uniqueOverlayIndices(isnan(uniqueOverlayIndices))=[];
  nUniqueIndices = length(whichBaseIndices);
  % keep aside to fill sparse matrix later
  allUniqueOverlayIndices(n+(1:nUniqueIndices)) = uniqueOverlayIndices;
  allWhichBaseIndices(n+(1:nUniqueIndices)) = whichBaseIndices;
  allInstances(n+(1:nUniqueIndices)) = i;
%   % for each overlay voxel found, set the corresponding base index
%   flat2volumeMap(uniqueOverlayIndices,i)=whichBaseIndices;
  % remove instances that were found from the overlay index map before going through the loop again
  overlayIndexMap(whichBaseIndices)=NaN;
  n = n+nUniqueIndices;
end
% fill the sparse matrix
flat2volumeMap = sparse(allUniqueOverlayIndices,allInstances,allWhichBaseIndices,prod(volumeDims),maxInstances);
mrCloseDlg(hWaitBar);
