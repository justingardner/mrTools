% function surf2volumeMap = inverseBaseCoordMap(coordsMap,volumeDims,<xform>)
%
%   Computes a correspondence map between each voxel in a given volume and corresponding
%   voxels of a flat map or surface and outputs it as a sparse matrix
%   By default this is computed for the volume on which the flat map or surface was based.
%   If an xform (rotation matrix) is provided, it is computed for the corresponding volume.
%   In either case, volumeDims gives the dimensions of the destination volume (after rotation).
%
%   To transform data from flat space to volume space, call :
%      volumeData = applyInverseBaseCoordMap(surf2volumeMap,volumeDims,surfData)
%
%   Taken out of combineTransformOverlays.m (22/07/2020)


function surf2volumeMap = inverseBaseCoordMap(coordsMap,volumeDims,xform)

if ieNotDefined('xform')
  xform = eye(4);
end
coordsMap = reshape(coordsMap,numel(coordsMap)/3,3);
coordsMap = (xform*[coordsMap';ones(1,size(coordsMap,1))])';
coordsMap = coordsMap(:,1:3);
coordsMap(all(~coordsMap,2),:)=NaN;
coordsMap = round(coordsMap);
coordsMap(any(coordsMap>repmat(volumeDims,size(coordsMap,1),1)|coordsMap<1,2),:)=NaN;
% convert volume coordinates to linear indices for manipulation ease
volIndexMap = sub2ind(volumeDims, coordsMap(:,1), coordsMap(:,2), coordsMap(:,3));
clearvars('coordsMap'); % save memory

% now make a coordinate map of which volume index correspond to which surface voxels/vertices indices
% (there will be several maps because each volume voxels might correspond to several surface voxel/vertex indices)

% first find the maximum number of surface points corresponding to a single volume voxel (this is modified from function 'unique')
% sort volume indices
[sortedVolIndices,whichSurfIndices] = sort(volIndexMap);
whichSurfIndices(isnan(sortedVolIndices)) = []; % remove NaNs
sortedVolIndices(isnan(sortedVolIndices)) = []; % remove NaNs
nSurfVoxels = numel(sortedVolIndices);
% find the first instance of each unique index (except the very first)
firstInstances = sortedVolIndices(1:end-1) ~= sortedVolIndices(2:end);
firstInstances = [true;firstInstances];
% compute the number of instances for each unique volume index
% (= number of different base indices for each unique volume index)
numberInstances = diff(find([firstInstances;true]));
maxInstances = max(numberInstances);
% number each instance of a unique index from 1 to number of instances
instanceIndices = ones(nSurfVoxels,1);
firstInstances(1) = false;
instanceIndices(firstInstances) = -numberInstances(1:end-1)+1;
instanceIndices = cumsum(instanceIndices);
% fill the sparse matrix
surf2volumeMap = sparse(sortedVolIndices,instanceIndices,whichSurfIndices,prod(volumeDims),maxInstances);
