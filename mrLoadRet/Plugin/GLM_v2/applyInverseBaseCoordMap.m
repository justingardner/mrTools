% function volumeData = applyInverseBaseCoordMap(surf2volumeMap,volumeDims,surfData)
%
%   Transforms data from surface or flattened cortical patch space
%   to volume space according to the mapping in surf2volumeMap.
%   The mapping must first be computed using function inverseBaseCoordMap:
%     (e.g. surf2volumeMap = inverseBaseCoordMap(baseCoordsMap,volumeDims,<xform>) )
%
%   Taken out of combineTransformOverlays.m (22/07/2020)
%

function volumeData = applyInverseBaseCoordMap(surf2volumeMap,volumeDims,surfData)

hWaitBar = mrWaitBar(-inf,'(applyInverseBaseCoordMap) Converting from surface to volume');

volumeData = zeros(volumeDims);
datapoints = zeros(volumeDims);

% first find the longest non-zero row of the sparse flat2volumeMap matrix,
% as it is often much longer than the other rows and so time can be saved by treating it differently
longestRow = find(surf2volumeMap(:,end));
for iRow = longestRow % in case there are several such rows (unlikely)
  volumeData(iRow) = volumeData(iRow) + sum(surfData(surf2volumeMap(iRow,:)));
  datapoints(iRow) = size(surf2volumeMap,2);
end
surf2volumeMap(longestRow,:) = 0;
surf2volumeMap(:,sum(surf2volumeMap>0)==0) = [];

% now do the rest colum by column
maxInstances = size(surf2volumeMap,2);
for i=1:maxInstances
  mrWaitBar( i/maxInstances, hWaitBar);
  thisBaseCoordsMap = surf2volumeMap(:,i);
  newData = surfData(thisBaseCoordsMap(logical(thisBaseCoordsMap)));
  indices = find(thisBaseCoordsMap);
  notNaN = ~isnan(newData);
  volumeData(indices(notNaN)) = volumeData(indices(notNaN)) + newData(notNaN);
  datapoints(indices(notNaN)) = datapoints(indices(notNaN)) + 1;
end
datapoints = reshape(datapoints,volumeDims);
volumeData = volumeData ./datapoints;
mrCloseDlg(hWaitBar);