% function volumeData = applyInverseBaseCoordMap(flat2volumeMap,flatData)
%
%   Transforms data from flattened cortical patch space to volume space
%   according to the mapping in flat2volumeMap.
%   The mapping must first be computed using function inverseBaseCoordMap:
%     (e.g. flat2volumeMap = inverseBaseCoordMap(baseCoordsMap,volumeDims,<xform>) )
%
%   Taken out of combineTransformOverlays.m (22/07/2020)
%

function volumeData = applyInverseBaseCoordMap(flat2volumeMap,volumeDims,flatData)

volumeData = zeros(volumeDims);
datapoints = zeros(volumeDims);
hWaitBar = mrWaitBar(-inf,'(applyInverseBaseCoordMap) Converting from surface to volume');
maxInstances = size(flat2volumeMap,2);
for i=1:maxInstances
  mrWaitBar( i/maxInstances, hWaitBar);
  thisBaseCoordsMap = full(flat2volumeMap(:,i));
  newData = flatData(thisBaseCoordsMap(logical(thisBaseCoordsMap)));
  indices = find(thisBaseCoordsMap);
  notNaN = ~isnan(newData);
  volumeData(indices(notNaN)) = volumeData(indices(notNaN)) + newData(notNaN);
  datapoints(indices(notNaN)) = datapoints(indices(notNaN)) + 1;
end
datapoints = reshape(datapoints,volumeDims);
volumeData = volumeData ./datapoints;
mrCloseDlg(hWaitBar);