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
datapoints=zeros(prod(volumeDims),1);
hWaitBar = mrWaitBar(-inf,'(applyInverseBaseCoordMap) Converting from flat to volume');
maxInstances = size(flat2volumeMap,2);
for i=1:maxInstances
  mrWaitBar( i/maxInstances, hWaitBar);
  thisBaseCoordsMap = full(flat2volumeMap(:,i));
  volumeData(logical(thisBaseCoordsMap)) = volumeData(logical(thisBaseCoordsMap)) + flatData(thisBaseCoordsMap(logical(thisBaseCoordsMap)));
  datapoints = datapoints+logical(thisBaseCoordsMap);
end
datapoints = reshape(datapoints,volumeDims);
volumeData = volumeData ./datapoints;
mrCloseDlg(hWaitBar);