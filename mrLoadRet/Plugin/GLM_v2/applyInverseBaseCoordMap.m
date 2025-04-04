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

nVolumes = size(surfData,4);
volumeData = zeros([volumeDims nVolumes]);

for iVolume = 1:nVolumes
  thisSurfData = surfData(:,:,:,iVolume);
  thisVolumeData = zeros(volumeDims);
  datapoints = zeros(volumeDims);
  thisSurf2volMap = surf2volumeMap;

  % first find the longest non-zero row of the sparse flat2volumeMap matrix,
  % as it is often much longer than the other rows and so time can be saved by treating it differently
  longestRow = find(thisSurf2volMap(:,end));
  for iRow = longestRow % in case there are several such rows (unlikely)
    thisVolumeData(iRow) = thisVolumeData(iRow) + sum(thisSurfData(thisSurf2volMap(iRow,:)));
    datapoints(iRow) = size(thisSurf2volMap,2);
  end
  thisSurf2volMap(longestRow,:) = 0;
  thisSurf2volMap(:,sum(thisSurf2volMap>0)==0) = [];

  % now do the rest colum by column
  maxInstances = size(thisSurf2volMap,2);
  for i=1:maxInstances
    mrWaitBar( ((iVolume-1)*maxInstances+i) / (maxInstances*nVolumes), hWaitBar);
    thisBaseCoordsMap = thisSurf2volMap(:,i);
    newData = thisSurfData(thisBaseCoordsMap(logical(thisBaseCoordsMap)));
    indices = find(thisBaseCoordsMap);
    notNaN = ~isnan(newData);
    thisVolumeData(indices(notNaN)) = thisVolumeData(indices(notNaN)) + newData(notNaN);
    datapoints(indices(notNaN)) = datapoints(indices(notNaN)) + 1;
  end
  datapoints = reshape(datapoints,volumeDims);
  volumeData(:,:,:,iVolume) = thisVolumeData ./datapoints;
end

mrCloseDlg(hWaitBar);
