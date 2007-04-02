% getRoiCoordinates.m
%
%      usage: scanCoords = getRoiCoordinates(view,roiNum,<scanNum>,<groupNum>)
%         by: justin gardner
%       date: 04/02/07
%    purpose: get roi coordinates in scan coordinates
%
function scanCoords = getRoiCoordinates(view,roiNum,scanNum,groupNum)

% check arguments
if ~any(nargin == [2 3 4])
  help getRoiCoordinates
  return
end

% get group and scan
if ieNotDefined('groupNum')
  groupNum = viewGet(view,'currentGroup');
end
if ieNotDefined('scanNum')
  scanNum = viewGet(view,'currentScan');
end

% get the scan transforms
scanXform = viewGet(view,'scanXform',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);

% get the roi transforms
roiXform = viewGet(view,'roiXform',roiNum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);
roiCoords = viewGet(view,'roiCoords',roiNum);

% Use xformROI to supersample the coordinates
scanCoords = round(xformROIcoords(roiCoords,inv(scanXform)*roiXform,roiVoxelSize,scanVoxelSize));

% return the unique ones
scanCoords = unique(scanCoords','rows')';
scanCoords = scanCoords(1:3,:);
