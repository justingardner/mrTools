% getRoiCoordinates.m
%
%      usage: scanCoords = getRoiCoordinates(view,roiNum,<scanNum>,<groupNum>)
%         by: justin gardner
%       date: 04/02/07
%    purpose: get roi coordinates in scan coordinates
%             if scanNum is 0, then will compute in base
%             coordinates. 
%             if roinum is a structure, works on the structure
%             rather than the roinum
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
if scanNum
  scanXform = viewGet(view,'scanXform',scanNum,groupNum);
  scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
else
  view = viewSet(view,'curGroup',groupNum);
  scanXform = viewGet(view,'baseXform');
  scanVoxelSize = viewGet(view,'baseVoxelSize');
end  

% get the roi transforms
if isstruct(roiNum)
  roiXform = roiNum.xform;
  roiVoxelSize = roiNum.voxelSize;
  roiCoords = roiNum.coords;
else
  roiXform = viewGet(view,'roiXform',roiNum);
  roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);
  roiCoords = viewGet(view,'roiCoords',roiNum);
end
% Use xformROI to supersample the coordinates
scanCoords = round(xformROIcoords(roiCoords,inv(scanXform)*roiXform,roiVoxelSize,scanVoxelSize));

% return the unique ones
scanCoords = unique(scanCoords','rows')';
if ~isempty(scanCoords)
  scanCoords = scanCoords(1:3,:);
end