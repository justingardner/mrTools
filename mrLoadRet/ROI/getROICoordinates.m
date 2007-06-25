% getROICoordinates.m
%
%      usage: scanCoords = getROICoordinates(view,roiNum,<scanNum>,<groupNum>)
%         by: justin gardner
%       date: 04/02/07
%    purpose: get roi coordinates in scan coordinates
%             if scanNum is 0, then will compute in base
%             coordinates. 
%             if roinum is a structure, works on the structure
%             rather than the roinum
%             if roinum is a string, will load the roi from
%             the directory
function scanCoords = getROICoordinates(view,roiNum,scanNum,groupNum)

scanCoords = [];
% check arguments
if ~any(nargin == [2 3 4])
  help getROICoordinates
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
  % use base xform if scanNum == 0
  view = viewSet(view,'curGroup',groupNum);
  scanXform = viewGet(view,'baseXform');
  scanVoxelSize = viewGet(view,'baseVoxelSize');
end  

% if roiNum is a string try to load it
if isstr(roiNum)
  roiname = fullfile(viewGet(view,'roidir'),fixBadChars(roiNum));
  roiname = sprintf('%s.mat',stripext(roiname));
  if ~isfile(roiname)
    disp(sprintf('(getROICoordinates) Could not find roi %s',roiNum));
    return
  end
  r = load(roiname);
  f = fieldnames(r);
  if length(f) >= 1
    roiNum = r.(f{1});
  else
    disp(sprintf('(getROICoordinates) Could not find roi in file %s',roiNum));
    return
  end
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

% check stuff
if (isempty(scanXform)) 
  disp(sprintf('(getRoiCoordinates) scanXform for %s:%i is empty',viewGet(view,'groupName',groupNum),scanNum));
  return
end
if (isempty(roiXform)) 
  disp(sprintf('(getRoiCoordinates) roiXform is empty'));
  return
end

% make sure we have normalized coordinates
if (size(roiCoords,1)==3)
  roiCoords(4,:) = 1;
end

% Use xformROI to supersample the coordinates
scanCoords = round(xformROIcoords(roiCoords,inv(scanXform)*roiXform,roiVoxelSize,scanVoxelSize));

% return the unique ones
scanCoords = unique(scanCoords','rows')';
scanCoords = scanCoords(1:3,:);

if ~isempty(scanCoords) && (scanNum ~= 0)

  % check scan dimensions
  scanDims = viewGet(view,'dims',scanNum,groupNum);

  % make sure we are inside scan dimensions
  xCheck = (scanCoords(1,:) >= 1) & (scanCoords(1,:) <= scanDims(1));
  yCheck = (scanCoords(2,:) >= 1) & (scanCoords(2,:) <= scanDims(2));
  sCheck = (scanCoords(3,:) >= 1) & (scanCoords(3,:) <= scanDims(3));

  % only return ones that are in bounds
  scanCoords = scanCoords(:,find(xCheck & yCheck & sCheck));
end

