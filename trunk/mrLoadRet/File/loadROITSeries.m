% loadROITSeries.m
%
%      usage: loadROITSeries(view,scanNum,<roiname>)
%         by: justin gardner
%       date: 03/22/07
%    purpose: load the time series for a roi
%
function rois = loadROITSeries(view,scanNum,roiname)

rois = {};

% check arguments
if ~any(nargin == [1 2 3])
  help loadROITSeries
  return
end

% get group and scan
groupNum = viewGet(view,'currentGroup');
if ~exist('scanNum','var')
  scanNum = viewGet(view,'currentScan');
end

% if there is no roi, ask the user to select
if ~exist('roiname','var')
  roiname = getPathStrDialog(viewGet(view,'roiDir'),'Choose one or more ROIs','*.mat','on');
elseif isstr(roiname)
  %make into a cell array
  tmp = roiname;
  roiname = {};
  roiname{1} = tmp;
end

% load the rois in turn
for roinum = 1:length(roiname)
  % check for file
  if ~isfile(roiname{roinum})
    disp(sprintf('(loadROITSeries) Could not find roi %s',roiname{roinum}));
  end
  % load the roi
  roi = load(roiname{roinum});
  roiFieldnames = fieldnames(roi);
  % get all the rois
  for roinum = 1:length(roiFieldnames)
    rois{end+1} = roi.(roiFieldnames{roinum});
    % convert to scan coordinates
    rois{end}.scanCoords = getROICoords(view,groupNum,scanNum,rois{end});
  end
end

  

function scanCoords = getROICoords(view,groupNum,scanNum,roi)

% get the scan transforms
scanXform = viewGet(view,'scanXform',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);

% Use xformROI to supersample the coordinates
scanCoords = round(xformROIcoords(roi.coords,inv(scanXform)*roi.xform,roi.voxelSize,scanVoxelSize));

% return the unique ones
scanCoords = unique(scanCoords','rows')';

keyboard