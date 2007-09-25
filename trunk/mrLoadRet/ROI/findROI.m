% findROI.m
%
%        $Id$
%      usage: findROI()
%         by: justin gardner
%       date: 09/25/07
%    purpose: looks for an ROI in the base anatomy (is called form mrLoadRetGUI)
%
function retval = findROI(view)

retval = [];

% check arguments
if ~any(nargin == [1])
  help findROI
  return
end

if ~isempty(viewGet(view,'baseCoordMap'))
  disp(sprintf('(findROI) Base anatomy cannot be flat to look for an ROI'));
  return
end
% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

% get the current scan number
scanNum = viewGet(view,'curScan');
% get roi coordinates
coords = getROICoordinates(view,roiNum,scanNum);

if isempty(coords)
  msgbox(sprintf('ROI %s has no voxels in the current anatomy',viewGet(view,'roiName',roiNum)));
  return
end

% get current slice
curSlice = viewGet(view,'curSlice');

% go find the ROI
sliceIndex = viewGet(view,'baseSliceIndex');

% find the closest slice
distanceToCurrentSlice = abs(coords(sliceIndex,:)-curSlice);
closestSlice = coords(sliceIndex,first(find(min(distanceToCurrentSlice)==distanceToCurrentSlice)));

% set the slice
mlrGuiSet(view.viewNum,'slice',closestSlice);
refreshMLRDisplay(view.viewNum);

