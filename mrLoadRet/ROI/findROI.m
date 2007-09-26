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

% get roi coordinates from the cache
roiCache = viewGet(view,'ROICache',viewGet(view,'currentROI'));
if isempty(roiCache)
  disp(sprintf('(findROI) ROI cache is emtpy (probably you are not viewing ROIs)'));
  return
end

if isempty(roiCache.s)
  msgbox(sprintf('ROI %s has no voxels in the current anatomy',viewGet(view,'roiName',roiNum)));
  return
end

% get current slice
curSlice = viewGet(view,'curSlice');

% find the closest slice
distanceToCurrentSlice = abs(roiCache.s-curSlice);
closestSlice = roiCache.s(first(find(min(distanceToCurrentSlice)==distanceToCurrentSlice)));

% set the slice
mlrGuiSet(view.viewNum,'slice',closestSlice);
refreshMLRDisplay(view.viewNum);

