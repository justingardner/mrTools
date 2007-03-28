% loadROITSeries.m
%
%      usage: rois = loadROITSeries(view,<roiname>,<scanNum>,<groupNum>)
%         by: justin gardner
%       date: 03/22/07
%    purpose: load the time series for a roi, without roiname
%             specified, brings up selection dialog. roiname
%             may be a cell array of rois. scanNum and groupNum
%             default to current scan/group
%
%        e.g.:
%
%v = newView('Volume')
%rois = loadROITseries(v);
%
function rois = loadROITSeries(view,roiname,scanNum,groupNum);

rois = {};

% check arguments
if ~any(nargin == [1 2 3 4])
  help loadROITSeries
  return
end

% get the roi directory
roidir = viewGet(view,'roidir');

% get group and scan
if ieNotDefined('groupNum')
  groupNum = viewGet(view,'currentGroup');
end
if ieNotDefined('scanNum')
  scanNum = viewGet(view,'currentScan');
end

% set the current group
view = viewSet(view,'currentGroup',groupNum);

% if there is no roi, ask the user to select
if ieNotDefined('roiname')
  roiname = getPathStrDialog(viewGet(view,'roiDir'),'Choose one or more ROIs','*.mat','on');
elseif isstr(roiname)
  %make into a cell array
  tmp = roiname;
  roiname = {};
  roiname{1} = tmp;
end

% load the rois in turn
for roinum = 1:length(roiname)
  % see if we have to past roi directory on
  if ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
    roiname{roinum} = fullfile(roidir,stripext(roiname{roinum}));
  end
  % check for file
  if ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
    disp(sprintf('(loadROITSeries) Could not find roi %s',roiname{roinum}));
  else
    % load the roi
    roi = load(roiname{roinum});
    roiFieldnames = fieldnames(roi);
    % get all the rois
    for roinum = 1:length(roiFieldnames)
      rois{end+1} = roi.(roiFieldnames{roinum});
      % convert to scan coordinates
      rois{end}.scanCoords = getROICoords(view,groupNum,scanNum,rois{end});
      % get x y and s in array form
      x = rois{end}.scanCoords(1,:);
      y = rois{end}.scanCoords(2,:);
      s = rois{end}.scanCoords(3,:);
      % set the n
      rois{end}.n = length(x);
      % load the tseries, voxel-by-voxel
      disppercent(-inf,sprintf('Loading tSeries for %s',rois{end}.name));
      % for now we always load by block, but if memory is an issue, we can
      % switch this if statement and load voxels indiviudally from file
      if 0
	% load each voxel time series indiviudally
	for voxnum = 1:rois{end}.n
	  rois{end}.tSeries(voxnum,:) = squeeze(loadTSeries(view,scanNum,s(voxnum),[],x(voxnum),y(voxnum)));
	  disppercent(voxnum/rois{end}.n);
	end
      else
	% load the whole time series as a block (i.e. a block including the min and max voxels)
	% this is usually faster then going back and loading each time series individually
	% but is more inefficient with memory
	tSeriesBlock = squeeze(loadTSeries(view,scanNum,[min(s) max(s)],[],[min(x) max(x)],[min(y) max(y)]));
	% now go through and pick out the voxels that we need.
	for voxnum = 1:rois{end}.n
	  rois{end}.tSeries(voxnum,:) = squeeze(tSeriesBlock(x(voxnum)-min(x)+1,y(voxnum)-min(y)+1,s(voxnum)-min(s)+1,:));
	  disppercent(voxnum/rois{end}.n);
	end
	clear tSeriesBlock;
      end
      disppercent(inf);
    end
  end
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert the roi coordinates into this scans coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scanCoords = getROICoords(view,groupNum,scanNum,roi)

% get the scan transforms
scanXform = viewGet(view,'scanXform',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);

% Use xformROI to supersample the coordinates
scanCoords = round(xformROIcoords(roi.coords,inv(scanXform)*roi.xform,roi.voxelSize,scanVoxelSize));

% return the unique ones
scanCoords = unique(scanCoords','rows')';

