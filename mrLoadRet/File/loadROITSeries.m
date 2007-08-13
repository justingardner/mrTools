function rois = loadROITSeries(view,roiname,scanList,groupNum,varargin);
% loadROITSeries.m
%
%      usage: rois = loadROITSeries(view,<roiname>,<scanList>,<groupNum>)
%         by: justin gardner
%       date: 03/22/07
%    purpose: load the time series for a roi, without roiname
%             specified, brings up selection dialog. roiname
%             may be a cell array. scanList and groupNum
%             default to current scan/group. If roiname is a roi
%             struct instead of a name then it will use that roi
%             instead of loading the roi from disk. Also, if roiname
%             is a number or cell array of numbers then it will use
%             the corresponding ROI from the view.
%        e.g.:
%
% v = newView('Volume')
% rois = loadROITSeries(v,[],1,1);
%
%             to load the roi coordinates, but not the time series
%       
% v = newView('Volume')
% rois = loadROITSeries(v,[],1,1,'loadType=none');
%
% see also tseriesROI

rois = {};

% check arguments
if nargin < 1
  help loadROITSeries
  return
end

% evaluate other arguments
eval(evalargs(varargin));

% no view specified
if ieNotDefined('view')
  view = newView('Volume');
end

% get the roi directory
roidir = viewGet(view,'roidir');

% get group and scan
if ieNotDefined('groupNum')
  groupNum = viewGet(view,'currentGroup');
end
groupName = viewGet(view,'groupName',groupNum);
if ieNotDefined('scanList')
  scanList = viewGet(view,'currentScan');
end

% set the current group
view = viewSet(view,'currentGroup',groupNum);

% if there is no roi, ask the user to select
if ieNotDefined('roiname')
  roiname = getPathStrDialog(viewGet(view,'roiDir'),'Choose one or more ROIs','*.mat','on');
end

%make into a cell array
roiname = cellArray(roiname);

% set the way to load the time series
% possible values are 'vox' which loads each
% time series voxel by voxel (slow but less memory
% intensive), or 'block' which loads the block
% of the image around the ROI and then subselects
% the voxels needed (default--fast). Note, both
% block and vox will return the same voxel time
% series, they just differ in how they access the data
% from disk. Set to 'none' to not load the time series.
if ieNotDefined('loadType')
  loadType = 'block';
end

% load the rois in turn
for roinum = 1:length(roiname)
  % see if we have to paste roi directory on
  if isstr(roiname{roinum}) && ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
    roiname{roinum} = fullfile(roidir,stripext(roiname{roinum}));
  end
  % check for file
  if isstr(roiname{roinum}) && ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
    disp(sprintf('(loadROITSeries) Could not find roi %s',roiname{roinum}));
    dir(fullfile(roidir,'*.mat'))
  elseif isnumeric(roiname{roinum}) && ((roiname{roinum} < 1) || (roiname{roinum} > viewGet(view,'numberOfROIs')))
    disp(sprintf('(loadROITSeries) No ROI number %i (number of ROIs = %i)',roiname{roinum},viewGet(view,'numberOfROIs')));
  else
    % load the roi, if the name is actually a struct
    % then assume it is an roi struct. if it is a number choose
    % from a loaded roi
    if isstr(roiname{roinum})
      roi = load(roiname{roinum});
    elseif isnumeric(roiname{roinum})
      roi = viewGet(view,'roi',roiname{roinum});
    else
      roi.(fixBadChars(roiname{roinum}.name)) = roiname{roinum};
    end
    roiFieldnames = fieldnames(roi);
    % get all the rois
    for roinum = 1:length(roiFieldnames)
      for scanNum = 1:length(scanList)
        % get current scan number
        scanNum = scanList(scanNum);
        rois{end+1} = roi.(roiFieldnames{roinum});
        % set a field in the roi for which scan we are collecting from
        rois{end}.scanNum = scanNum;
        rois{end}.groupNum = groupNum;
        % convert to scan coordinates
        rois{end}.scanCoords = getROICoordinates(view,rois{end},scanNum,groupNum);
        % if there are no scanCoords then set to empty and continue
        if isempty(rois{end}.scanCoords)
          rois{end}.n = 0;
          rois{end}.tSeries = [];
          continue;
        end
        % get x y and s in array form
        x = rois{end}.scanCoords(1,:);
        y = rois{end}.scanCoords(2,:);
        s = rois{end}.scanCoords(3,:);
        % set the n
        rois{end}.n = length(x);
        % load the tseries, voxel-by-voxel
        disppercent(-inf,sprintf('Loading tSeries for %s from %s: %i',rois{end}.name,groupName,scanNum));
        % for now we always load by block, but if memory is an issue, we can
        % switch this if statement and load voxels indiviudally from file
        if strcmp(loadType,'vox')
          % load each voxel time series indiviudally
          for voxnum = 1:rois{end}.n
            rois{end}.tSeries(voxnum,:) = squeeze(loadTSeries(view,scanNum,s(voxnum),[],x(voxnum),y(voxnum)));
            disppercent(voxnum/rois{end}.n);
          end
	  disppercent(inf);
        elseif strcmp(loadType,'block');
          % load the whole time series as a block (i.e. a block including the min and max voxels)
          % this is usually faster then going back and loading each time series individually
          % but is more inefficient with memory
          tSeriesBlock = loadTSeries(view,scanNum,[min(s) max(s)],[],[min(x) max(x)],[min(y) max(y)]);
          % now go through and pick out the voxels that we need.
          for voxnum = 1:rois{end}.n
            rois{end}.tSeries(voxnum,:) = squeeze(tSeriesBlock(x(voxnum)-min(x)+1,y(voxnum)-min(y)+1,s(voxnum)-min(s)+1,:));
            disppercent(voxnum/rois{end}.n);
          end
          clear tSeriesBlock;
	  disppercent(inf);
        else
	  disppercent(inf);
	  disp(sprintf('(loadROITSeries) Not loading time series (loadType=%s)',loadType));
	end
      end
    end
  end
end
if length(rois) == 1
  rois = rois{1};
end

