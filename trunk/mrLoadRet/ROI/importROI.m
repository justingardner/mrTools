% importROI.m
%
%      usage: importROI(view,pathStr)
%         by: justin gardner
%       date: 03/16/07
%    purpose: import roi from mrLoadRet 3.1 version to 4.5
%
function view = importROI(view,pathStr)

% check arguments
if ~any(nargin == [1 2])
  help importROI
  return
end

mrGlobals;
% Complete pathStr
if ieNotDefined('pathStr')
  % start in an roi directory
  %startPathStr = fullfile(viewGet(view,'viewType'),'ROIs');
  startPathStr = mrGetPref('importROIPath');
  if isempty(startPathStr)
    startPathStr = 'Inplane/ROIs';
  end
  if ~isdir(startPathStr),startPathStr='';,end
  % get the user defined path
  pathStr = getPathStrDialog(startPathStr,'Choose roi files to import','*.mat','on');
end
if isempty(pathStr),disp('No ROI selected');,return,end
mrSetPref('importROIPath',fileparts(pathStr{1}));

% get some info
baseNum = viewGet(view,'currentBase');
xform = viewGet(view,'basexform',baseNum);
voxelSize = viewGet(view,'baseVoxelSize',baseNum);
baseDims = viewGet(view,'baseDims',baseNum);

for roinum = 1:length(pathStr)
  % try to load the roi
  l = load(pathStr{roinum});
  if isfield(l,'ROI')
    clear ROI;
    ROI.name = l.ROI.name;
    ROI.viewType = view.viewType;
    ROI.color = l.ROI.color;
    if isfield(l.ROI,'viewType') && strcmp(l.ROI.viewType,'Gray')
      % not sure why gray rois are different from inplane but
      % this seems to work in conversion
      ROI.coords(1,:) = l.ROI.coords(3,:);
      ROI.coords(2,:) = baseDims(2)-l.ROI.coords(2,:)+1;
      ROI.coords(3,:) = baseDims(3)-l.ROI.coords(1,:)+1;
    else
      % there is just an x/y flip for the inplane ROIs
      ROI.coords(1,:) = l.ROI.coords(2,:);
      ROI.coords(2,:) = l.ROI.coords(1,:);
      ROI.coords(3,:) = l.ROI.coords(3,:);
    end
    ROI.coords(4,:) = 1;
    ROI.xform = xform;
    ROI.voxelSize = voxelSize;
    ROI.date = datestr(now);
    % Add it to the view
    view = viewSet(view,'newROI',ROI);
    %ROI.coords
  end
end

if exist('ROI','var')
  ROInum = viewGet(view,'ROInum',ROI.name);
  if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
    view = viewSet(view,'prevROIcoords',[]);
  end
  refreshMLRDisplay(viewGet(view,'viewNum'));
end
