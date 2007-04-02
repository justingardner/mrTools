% importROI.m
%
%      usage: importROI(view,pathStr)
%         by: justin gardner
%       date: 03/16/07
%    purpose: import pathStr
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
  if isfield(MLR.prefs,'importROIpath')
    startPathStr = MLR.prefs.importROIpath;
  else
    startPathStr = 'Inplane/ROIs';
  end
  if ~isdir(startPathStr),startPathStr='';,end
  % get the user defined path
  pathStr = getPathStrDialog(startPathStr,'Choose roi files to import','*.mat','on');
end
if isempty(pathStr),disp('No ROI selected');,return,end

% remember the path for next
MLR.prefs.importROIpath = fileparts(pathStr{1});

% get some info
baseNum = viewGet(view,'currentBase');
xform = viewGet(view,'basexform',baseNum);
voxelSize = viewGet(view,'baseVoxelSize',baseNum);

for roinum = 1:length(pathStr)
  % try to load the roi
  l = load(pathStr{roinum});
  if isfield(l,'ROI')
    clear ROI;
    ROI.name = l.ROI.name;
    ROI.viewType = view.viewType;
    ROI.color = l.ROI.color;
    ROI.coords(1,:) = l.ROI.coords(2,:);
    ROI.coords(2,:) = l.ROI.coords(1,:);
    ROI.coords(3,:) = l.ROI.coords(3,:);
    ROI.coords(4,:) = 1;
    ROI.xform = xform;
    ROI.voxelSize = voxelSize;
    ROI.date = datestr(now);
    % Add it to the view
    view = viewSet(view,'newROI',ROI);
  end
end

if exist('ROI','var')
  ROInum = viewGet(view,'ROInum',ROI.name);
  if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
    view = viewSet(view,'prevROIcoords',[]);
  end
end
