function view = newROI(view,name,select,color,xform,voxelSize)

% function view = newROI(view,[name],[select],[color],[xform],[voxelSize])
%
% Makes new empty ROI, adds it to view.ROIs, and selects it.
%
% name: name (string) for the ROI.
%   Default: 'ROI<number>' where <number> is chosen (to avoid conflicts with
%   other ROI names) as one plus the number of previously existing ROIs.
% select: if non-zero, chooses the new ROI as the selectedROI
%    Default: 1.
% color: sets color for drawing the ROI.
%    Default: 'b'
% xform: 4x4 transform from functionals to ROI. 
%    Default:  uses the xform from the current base anatomy
% voxelSize: 3-vector specifying the voxel size of the ROI
%    Default: uses the voxel size of the current base anatomy
%
% djh, 7/2005 (modified from mrLoadRet-3.1)

if isempty(viewGet(view,'curBase')) & ieNotDefined('xform') & ieNotDefined('voxelSize')
  mrErrorDlg('You must load a base anatomy before creating an ROI.');
end

if ieNotDefined('name')
  name=sprintf('ROI%.0f',(length(view.ROIs)+1));
end
if ieNotDefined('select')
  select = 1;
end
if ieNotDefined('color')
  color = 'blue';
end
if ieNotDefined('xform')
  baseNum = viewGet(view,'currentBase');
  xform = viewGet(view,'basexform',baseNum);
end
if ieNotDefined('voxelSize')
  baseNum = viewGet(view,'currentBase');
  voxelSize = viewGet(view,'baseVoxelSize',baseNum);
end

% Set required fields. Additional (optional) optional fields are set by
% isroi which is called by viewSet newROI.
ROI.name = name;
ROI.viewType = view.viewType;
ROI.color = color;
ROI.xform = xform;
ROI.voxelSize = voxelSize;

% Add it to the view
view = viewSet(view,'newROI',ROI);

% Select it and reset view.prevCoords
if select
  ROInum = viewGet(view,'ROInum',name);
  if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
    view = viewSet(view,'prevROIcoords',[]);
  end
end

return;
