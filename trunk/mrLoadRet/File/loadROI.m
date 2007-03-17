function view = loadROI(view,filename,select)
%
% view = loadROI(view,[filename],[select])
%
% Loads an ROI and adds it to view.ROIs.
%
% If filename is not specified, prompts user to select a file. If filename
% is specified, it loads from view.subdir/ROIs/name.mat. Filename can
% be a string specifying an ROI file or it can be a cell array of
% filenames to load multiple ROIs at once.
%
% select: If non-zero, selects the loaded ROI as the current (selected)
% ROI. Default: 1.
%
% The file must contain a structure or structures, each with the following
% fields:
% - name: string
% - viewType: 'Volume', 'Surface', or 'Flat'
% - color: 
% - coords: 4xN array of ROI coordinates
% - xform: 4x4 homogeneous transform matrix that maps roi coordinates to base coordinate frame. 
% - voxelSize: 3-vector specifying the ROI voxel size.
%
% In addition the 'name' field is set to the variable name to ensure
% consistency.
%
% djh, 1/9/98
% 8/2005, djh, updated to mrLoadRet-4.0

mrGlobals

if ieNotDefined('select')
  select = 1;
end

% Path to overlays
startPathStr = viewGet(view,'roiDir');

% Complete pathStr
if ieNotDefined('filename')
  pathStr = getPathStrDialog(startPathStr,'Choose one or more ROIs','*.mat','on');
else
  if iscell(filename)
    pathStr = cell(size(filename));
    for p=1:length(pathStr)
      pathStr{p} = fullfile(startPathStr,[filename{p},'.mat']);
    end
  else
    pathStr = {fullfile(startPathStr,[filename,'.mat'])};
  end
end
if isempty(pathStr),return,end
if ~iscell(pathStr)
    pathStr = {pathStr};
end

% Load the file. Loop through the variables that were loaded and add
% each of them as a new ROI, setting roi.fieldnames as we go.
for p = 1:length(pathStr)
  if exist(pathStr{p},'file')
    s = load(pathStr{p});
    varNames = fieldnames(s);
    roi = eval(['s.',varNames{1}]);
    roi.name = varNames{1};
    % Add it to the view
    view = viewSet(view,'newROI',roi);
    % Select it and reset view.prevCoords
    if select
      ROInum = viewGet(view,'numberofROIs');
      if (ROInum > 0)
	view = viewSet(view,'currentROI',ROInum);
      end
    end
  else
    mrWarnDlg(['ROI ',pathStr{p},' not found.']);
  end
end

return;

% Test/debug
view = loadROI(MLR.views{1},'ROI1');
view = loadROI(MLR.views{1},{'ROI1','ROI2'});
view = loadROI(MLR.views{1});

