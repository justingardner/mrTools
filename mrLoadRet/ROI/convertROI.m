% convertROI.m
%
%        $Id$
%      usage: convertROI()
%         by: justin gardner
%       date: 10/15/07
%    purpose: convert ROIs to current base anatomy xform/voxelSize
%
function v = convertROI(v)

% check arguments
if ~any(nargin == [1])
  help convertROI
  return
end

% base xform and voxel size
baseXform = viewGet(v,'baseXform');
baseVoxelSize = viewGet(v,'baseVoxelSize');

% number of rois
numrois = viewGet(v,'numberofrois');
% put up a dialog with rois to delete
roinames = viewGet(v,'roiNames');
paramsDialog = {};
for roinum = 1:length(roinames)
  % roi xform and voxel size
  roiXform = viewGet(v,'roiXform',roinum);
  roiVoxelSize = viewGet(v,'roiVoxelSize',roinum);
  % set help info
  if roinum == 1
    helpinfo = sprintf('Convert ROI %i: %s. Selcting this will convert the ROI to the current base anatomy xform permanently. This is not normally necessary, but you may want to do this if you originally defined the ROI on a low resolution scan and want to have the ROI represented in a higher resolution.',roinum,roinames{roinum});
  else
    helpinfo = sprintf('Convert ROI %i: %s',roinum,roinames{roinum});
  end
  % only add the roi if it does not mach current base anatomy
  if ~isequal(roiXform,baseXform) || ~isequal(baseVoxelSize,roiVoxelSize)
    paramsDialog{end+1} = {fixBadChars(roinames{roinum}),0,'type=checkbox',helpinfo};
    paramsDialog{end+1} = {sprintf('%s_voxelSize',fixBadChars(roinames{roinum})),roiVoxelSize,'editable=0',sprintf('Current voxel size for roi %s',roinames{roinum})};
  end
end

% put up dialog
if ~isempty(paramsDialog)
  params = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert to [%s] resolution',num2str(baseVoxelSize)));
else
  mrWarnDlg('(convertROI) No ROIs need conversion to current base xform');
  return
end

needToRefresh = 0;
% now go through and do conversion
if ~isempty(params)
  % now go through and delete anything the user selected
  for roinum = 1:length(roinames)
    if isfield(params,fixBadChars(roinames{roinum})) && params.(fixBadChars(roinames{roinum}))
      % get the roi
      roi = viewGet(v,'ROI',roinum);
      roiBaseCoords = getROIBaseCoords(v,roinum,baseXform,baseVoxelSize);
      if ~isempty(roiBaseCoords)
	disp(sprintf('(convertROI) Converting ROI %i:%s',roinum,roinames{roinum}));
	% valid roi coordinates, change the roi
	roi.coords = roiBaseCoords;
	roi.xform = baseXform;
	roi.voxelSize = baseVoxelSize;
	% remove the roi
	v = viewSet(v,'deleteROI',roinum);
	% and add it back with the new coordinates
	v = viewSet(v,'newROI',roi);
	needToRefresh = 1;
      else
	mrWarnDlg(sprintf('(convertROI) ROI %i:%s has empty coordinates in transformation, skipping...',roinum,roinames{roinum}));
      end
      %      v = viewSet(v,'deleteROI',viewGet(v,'roinum',roinames{roinum}));
    end
  end
  if needToRefresh
    refreshMLRDisplay(viewGet(v,'viewNum'));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getROIBaseCoords   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiBaseCoords = getROIBaseCoords(v,roinum,baseXform,baseVoxelSize)

% look in cache for converted coordinates
roiCache = viewGet(v,'ROICache',roinum);
if isempty(roiCache)
  disppercent(-inf,sprintf('Computing ROI base coordinates for %i:%s',roinum,viewGet(v,'roiName',roinum)));
  %viewGet
  roiCoords = viewGet(v,'roiCoords',roinum);
  roiXform = viewGet(v,'roiXform',roinum);
  roiVoxelSize = viewGet(v,'roiVoxelSize',roinum);
  if ~isempty(roiCoords) & ~isempty(roiXform) & ~isempty(baseXform)
    % Use xformROI to supersample the coordinates
    roiBaseCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));
  else
    roiBaseCoords = [];
  end
else
  roiBaseCoords = roiCache.roiBaseCoords;
end
