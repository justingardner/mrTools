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
baseSformCode = viewGet(v,'baseSformCode');
baseVol2mag = viewGet(v,'baseVol2mag');
baseVol2tal = viewGet(v,'baseVol2mag');
baseVoxelSize = viewGet(v,'baseVoxelSize');

% number of rois
numrois = viewGet(v,'numberofrois');
% put up a dialog with rois to delete
roinames = viewGet(v,'roiNames');
paramsDialog = {};
paramsDialog{end+1} = {'convertType',{'Convert','Adopt xform'},'Convert will convert the ROI to the current base anatomy xform permanently. This is not normally necessary (as ROIs are always converted on the fly to the base anatomy), but you may want to do this if you originally defined the ROI on a low resolution scan and want to have the ROI represented in a higher resolution. ''Adopt xform'' is only necessary in even more rare cases. This does not convert the roi voxels, but simply adopts the xform and voxel size of the base anatomy. If for example you changed the qform of your base anatomy you might need to use this.'};
for roinum = 1:length(roinames)
  % roi xform and voxel size
  roiXform = viewGet(v,'roiXform',roinum);
  roiVoxelSize = viewGet(v,'roiVoxelSize',roinum);
  % set help info
  if roinum == 1
    helpinfo = sprintf('Convert ROI %i: %s',roinum,roinames{roinum});
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
if length(paramsDialog)>1
  params = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert to [%0.2g %0.2g %0.2g] resolution',baseVoxelSize(1),baseVoxelSize(2),baseVoxelSize(3)));
else
  mrWarnDlg('(convertROI) No ROIs need conversion to current base xform');
  return
end

currentROIName = viewGet(v,'roiname');
needToRefresh = 0;
% now go through and do conversion
if ~isempty(params)
  % now go through and delete anything the user selected
  for roinum = 1:length(roinames)
    if isfield(params,fixBadChars(roinames{roinum})) && params.(fixBadChars(roinames{roinum}))
      % get the roi
      thisroinum = viewGet(v,'roinum',roinames{roinum});
      roi = viewGet(v,'ROI',thisroinum);
      roiBaseCoords = getROICoordinates(v,thisroinum,0);
      if ~isempty(roiBaseCoords)
	if strcmp(params.convertType,'Convert')
	  disp(sprintf('(convertROI) Converting ROI %i:%s',roinum,roinames{roinum}));
	  roi.coords = roiBaseCoords;
	else
	  disp(sprintf('(convertROI) Adopting base xform for ROI %i:%s',roinum,roinames{roinum}));
	end
	% set roi xform and voxel size
	roi.sformCode = baseSformCode;
	roi.xform = baseXform;
	% now set the vol2mag and vol2tal correctly
	roi.vol2mag = baseVol2mag;
	roi.vol2tal = baseVol2tal;
	roi.voxelSize = baseVoxelSize;
	% remove the roi
	v = viewSet(v,'deleteROI',thisroinum);
	% and add it back with the new coordinates
	v = viewSet(v,'newROI',roi);
	needToRefresh = 1;
      else
	mrWarnDlg(sprintf('(convertROI) ROI %i:%s has empty coordinates in transformation, skipping...',roinum,roinames{roinum}));
      end
    end
  end
  if needToRefresh
    % select the same current roi
    v = viewSet(v,'curROI',viewGet(v,'roinum',currentROIName));
    refreshMLRDisplay(viewGet(v,'viewNum'));
  end
end

