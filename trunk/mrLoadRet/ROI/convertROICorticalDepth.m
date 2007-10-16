% convertROICorticalDepth.m
%
%        $Id$
%      usage: convertROICorticalDepth()
%         by: justin gardner
%       date: 10/15/07
%    purpose: used to extend or restrict ROI coordinates across
%             cortical depths
%
function v = convertROICorticalDepth(v)

% check arguments
if ~any(nargin == [1])
  help convertROI
  return
end

% number of rois
numrois = viewGet(v,'numberofrois');
if numrois == 0
  mrWarnDlg('(convertROICorticalDepth) No currently loaded ROIs');
  return
end

% get cortical depth
corticalDepth = viewGet(v,'corticalDepth');

% put up some parameter choices
paramsInfo = {};
paramsInfo{end+1} = {'conversionType',{'Project through depth','Restrict to reference depth'},'If you set project through depth, then this will add all the voxels from each cortical depth that are in the same position as the ones at the reference depth. If you set to restrict to reference depth, this will remove any voxels that are not on the reference depth (note that you will still see some voxels on other depths, but those are voxels that exist at the reference depth--also, voxels that do not exist on this flat map will not be affected)'};
paramsInfo{end+1} = {'referenceDepth',corticalDepth,'min=0','max=1','incdec=[-0.1 0.1]','The cortical depth to start from'};
paramsInfo{end+1} = {'minDepth',0,'min=0','max=1','incdec=[-0.1 0.1]','The start depth'};
paramsInfo{end+1} = {'depthStep',0.1,'min=0','max=1','incdec=[-0.1 0.1]','The depth step (i.e. we will go from minDepth:depthStep:maxDepth (skipping the reference depth), including or excluding voxels'};
paramsInfo{end+1} = {'maxDepth',max(0.5,corticalDepth),'min=0','max=1','incdec=[-0.1 0.1]','The end depth'};

params = mrParamsDialog(paramsInfo,'ROI cortical depth conversion');
if isempty(params),return,end

% now select rois
% put up a dialog with rois to select
roinames = viewGet(v,'roiNames');
paramsDialog = {};
for roinum = 1:length(roinames)
  helpinfo = sprintf('Convert cortical depth of ROI %i: %s',roinum,roinames{roinum});
  paramsDialog{end+1} = {fixBadChars(roinames{roinum}),0,'type=checkbox',helpinfo};
end

% put up dialog
whichROI = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert cortical depth'));

% get base info
baseXform = viewGet(v,'baseXform');
baseVoxelSize = viewGet(v,'baseVoxelSize');
baseCoordMap = viewGet(v,'baseCoordMap',[],params.referenceDepth);
baseDims = baseCoordMap.dims;
flatDims = viewGet(v,'baseDims');
baseCoordMap = round(baseCoordMap.coords);
referenceBaseCoordMap = mrSub2ind(baseDims,baseCoordMap(:,:,:,1),baseCoordMap(:,:,:,2),baseCoordMap(:,:,:,3));
referenceBaseCoordMap = referenceBaseCoordMap(:);

currentROI = viewGet(v,'currentROI');
% now go through and do conversion
if ~isempty(whichROI)
  needToRefresh = 0;
  % now go through and delete anything the user selected
  for roinum = 1:length(roinames)
    if whichROI.(fixBadChars(roinames{roinum}))
      needToRefresh = 1;
      disppercent(-inf,sprintf('(convertROICorticalDepth) Processing ROI %i:%s',roinum,roinames{roinum}));
      % get the roi
      v = viewSet(v,'curROI',roinum);
      roi = viewGet(v,'ROI');
      % get the roiBaseCoords
      roiBaseCoords = getROIBaseCoords(v,roinum,baseXform,baseVoxelSize);
      roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
      % now find which baseCoords are in the current roi
      isInROI = ismember(referenceBaseCoordMap,roiBaseCoordsLinear);
      % make sure to keep the voxels at the reference depth
      roiBaseCoordsReferenceLinear = roiBaseCoordsLinear(ismember(roiBaseCoordsLinear,referenceBaseCoordMap));
      % now get each cortical depth, and add/remove voxels
      for corticalDepth = params.minDepth:params.depthStep:params.maxDepth
	% get the coordinates at this depth
	baseCoordMap = viewGet(v,'baseCoordMap',[],corticalDepth);
	baseCoordMap = round(baseCoordMap.coords);
	baseCoordMap = mrSub2ind(baseDims,baseCoordMap(:,:,:,1),baseCoordMap(:,:,:,2),baseCoordMap(:,:,:,3));
	baseCoordMap = baseCoordMap(:);
	% add the coordinates to our list
	roiBaseCoordsLinear = union(roiBaseCoordsLinear,baseCoordMap(isInROI));
      end
      roiBaseCoordsLinear = roiBaseCoordsLinear(~isnan(roiBaseCoordsLinear));
      % now convert back to regular coords
      roiBaseCoords = [];
      [roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsLinear);
      roiBaseCoords(4,:) = 1;
      % set the coordinates
      if strcmp(params.conversionType,'Project through depth')
	% add them to the ROI
	v = modifyROI(v,roiBaseCoords,baseXform,baseVoxelSize,1);
      else
	% get current coords
	curROICoords = viewGet(v,'roiCoords',roinum);
	% remove them from the ROI
	v = modifyROI(v,roiBaseCoords,baseXform,baseVoxelSize,0);
	% but make sure we have the voxels at the reference depth
	roiBaseCoords = [];
	[roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsReferenceLinear);
	roiBaseCoords(4,:) = 1;
	v = modifyROI(v,roiBaseCoords,baseXform,baseVoxelSize,1);
	% and save for undo (note we do this instead of allowing
	% modifyROI to do it since we have called modifyROI twice)
	v = viewSet(v,'prevROIcoords',curROICoords);
      end
      disppercent(inf);
    end
  end
  v = viewSet(v,'currentROI',currentROI);
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
