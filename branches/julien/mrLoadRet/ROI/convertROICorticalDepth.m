% convertROICorticalDepth.m
%
%        $Id$
%      usage: convertROICorticalDepth()
%         by: justin gardner
%       date: 10/15/07
%    purpose: used to extend or restrict ROI coordinates across
%             cortical depths
%
%  12/8/08 Modified by Taosheng Liu to take params. If params is set, GUI
%  will not show for setting params, also it assumes then all ROIs
%  associated with a view will be converted.

function [v params] = convertROICorticalDepth(v,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4])
  help convertROICorticalDepth
  return
end

eval(evalargs(varargin,[],[],{'justGetParams','defaultParams'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
% number of rois
numrois = viewGet(v,'numberofrois');
if numrois == 0
  mrWarnDlg('(convertROICorticalDepth) No currently loaded ROIs');
  return
end
roinames = viewGet(v,'roiNames');

if ieNotDefined('params')
  % get cortical depth
  corticalDepth = viewGet(v,'corticalDepth');
  paramsInfo = {};
  paramsInfo{end+1} = {'conversionType',{'Project through depth','Restrict to reference depth'},'type=popupmenu','If you set project through depth, then this will add all the voxels from each cortical depth that are in the same position as the ones at the reference depth. If you set to restrict to reference depth, this will remove any voxels that are not on the reference depth (note that you will still see some voxels on other depths, but those are voxels that exist at the reference depth--also, voxels that do not exist on this flat map will not be affected)'};
  paramsInfo{end+1} = {'referenceDepth',corticalDepth,'min=0','max=1','incdec=[-0.1 0.1]','The cortical depth to start from'};
  paramsInfo{end+1} = {'minDepth',0,'min=0','max=1','incdec=[-0.1 0.1]','The start depth'};
  paramsInfo{end+1} = {'depthStep',0.1,'min=0','max=1','incdec=[-0.1 0.1]','The depth step (i.e. we will go from minDepth:depthStep:maxDepth (skipping the reference depth), including or excluding voxels'};
  paramsInfo{end+1} = {'maxDepth',max(1,corticalDepth),'min=0','max=1','incdec=[-0.1 0.1]','The end depth'};
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    % put up some parameter choices
    params = mrParamsDialog(paramsInfo,'ROI cortical depth conversion');
  end
  % now select rois
  % put up a dialog with rois to select
  paramsDialog = {};
  for roinum = 1:length(roinames)
    helpinfo = sprintf('Convert cortical depth of ROI %i: %s',roinum,roinames{roinum});
    paramsDialog{end+1} = {fixBadChars(roinames{roinum}),0,'type=checkbox',helpinfo};
  end
  paramsDialog{end+1} = {'all',0,'type=checkbox','Select all ROIs'};
  % put up dialog
  whichROI = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert cortical depth'));
else
  disp('(convertROICorticalDepth) coverting all ROIs in the view');
  whichROI.all=1;
end
if isempty(params),return,end
% just return parameters
if justGetParams, return, end

% get base info
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
  % now go through and convert anything the user selected
  for roinum = 1:length(roinames)
    if whichROI.all || whichROI.(fixBadChars(roinames{roinum})) 
      needToRefresh = 1;
      disppercent(-inf,sprintf('(convertROICorticalDepth) Processing ROI %i:%s',roinum,roinames{roinum}));
      % get the roi
      v = viewSet(v,'curROI',roinum);
      roi = viewGet(v,'ROI');
      % get the roiBaseCoords
      roiBaseCoords = getROICoordinates(v,roinum,0);
      if isempty(roiBaseCoords)
        disppercent(inf);
        disp(sprintf('(convertROICorticalDepth) %s has no coordinates on this flat',roinames{roinum}));
        continue;
      end
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
      base2roi = viewGet(v,'base2roi');
      if strcmp(params.conversionType,'Project through depth')
        % add them to the ROI
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,1);
      else
        % get current coords
        curROICoords = viewGet(v,'roiCoords',roinum);
        % remove them from the ROI
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
        % but make sure we have the voxels at the reference depth
        roiBaseCoords = [];
        [roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsReferenceLinear);
        roiBaseCoords(4,:) = 1;
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,1);
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

