% convertROICorticalDepth.m
%
%        $Id$
%      usage: convertROICorticalDepth(v, params, <'justGetParams=0/1','defaultParams=0/1', 'roiList', roiList>)
%         by: justin gardner
%       date: 10/15/07
%    purpose: used to extend or restrict ROI coordinates across
%             cortical depths
%
%  12/8/08 Modified by Taosheng Liu to take params. If params is set, GUI
%  will not show for setting params, also it assumes then all ROIs
%  associated with a view will be converted.
%
%             To just get a default parameter structure:
% 
%             v = newView;
%             [v params] = convertROICorticalDepth(v,[],'justGetParams=1');
%             [v params] = convertROICorticalDepth(v,[],'justGetParams=1','defaultParams=1');
%             [v params] = convertROICorticalDepth(v,[],'justGetParams=1','defaultParams=1','roiList=[1 2]');


function [v params] = convertROICorticalDepth(v,params,varargin)

% check arguments
if nargin < 1
  help convertROICorticalDepth
  return
end

% optional arguments
getArgs(varargin);

if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('distanceThreshold'), distanceThreshold = 2; end %distance threshold (in mm) to exclude

% number of rois
numrois = viewGet(v,'numberofrois');
if numrois == 0
  mrWarnDlg('(convertROICorticalDepth) No currently loaded ROIs');
  return
end

if ieNotDefined('params')
  askForParams = 1;
  % get cortical depth
  corticalDepth = viewGet(v,'corticalDepth');
  referenceDepth= mean(corticalDepth);
  if length(corticalDepth)==2 && corticalDepth(1)~=corticalDepth(2)
    minDepth = corticalDepth(1);
    maxDepth = corticalDepth(2);
  else
    minDepth = 0;
    maxDepth = 1;
  end
  depthStep = 1/(mrGetPref('corticalDepthBins')-1);
  incdecString = sprintf('incdec=[-%f %f]',depthStep,depthStep);
  while askForParams
    paramsInfo = {};
    paramsInfo{end+1} = {'conversionType',{'Project through depth','Restrict to reference depth'},'type=popupmenu','If you set project through depth, then this will add all the voxels from each cortical depth that are in the same position as the ones at the reference depth. If you set to restrict to reference depth, this will remove any voxels that are not on the reference depth (note that you will still see some voxels on other depths, but those are voxels that exist at the reference depth--also, voxels that do not exist on this flat map will not be affected)'};
    paramsInfo{end+1} = {'referenceDepth',referenceDepth,'min=0','max=1',incdecString,'The cortical depth to start from'};
    paramsInfo{end+1} = {'minDepth',minDepth,'max=1',incdecString,'The minimum depth. Negative values will extend the ROI into white matter'};
    paramsInfo{end+1} = {'depthStep',depthStep,'min=0','max=1',incdecString,'The depth step (i.e. we will go from minDepth:depthStep:maxDepth (skipping the reference depth), including or excluding voxels'};
    paramsInfo{end+1} = {'maxDepth',maxDepth,'min=0','max=1',incdecString,'The maximum depth'};
    paramsInfo{end+1} = {'excludeOtherVoxels',1,'type=checkbox','If ROI voxels exist oustide the projected surface, they will be removed. Uncheck to keep them. This option is ignored if restriction is selected'};
    paramsInfo{end+1} = {'allowProjectionThroughSulci',1,'type=checkbox','Voxels will be kept even if they also belong to another part of the cortical surface through a sulcus. Uncheck to exclude these voxels. Note that voxels projected to another part of the cortex through white matter (for instance in case minDepth is negative) will be excluded too. If the ROI is projected based on a flat map, only voxels on the flat map, not the whole surface, will be excluded. This option is ignored if restriction is selected'};
    if defaultParams
      params = mrParamsDefault(paramsInfo);
    else
      % put up some parameter choices
      params = mrParamsDialog(paramsInfo,'ROI cortical depth conversion');
    end
    % Abort if params empty
    if ieNotDefined('params'),return,end
    
    if 0
      %checks on params here if needed
    else
      askForParams = 0;
      % now select rois
      % put up a dialog with rois to select
      if defaultParams
        params.roiList = viewGet(v,'curROI');
      else
        params.roiList = selectInList(v,'rois');
        if isempty(params.roiList)
          askForParams = 1;
        end
      end
    end
  end
end
if isempty(params),return,end

if ~ieNotDefined('roiList')
  params.roiList = roiList;
end

% just return parameters
if justGetParams, return, end

%remember what ROIs were selected in the view for later
currentROI = viewGet(v,'currentROI');
% now go through and convert anything the user selected
for roinum = params.roiList
  roiName = viewGet(v,'roiname', roinum);
  mlrDispPercent(-inf,sprintf('(convertROICorticalDepth) Processing ROI %i:%s',roinum,roiName));
  % get the roi
  v = viewSet(v,'curROI',roinum);
  % now try to figure out what base this was created on
  roiCreatedOnBase = viewGet(v,'roiCreatedOnBase',roiName);
  if isempty(roiCreatedOnBase)
    fprintf('(convertROICorticalDepth) Converting %s based on base:%s because roiCreatedOnBase has not been set.\n',roiName,viewGet(v,'baseName'));
    baseNum = viewGet(v,'curBase');
  else
    % get the basenumber for the base that this was created on
    baseNum = viewGet(v,'baseNum',roiCreatedOnBase);
    if isempty(baseNum)
      fprintf('(convertROICorticalDepth) Converting %s based on base:%s because base:%s which this roi was created on is not loaded\n',roiName,viewGet(v,'baseName'),roiCreatedOnBase);
      baseNum = viewGet(v,'curBase');
    end
    if viewGet(v,'basetype',baseNum)==0
      fprintf('(convertROICorticalDepth) Converting %s based on base:%s because base:%s which this roi was created on is not a surface or a flat map\n',roiName,viewGet(v,'baseName'),roiCreatedOnBase);
      baseNum = viewGet(v,'curBase');
    end
  end
  % check the base type to see if it's compatible with the current implementation of params.allowProjectionThroughSulci
  if ~params.allowProjectionThroughSulci
    if viewGet(v,'basetype',baseNum) == 2
      mrWarnDlg('(convertROICorticalDepth) Unchecking allowProjectionThroughSulci parameters is not yet implemented for surface bases')
      return
    end
  end
  % get the roi transformation in order to set the coordinates later
  base2roi = viewGet(v,'base2roi',roinum,baseNum);
  % get the roiBaseCoords
  roiBaseCoords = getROICoordinates(v,roinum,[],[],'baseNum',baseNum);
  if isempty(roiBaseCoords)
    mlrDispPercent(inf);
    mrWarnDlg(sprintf('(convertROICorticalDepth) %s has no coordinates on this flat',roiName));
    continue;
  end
  nVoxelsOriginalROI = size(roiBaseCoords,2);
  % get base info, including (rounded) base coordinates corresponding to the reference cortical depth
  baseVoxelSize = viewGet(v,'baseVoxelSize',baseNum);
  baseCoordMap = viewGet(v,'baseCoordMap',baseNum,params.referenceDepth);
  mapDims = size(baseCoordMap.coords);
  baseDims = baseCoordMap.dims;
  baseCoordMap = round(baseCoordMap.coords);
  referenceBaseCoordMap = mrSub2ind(baseDims,baseCoordMap(:,:,:,1),baseCoordMap(:,:,:,2),baseCoordMap(:,:,:,3));
  referenceBaseCoordMap = referenceBaseCoordMap(:);
  % get roi linear coordinates
  roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
  % now find which baseCoords are in the current roi at the reference depth
  [isInROI, roiInBase] = ismember(referenceBaseCoordMap,roiBaseCoordsLinear);
  % get the roi base coordinates that are found in base at the reference depth
  roiInBase = unique(setdiff(roiInBase,0));
  % (note that here we could have used ismember(roiBaseCoordsLinear,referenceBaseCoordMap) instead, which is perhaps easier to understand)
  
  % if we don't find most of the coordinates, then
  % probably good to complain and give up
  if (length(roiInBase)/length(roiBaseCoordsLinear)) < 0.1
    mlrDispPercent(inf);
    mrWarnDlg(sprintf('(convertROICorticalDepth) !!! %s has less than %0.0f%% coordinates on surface %s. Perhaps you need to load the base that it was orignally created on. !!!',roiName,ceil(100*(length(roiInBase)/length(roiBaseCoordsLinear))),viewGet(v,'baseName',baseNum)));
    continue;
  end
  % make sure to keep the voxels at the reference depth
  roiBaseCoordsReferenceLinear = roiBaseCoordsLinear(ismember(roiBaseCoordsLinear,referenceBaseCoordMap));
  % (Note that we could have used roiBaseCoordsLinear(roiInBase) instead here)

  % if excluding voxels that belong to two distant locations of the cortex, we need to compute the shortest distance of all elements
  % of the flat map or surface to the ROI, in flat/surface space. (This is not yet implemented for surfaces and would require computing the Dijkstra distance)
  if ~params.allowProjectionThroughSulci
    [flatCoordsX, flatCoordsY] = meshgrid(1:mapDims(1),1:mapDims(2)); %compute the coordinates of the points on the flat map (i.e. in flat space). 
    % (Distance calculations could be done using the actual surface locations, but this would require computing Dijkstra distances.
    % Easier like this and sufficient for our purposes, until the same is implemented for surfaces)
    % find coordinates of the ROI on the flat map
    flatCoordsRoiX = flatCoordsX(isInROI);
    flatCoordsRoiY = flatCoordsY(isInROI);
    % for each pixel of the flat map that corresponds to a base voxels, but that does not belong to the ROI,
    % compute its shortest distance to any pixel in the ROI (in flat space)
    minDistanceToROI = zeros(mapDims(1)*mapDims(2),1);
    for iPixel = find(~isInROI & ~isnan(referenceBaseCoordMap))'
      minDistanceToROI(iPixel) = min(sqrt((flatCoordsX(iPixel) - flatCoordsRoiX).^2 + (flatCoordsY(iPixel) - flatCoordsRoiY).^2));
    end
    % now find flat map pixels that are a minimum distance from any pixel in the ROI.
    % first compute approximate pixel size (based on surface coordinates only at the reference depth)
    % separate x,y and z coordinates of flat map in base space and convert to mm
    xBaseCoordMap = baseVoxelSize(1)*baseCoordMap(:,:,1,1);
    yBaseCoordMap = baseVoxelSize(2)*baseCoordMap(:,:,1,2);
    zBaseCoordMap = baseVoxelSize(3)*baseCoordMap(:,:,1,3);
    % remove pixels that do not index a location on the surface
    xBaseCoordMap(isnan(referenceBaseCoordMap))=NaN;
    yBaseCoordMap(isnan(referenceBaseCoordMap))=NaN;
    zBaseCoordMap(isnan(referenceBaseCoordMap))=NaN;
    %compute pixel size in pixels
    pixelSize(1) = nanmean(nanmean(sqrt(diff(xBaseCoordMap,1,1).^2 + diff(yBaseCoordMap,1,1).^2 + diff(zBaseCoordMap,1,1).^2)));
    pixelSize(2) = nanmean(nanmean(sqrt(diff(xBaseCoordMap,1,2).^2 + diff(yBaseCoordMap,1,2).^2 + diff(zBaseCoordMap,1,2).^2)));
    %compute distance threshold in pixels 2 based on approximate pixel size
    isFarEnoughFromROI = minDistanceToROI > (distanceThreshold / min(pixelSize));
  end
  
  if strcmp(params.conversionType,'Project through depth')
    roiBaseCoordsLinear=[];
    % now get each cortical depth, and add/remove voxels
    corticalDepths = params.minDepth:params.depthStep:params.maxDepth;
    % (negative cortical depths mean that the flat/surface base (and subsequently the ROI) will be extended into white matter
    baseCoordMap = viewGet(v,'baseCoordMap',baseNum,corticalDepths);
    for iDepth = 1:size(baseCoordMap.coords,5)
      % get the (rounded) base coordinates at this depth
      baseCoords = round(baseCoordMap.coords(:,:,:,:,iDepth));
      baseCoords = mrSub2ind(baseDims,baseCoords(:,:,:,1),baseCoords(:,:,:,2),baseCoords(:,:,:,3));
      baseCoords = baseCoords(:);
      % find the baseCoords that are in the ROI at this depth and add them to our list
      roiBaseCoordsLinear = union(roiBaseCoordsLinear,baseCoords(isInROI));
    end
    roiBaseCoordsLinear = roiBaseCoordsLinear(~isnan(roiBaseCoordsLinear));
    % exclude any projected voxel that ends up in a different part of cortex either through a sulcus or through white matter
    if ~params.allowProjectionThroughSulci
      % get the (rounded) base coordinates at all depths between 0 and 1
      baseCoords = round(baseCoordMap.coords(:,:,:,:,corticalDepths>=0 & corticalDepths<=1));
      baseCoords = mrSub2ind(baseDims,baseCoords(:,:,:,1,:),baseCoords(:,:,:,2,:),baseCoords(:,:,:,3,:));
      baseCoords = reshape(baseCoords,size(baseCoords,1)*size(baseCoords,2)*size(baseCoords,3),size(baseCoords,5));
      % exclude any voxel that is also part of the flat map (or surface) at a location away form the ROI
      roiBaseCoordsLinear = setdiff(roiBaseCoordsLinear,baseCoords(isFarEnoughFromROI,:));
    end
    % now convert back to regular coords
    additionalRoiBaseCoords = [];
    [additionalRoiBaseCoords(1,:), additionalRoiBaseCoords(2,:), additionalRoiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsLinear);
    additionalRoiBaseCoords(4,:) = 1;
    %clear all existing voxels if we're not keeping voxels outside the projection
    if params.excludeOtherVoxels
      % remove everything from the ROI
      roiBaseCoords(4,:) = 1;
      v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
    end
    % add the projected voxels to the ROI
    v = modifyROI(v,additionalRoiBaseCoords,base2roi,baseVoxelSize,1);
  else
    % get current coords
    curROICoords = viewGet(v,'roiCoords',roinum);
    % remove them from the ROI
    v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
    % but make sure we have the voxels at the reference depth
    additionalRoiBaseCoords = [];
    [additionalRoiBaseCoords(1,:), additionalRoiBaseCoords(2,:), additionalRoiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsReferenceLinear);
    additionalRoiBaseCoords(4,:) = 1;
    v = modifyROI(v,additionalRoiBaseCoords,base2roi,baseVoxelSize,1);
    % and save for undo (note we do this instead of allowing
    % modifyROI to do it since we have called modifyROI twice)
    v = viewSet(v,'prevROIcoords',curROICoords);
  end
  mlrDispPercent(inf);
  fprintf(1,'(convertROICorticalDepth) Number of voxels in original ROI: %d\t Number of voxels in modified ROI: %d\n',nVoxelsOriginalROI,size(additionalRoiBaseCoords,2));
end
v = viewSet(v,'currentROI',currentROI);

