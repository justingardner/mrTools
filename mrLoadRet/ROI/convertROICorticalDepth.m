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

% number of rois
numrois = viewGet(v,'numberofrois');
if numrois == 0
  mrWarnDlg('(convertROICorticalDepth) No currently loaded ROIs');
  return
end

if ieNotDefined('params')
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
  paramsInfo = {};
  paramsInfo{end+1} = {'conversionType',{'Project through depth','Restrict to reference depth'},'type=popupmenu','If you set project through depth, then this will add all the voxels from each cortical depth that are in the same position as the ones at the reference depth. If you set to restrict to reference depth, this will remove any voxels that are not on the reference depth (note that you will still see some voxels on other depths, but those are voxels that exist at the reference depth--also, voxels that do not exist on this flat map will not be affected)'};
  paramsInfo{end+1} = {'referenceDepth',referenceDepth,'min=0','max=1',incdecString,'The cortical depth to start from'};
  paramsInfo{end+1} = {'minDepth',minDepth,'min=0','max=1',incdecString,'The start depth'};
  paramsInfo{end+1} = {'depthStep',depthStep,'min=0','max=1',incdecString,'The depth step (i.e. we will go from minDepth:depthStep:maxDepth (skipping the reference depth), including or excluding voxels'};
  paramsInfo{end+1} = {'maxDepth',maxDepth,'min=0','max=1',incdecString,'The end depth'};
  paramsInfo{end+1} = {'excludeOtherVoxels',1,'type=checkbox','If ROI voxels exist oustide the projected surface, they will be remove. Uncheck to keep them. this option is ignored if restriction is selected'};
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    % put up some parameter choices
    params = mrParamsDialog(paramsInfo,'ROI cortical depth conversion');
  end
  % now select rois
  % put up a dialog with rois to select
  if defaultParams
    params.roiList = viewGet(v,'curROI');
  else
    params.roiList = selectInList(v,'rois');
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
  disppercent(-inf,sprintf('(convertROICorticalDepth) Processing ROI %i:%s',roinum,roiName));
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
  % get the roi transformation in order to set the coordinates later
  base2roi = viewGet(v,'base2roi',roinum,baseNum);
  % get the roiBaseCoords
  roiBaseCoords = getROICoordinates(v,roinum,[],[],'baseNum',baseNum);
  if isempty(roiBaseCoords)
    disppercent(inf);
    mrWarnDlg(sprintf('(convertROICorticalDepth) %s has no coordinates on this flat',roiName));
    continue;
  end
  nVoxelsOriginalROI = size(roiBaseCoords,2);
  % get base info
  baseVoxelSize = viewGet(v,'baseVoxelSize',baseNum);
  baseCoordMap = viewGet(v,'baseCoordMap',baseNum,params.referenceDepth);
  baseDims = baseCoordMap.dims;
  baseCoordMap = round(baseCoordMap.coords);
  referenceBaseCoordMap = mrSub2ind(baseDims,baseCoordMap(:,:,:,1),baseCoordMap(:,:,:,2),baseCoordMap(:,:,:,3));
  referenceBaseCoordMap = referenceBaseCoordMap(:);
  % get roi linear coordinates
  roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
  % now find which baseCoords are in the current roi
  [isInROI roiInBase] = ismember(referenceBaseCoordMap,roiBaseCoordsLinear);
  % get the roi base coordinates that are found in base
  roiInBase = unique(setdiff(roiInBase,0));
  % if we don't find most of the coordinates, then
  % probably good to complain and give up
  if (length(roiInBase)/length(roiBaseCoordsLinear)) < 0.1
    disppercent(inf);
    mrWarnDlg(sprintf('(convertROICorticalDepth) !!! %s has less than %0.0f%% coordinates on surface %s. Perhaps you need to load the base that it was orignally created on. !!!',roiName,ceil(100*(length(roiInBase)/length(roiBaseCoordsLinear))),viewGet(v,'baseName',baseNum)));
    continue;
  end
  % make sure to keep the voxels at the reference depth
  roiBaseCoordsReferenceLinear = roiBaseCoordsLinear(ismember(roiBaseCoordsLinear,referenceBaseCoordMap));

  if strcmp(params.conversionType,'Project through depth')
    %clear all voxels if we're not keeping voxels outside the projection
    if params.excludeOtherVoxels
      % remove everything from the ROI
      roiBaseCoords(4,:) = 1;
      v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
    end
    roiBaseCoordsLinear=[];
    % now get each cortical depth, and add/remove voxels
    corticalDepths = params.minDepth:params.depthStep:params.maxDepth;
    baseCoordMap = viewGet(v,'baseCoordMap',baseNum,corticalDepths);
    for iDepth = 1:size(baseCoordMap.coords,5)
      % get the coordinates at this depth
      baseCoords = round(baseCoordMap.coords(:,:,:,:,iDepth));
      baseCoords = mrSub2ind(baseDims,baseCoords(:,:,:,1),baseCoords(:,:,:,2),baseCoords(:,:,:,3));
      baseCoords = baseCoords(:);
      % add the coordinates to our list
      roiBaseCoordsLinear = union(roiBaseCoordsLinear,baseCoords(isInROI));
    end
    roiBaseCoordsLinear = roiBaseCoordsLinear(~isnan(roiBaseCoordsLinear));
    % now convert back to regular coords
    roiBaseCoords = [];
    [roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsLinear);
    roiBaseCoords(4,:) = 1;
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
  fprintf(1,'(convertROICorticalDepth) Number of voxels in original ROI: %d\t Number of voxels in modified ROI: %d\n',nVoxelsOriginalROI,size(roiBaseCoords,2));
end
v = viewSet(v,'currentROI',currentROI);

