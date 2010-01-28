% searchForVoxel.m
%
%      usage: searchForVoxel()
%         by: justin gardner
%       date: 11/08/07
%    purpose: interrogator function that displays desired voxel in
%             a different view
%
function retval = searchForVoxel(v,overlayNum,scan,x,y,s,roi,varargin)

% get the view number
viewNum = viewGet(v,'viewNum');

% get parameters
paramsInfo = {};
scanDims = viewGet(v,'scanDims');
baseName = viewGet(v,'baseName');
sliceOrientations = {'default','axial','coronal','sagittal'};
baseNames = putOnTopOfList(baseName,viewGet(v,'baseNames'));
paramsInfo{end+1} = {'baseName',baseNames,'Choose which base to find voxel on'};
paramsInfo{end+1} = {'baseOrientation',sliceOrientations,'Choose which slice orientation to view voxel in'};
paramsInfo{end+1} = {'x',x,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(1)),'Choose X coordinate to look for'};
paramsInfo{end+1} = {'y',y,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(2)),'Choose Y coordinate to look for'};
paramsInfo{end+1} = {'s',s,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(3)),'Choose S coordinate to look for'};
paramsInfo{end+1} = {'color',color2RGB,'Choose color to display voxel in'};
paramsInfo{end+1} = {'maxSearchRadius',5,'type=numeric','If there is no exact match, then will display the closest voxel within this search radius. Set smaller to force only display of more exact matches. Set higher if you want to allow matches that are farther away.'};
params = mrParamsDialog(paramsInfo,'Look for voxel');
if isempty(params)
  % if user hit cancel, then refresh
  % display to remove any marked voxel and return
  refreshMLRDisplay(viewNum);
  return
end

% switch to the chosen base
if ~strcmp(baseName,params.baseName)
  v = viewSet(v,'curBase',viewGet(v,'baseNum',params.baseName));
end

% now switch to the appropriate orientation and slice
baseType = viewGet(v,'baseType');
if baseType == 0
  % check orientation
  whichSliceOrientation = find(strcmp(params.baseOrientation,sliceOrientations));
  % if user didn't ask for the default, set orientation apporpriately
  if whichSliceOrientation ~= 1
    v = viewSet(v,'sliceOrientation',whichSliceOrientation-1);
  end
  % first get which dimension we are looking at
  baseSliceIndex = viewGet(v,'baseSliceIndex');
  % now find voxel in base coordinates
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');
  % now get the slice to switch to
  sliceNum = baseVoxel(baseSliceIndex);
  % make sure the asked for slice exists
  baseDims = viewGet(v,'baseDims');
  if ((sliceNum > 0) && (sliceNum <= baseDims(baseSliceIndex)))
    % if so, switch to it
    viewSet(v,'curSlice',sliceNum);
  end
end

% refresh the display
refreshMLRDisplay(viewNum);
hold on

% get the refreshed view
v = viewGet([],'view',viewNum);
  
% and then display the voxel
if baseType <= 1
  % find voxel in base coordinates
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');

  % find the coordinates of the view
  baseCoords = round(viewGet(v,'cursliceBaseCoords'));

  % compute distance to matching voxel
  matchDistance = sqrt((baseCoords(:,:,1) - baseVoxel(1)).^2 + (baseCoords(:,:,2) - baseVoxel(2)).^2 + (baseCoords(:,:,3) - baseVoxel(3)).^2);
  minMatchDistance = min(matchDistance(:));
  
  % if there is no match within 5 voxels, give up
  if minMatchDistance > params.maxSearchRadius
    mrWarnDlg(sprintf('(searchForVoxel) No voxel within %f of searched for voxel. Closest voxel is within a radius of %f',params.maxSearchRadius,minMatchDistance));
    hold off
    return
  end

  % display info
  disp(sprintf('(searchForVoxel) Scan voxel [%i %i %i] projects to [%i %i %i] in %s. Displaying %i voxels within a radius of %f.',params.x,params.y,params.s,baseVoxel(1),baseVoxel(2),baseVoxel(3),params.baseName,length(find(matchDistance==minMatchDistance)),minMatchDistance));

  [y x] = find(matchDistance==minMatchDistance);
  for i = 1:length(x)
    plot(x,y,'.','Color',color2RGB(params.color));
  end
else
  % if we need to display on a surface,
  % first get base coord of point
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');
  
  % and plot
  plot3(baseVoxel(1),baseVoxel(2),baseVoxel(3),'.','Color',color2RGB(params.color),'MarkerSize',30);
end

hold off
