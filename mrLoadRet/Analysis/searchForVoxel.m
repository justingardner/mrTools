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

% ask voxel number
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
params = mrParamsDialog(paramsInfo,'Look for voxel');
if isempty(params)
  % refresh and retrun
  refreshMLRDisplay(viewNum);
  return
end


% switch to the chosen base
if ~strcmp(baseName,params.baseName)
  v = viewSet(v,'curBase',viewGet(v,'baseNum',params.baseName));
end

% now switch to the appropriate orientation and slice
if viewGet(v,'baseType') == 0
  % check orientation
  whichSliceOrientation = find(strcmp(params.baseOrientation,sliceOrientations));
  % if user didn't ask for the default, set orientation apporpriately
  if whichSliceOrientation ~= 1
    v = viewSet(v,'sliceOrientation',whichSliceOrientation-1);
  end
  % first get which dimension we are looking at
  baseSliceIndex = viewGet(v,'baseSliceIndex');
  % now find voxel in base coordinates
  baseXform = viewGet(v,'baseXform');
  scanXform = viewGet(v,'scanXform',viewGet(v,'curScan'));
  shiftXform = shiftOriginXform;
  baseVoxel = round(inv(shiftXform)*inv(baseXform)*scanXform*shiftXform*[params.x params.y params.s 1]');
  % now get the slice to switch to
  sliceNum = baseVoxel(baseSliceIndex);
  % make sure the asked for slice exists
  baseDims = viewGet(v,'baseDims');
  if ((sliceNum > 0) && (sliceNum <= baseDims(baseSliceIndex)))
    % if so, switch to it
    mlrGuiSet(viewNum,'slice',sliceNum);
  end
end
% refresh the display
refreshMLRDisplay(viewNum);
hold on

% get the refreshed view
v = viewGet([],'view',viewNum);
  
% get coords of chosen voxel
overlayCoords = round(viewGet(v,'cursliceOverlayCoords'));
match = (overlayCoords(:,:,1) == params.x) & (overlayCoords(:,:,2) == params.y) & (overlayCoords(:,:,3) == params.s);

% no match, give up
if isempty(find(match))
  disp(sprintf('(searchForVoxel) No matching voxel'));
  return
end

% and then display the voxel
[y x] = find(match);
for i = 1:length(x)
  plot(x,y,'.','Color',color2RGB(params.color));
end

