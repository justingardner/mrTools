% baseInfo.m
%
%      usage: baseInfo(view)
%         by: justin gardner
%       date: 09/28/07
%    purpose: 
%
function baseInfo(view)

% check arguments
if ~any(nargin == [1])
  help baseInfo
  return
end

% get base info
scanNum = viewGet(view,'curScan');
groupNum = viewGet(view,'curGroup');
baseDims = viewGet(view,'baseDims');
baseQform = viewGet(view,'baseqform');
baseSform = viewGet(view,'baseXform');
baseVolPermutation = viewGet(view,'baseVolPermutation');
baseVoxelSize = viewGet(view,'baseVoxelSize');
baseName = viewGet(view,'baseName');
baseCoordMap = viewGet(view,'baseCoordMap');
baseType = viewGet(view,'baseType');

% set parameters
paramsInfo = {{'baseName',baseName,'editable=0','The name of the base anatomy'},...
    {'voxelSize',baseVoxelSize,'editable=0','Voxel dimensions in mm'},...
    {'baseDims',baseDims,'editable=0','Dimensions of base anatomy'},...
    {'qform',baseQform,'editable=0','Qform matrix specifies the transformation to the scanner coordinate frame'},...
    {'sform',baseSform,'editable=0','Sform matrix is set by mrAlign and usually specifies the transformation to base coordinate system'},...
    {'baseType',baseType,'editable=0','Type of base. 0 = inplane. 1 = flat. 2 = surface'}};

% add baseCoordMap info for flat files
if baseType == 1
  paramsInfo{end+1} = {'flatDir',baseCoordMap.flatDir,'editable=0','Directory from which this flat map was originally created'};
  paramsInfo{end+1} = {'flatFileName',baseCoordMap.flatFileName,'editable=0','Name of original off file from which this flat map was created'};
  paramsInfo{end+1} = {'innerFileName',baseCoordMap.innerFileName,'editable=0','Name of inner mesh (aka gray matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'outerFileName',baseCoordMap.outerFileName,'editable=0','Name of outer mesh (aka white matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'curvFileName',baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this flat map was created'};
  paramsInfo{end+1} = {'anatFileName',baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this flat map was taken'};
  paramsInfo{end+1} = {'viewFlatOnSurface',[],'type=pushbutton','buttonString=View flat on surface','callback',@viewFlatOnSurface,'passParams=1','Click to view flat on the surface meshes'};
end

% add info for surfaces
if baseType == 2
  paramsInfo{end+1} = {'inner',baseCoordMap.innerFileName,'editable=0','Inner surface'};
  paramsInfo{end+1} = {'innerCoords',baseCoordMap.innerCoordsFileName,'editable=0','Inner surface coordinates'};
  paramsInfo{end+1} = {'outer',baseCoordMap.outerFileName,'editable=0','Outer surface'};
  paramsInfo{end+1} = {'outerCoords',baseCoordMap.outerCoordsFileName,'editable=0','Outer surface coordinates'};
  paramsInfo{end+1} = {'curv',baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this surface was created'};
  paramsInfo{end+1} = {'anatomy',baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this surface was taken'};
end

% bring up dialog
mrParamsDialog(paramsInfo,'Base anatomy information');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   viewFlatOnSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = viewFlatOnSurface(params,varargin);

thispwd = pwd;
if isdir(params.flatDir)
  cd(params.flatDir);
else
  mrWarnDlg(sprintf('Directory %s does not exist, please find the anatomy folder',params.flatDir));
  pathStr = uigetdir(mrGetPref('volumeDirectory','Find anatomy directory'));
  if pathStr == 0,return,end
  cd(pathStr);
end

mrFlatViewer(params.flatFileName,params.outerFileName,params.innerFileName,params.curvFileName,params.anatFileName);
cd(thispwd);