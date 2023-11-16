%
% function thisView = makeROIsFromSurfaces(thisView,params,<'justGetParams'>)
%
%   Purpose: make ROIs based on cortical surfaces. Two non-overlapping ROIs will be made:
%               - one white-matter ROI with voxels within the inner (gw/wm boundary) surface,
%               - one grey matter ROI with voxels between the innner and outer (pial) surfaces
%   
%   Usage:   [~,params] = makeROIsFromSurfaces(thisView,[],'justGetParams') % to get default parameters
%            % then set the parameters and use them to create the ROIs, e.g.:
%            params.roiSpace = 'current scan'
%            thisView = makeROIsFromSurfaces(thisView,params); % make the ROIs and add them to the view
%            refreshMLRDisplay(thisView); % display the ROIs in the GUI
%
%   Params: - names: name of the ROIs. If left empty, the names will be self-explanatory
%           - colors: color of the ROIs
%           - surfaceName: name of the base containing the surfaces (default: current base or next available surface base in the view)
%           - roiSpace: space in which the roi coordinates will be saved: options are 'surface base': the
%               base corresponding to the surface coordinates (default), 'current base' or 'current scan'.
%           - tolerance: this is a parameter that is directly passed on to function inpolyhedron.m, and determines how far the center of
%               a voxel can be from (outside) the suface and still be considered to be located within an closed surface (default = 0).
%               It's unclear in what units this parameter is expressed, but a larger value will result in more voxels being included.
%               With the default value, voxels that have at least 50% of their volume within the surface should be included. However
%               this is not exact because inpolyhedron.m does not take the voxels size into account. It can also take negative values.
%

function [thisView, params] = makeROIsFromSurfaces(thisView,params,varargin)

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
% if ieNotDefined('defaultParams'),defaultParams = 0;end

if ~ismember(nargin,[1 2 3])
  help('makeROIsFromSurfaces');
  return
end

if ieNotDefined('thisView')
  mrWarnDlg('(makeSphereROIatCoords) No valid mrLoadRet view was provided');
  return
end

% set default parameters
if ieNotDefined('params')
  params = struct;
end

if fieldIsNotDefined(params,'colors')
  params.colors = {'green','orange'}; % color of the ROIs
end
if fieldIsNotDefined(params,'roiSpace')
  params.roiSpace = 'surface base'; % space in which the coordinates are specified: options are 'surface base': the base corresponding to the surface coordinates (default), 'current base' or 'current scan'.
end
if fieldIsNotDefined(params,'surfaceName')
  params.surfaceName = '';
  for iBase = [viewGet(thisView,'curbase'):viewGet(thisView,'numbase') viewGet(thisView,'curbase')-1:-1:1]
    if viewGet(thisView,'baseType',iBase) == 2
      params.surfaceName = viewGet(thisView,'baseName',iBase);
      break;
    end
  end
  if isempty(params.surfaceName)
    mrWarnDlg('(makeROIsFromSurfaces) Could not find any surface base in view');
  end
end
if fieldIsNotDefined(params,'names') % name of the ROIs.
  [~,roiBaseName] = fileparts(params.surfaceName); % remove .off extension
  roiBaseName = fixBadChars(roiBaseName);
  params.names{1} = [roiBaseName '_gm'];
  params.names{2} = [roiBaseName '_wm'];
end
if fieldIsNotDefined(params,'tolerance')
  params.tolerance = 0;
end

if justGetParams
  return;
end

surfBaseNum = viewGet(thisView,'basenum',params.surfaceName);
if isempty(surfBaseNum)
  mrWarnDlg(sprintf('(makeSphereROIatCoords) Could not find base %s in view.',params.surfaceName));
  return
elseif viewGet(thisView,'baseType',surfBaseNum) ~= 2
  mrWarnDlg(sprintf('(makeSphereROIatCoords) %s is not a surface base.',params.surfaceName));
  return
end

surfBaseCoordMap = viewGet(thisView,'baseCoordMap',surfBaseNum); % get the surface vertex coords and other surface info

% get the xform from the surface base space to the ROI space, as weel as voxel size in that space
switch(params.roiSpace)
  case 'surface base' % if the ROI and surface space are the same
    surf2roi = eye(4); % the xform is the identity matric
    roi.xform = viewGet(thisView,'basexform',surfBaseNum);
    roiVolDims = surfBaseCoordMap.dims;
    roi.voxelSize = viewGet(thisView,'baseVoxelSize',surfBaseNum);

  case 'current base'
    surf2roi = viewGet(thisView,'base2base',[],surfBaseNum);
    roi.xform = viewGet(thisView,'basexform');
    roi.voxelSize = viewGet(thisView,'baseVoxelSize'); % check that this works for flat maps?
    switch(viewGet(thisView,'baseType'))
      case 0
        roiVolDims = viewGet(thisView,'baseDims');
      case {1,2}
        baseCoordMap = viewGet(thisView,'baseCoordMap');
        roiVolDims = baseCoordMap.dims;
        if viewGet(thisView,'baseType') ==  1
          keyboard; % this function hasn't been tested with flat maps
        end
    end

  case 'current scan'
    surf2roi = viewGet(thisView,'base2scan',[],[],surfBaseNum);
    roi.xform = viewGet(thisView,'scanxform');
    roiVolDims = viewGet(thisView,'scanDims');
    roi.voxelSize = viewGet(thisView,'scanVoxelSize');

  otherwise
    mrWarnDlg(sprintf('(makeSphereROIatCoords) Unknow roiSpace parameter %s.',params.roiSpace));
    return
end

% Identify voxels within the inner surface (WM) using inpolyhedron function (obtained from Matlab Exchange).
% First, determine whether surface normals will be pointing out or in after transformation to the ROI space
% I think that Freesurfer (surfRelax?)'s convention is that they are pointing out, but if the transformation
% to ROI space does not preserve orientation (e.g. functional space can have left/right flipped), then that
% convention ends up reversed. I'm not sure how to determine the original mesh's convention, or even if it's
% possible (I suspect not) but, if we assume that the normals point up in the original mesh, then we can look
% at the determinant of the xfrom from surface to ROI space to see whether it preserves orientation or not.
% A negative determinant means that the convention has been flipped.
flipNormals = det(surf2roi)<0;
% now convert the surfaces to ROI space and find the voxels within the inner surface
nVtcs = size(surfBaseCoordMap.innerCoords,2);
innerCoords = surf2roi * [permute(surfBaseCoordMap.innerCoords,[4 2 1 3]); ones(1,nVtcs)]; % convert inner surface coordinates to ROI space
innerSurf.vertices = innerCoords(1:3,:)';
innerSurf.faces = surfBaseCoordMap.tris;
wmCoordsMask = inpolyhedron(innerSurf,1:roiVolDims(2),1:roiVolDims(1),1:roiVolDims(3),'flipnormals',flipNormals,'tol',params.tolerance);
wmCoordsMask = permute(wmCoordsMask,[2 1 3]); % permute X and Y because inpolyhedron uses the meshgrid convention

% do the same for the outer surface
outerCoords = surf2roi * [permute(surfBaseCoordMap.outerCoords,[4 2 1 3]); ones(1,nVtcs)];
outerSurf.vertices = outerCoords(1:3,:)';
outerSurf.faces = surfBaseCoordMap.tris;
gmCoordsMask = inpolyhedron(outerSurf,1:roiVolDims(2),1:roiVolDims(1),1:roiVolDims(3),'flipnormals',flipNormals,'tol',params.tolerance);
gmCoordsMask = permute(gmCoordsMask,[2 1 3]); % permute X and Y because inpolyhedron uses the meshgrid convention
gmCoordsMask(wmCoordsMask) = false; % subtract WM voxels to get on GM voxels

% make ROI structures and add them to the view
roi.createdBy = 'makeROIsFromSurfaces';
[~, roi] = isroi(roi); % add all the optional fields
% gm ROI
roi.color = params.colors{1};
roi.name = params.names{1};
roi.coords = ones(4,nnz(gmCoordsMask));
[roi.coords(1,:), roi.coords(2,:), roi.coords(3,:)] = ind2sub(roiVolDims,find(gmCoordsMask));
thisView= viewSet(thisView,'newROI',roi);
% wm ROI
roi.color = params.colors{2};
roi.name = params.names{2};
roi.coords = ones(4,nnz(wmCoordsMask));
[roi.coords(1,:), roi.coords(2,:), roi.coords(3,:)] = ind2sub(roiVolDims,find(wmCoordsMask));
thisView= viewSet(thisView,'newROI',roi);

