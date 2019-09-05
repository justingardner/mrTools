% mlrMakeROIFromSurfaceVertices.m
%
%      usage: roi = mlrMakeROIFromSurfaceVertices(vertices, surfaceNames)
%         by: justin gardner
%       date: 09/05/19
%    purpose: Makes a volume ROI (i.e. list of voxels) from surface vertices.
%             vertices is an array of nVertices that has a 1 or 0 for each vertex to indicate whether
%               it should be included in the roi. This program will check to make sure that the
%               number of vertices in the array matches the surface
%             surfaceNames is a struct similar to the one the mrSurfViewer returns with added field path for
%               where the surfaces should be located. 
%
function roi = mlrMakeROIFromSurfaceVertices(vertexNums,surfaceNames,varargin)

% check arguments
if nargin < 2
  help mlrMakeROIFromSurfaceVertices
  return
end

% empty roi to start with
roi = [];

% parse arguments
getArgs(varargin,{'corticalDepth',[0:0.1:1]});

% load the surfaces
base = importSurfaceOFF(surfaceNames);

% make sure that the number of vertices match
nVertices = length(vertexNums);
if nVertices ~= size(base.coordMap.innerVtcs,1)
  disp(sprintf('(mlrMakeROIFromSurfaceVertices) Number of vertices %i does not match the number in base %s (%i)',nVertices,fullfile(base.coordMap.path,base.coordMap.innerCoordsFileName),size(base.coordMap.innerVtcs,1)));
  return
end

% make an roi
roi.name = 'roi';
roi.voxelSize = base.hdr.pixdim(2:4);
% get xform - either sform if set, or qform
if ~isempty(base.hdr.sform44) && (base.hdr.sform_code ~= 0)
  roi.xform = base.hdr.sform44;
else
  roi.xform = base.hdr.qform44;
end
roi.createdBy = 'mlrMakeROIFromSurfaceVertices';
roi.surface.base = base;
roi.surface.vertexNums = vertexNums;
[tf roi] = isroi(roi);
if ~tf,keyboard,end

% get inner and outer coordinates
innerCoords = squeeze(base.coordMap.innerCoords(1,find(vertexNums),1,:));
outerCoords = squeeze(base.coordMap.outerCoords(1,find(vertexNums),1,:));

% cycle over cotical depth
coords = [];
for iCorticalDepth = corticalDepth
  % get coordinates for this cortical depth
  coords = [coords;round(iCorticalDepth * innerCoords + (1-iCorticalDepth)*outerCoords)];
end

% get unique coords
dims(1) = base.coordMap.dims(2);
dims(2) = base.coordMap.dims(1);
dims(3) = base.coordMap.dims(3);
coordsIndex = unique(sub2ind(dims,coords(:,1),coords(:,2),coords(:,3)));
[uniqueCoords(:,1) uniqueCoords(:,2) uniqueCoords(:,3)] = ind2sub(dims,coordsIndex);

% put these coords into the roi
roi.coords = uniqueCoords';


