% loadSurface.m
%
%      usage: v = loadSurface(v)
%         by: justin gardner
%       date: 10/24/07
%    purpose: 
%
function base = loadSurface

% check arguments
if ~any(nargin == [0])
  help loadSurface
  return
end
base = [];

% Open dialog box to have user choose the file
startPathStr = mrGetPref('volumeDirectory');
filterspec = {'*.off','Off Surface file (*.off)'};
title = 'Choose outer surface file';
pathStr = getPathStrDialog(startPathStr,title,filterspec,'off');

% Aborted
if ieNotDefined('pathStr'),return,end

% get surface name using mrSurfViewer
[filepath filename] = fileparts(pathStr);
thispwd = pwd;
cd(filepath);
params = mrSurfViewer(filename);

% Aborted
if isempty(params),cd(thispwd);return,end

% Create the base
base.hdr = cbiReadNiftiHeader(params.anatomy);
base.name = filename;
base.permutationMatrix = getPermutationMatrix(base.hdr);
base.data(1,:,1) = loadVFF(params.curv);
base.range = [-1.5 1.5];
base.clip = [-1.5 1.5];

% get the vol2mag and vol2tal fields from the volume anatomy
% do it using a subfunction so don't confuse the two 'base' structure variables
anatomyMatFile = sprintf('%s.mat',stripext(params.anatomy));
base.vol2tal = getBaseField(anatomyMatFile,'vol2tal');
base.vol2mag = getBaseField(anatomyMatFile,'vol2mag');

% load the inner and outer coords
innerSurface = loadSurfOFF(params.innerSurface);
outerSurface = loadSurfOFF(params.outerSurface);
if strcmp(params.innerCoords,'Same as surface')
  inner = innerSurface;
else
  inner = loadSurfOFF(params.innerCoords);
end
if strcmp(params.outerCoords,'Same as surface')
  outer = outerSurface;
else
  outer = loadSurfOFF(params.outerCoords);
end

% save names
base.coordMap.innerFileName = params.innerSurface;
base.coordMap.innerCoordsFileName = params.innerCoords;
base.coordMap.outerFileName = params.outerSurface;
base.coordMap.outerCoordsFileName = params.outerCoords;
base.coordMap.curvFileName = params.curv;
base.coordMap.anatFileName = params.anatomy;
% save coords
base.coordMap.innerCoords(1,:,1,1)  = inner.vtcs(:,2);
base.coordMap.innerCoords(1,:,1,2)  = inner.vtcs(:,1);
base.coordMap.innerCoords(1,:,1,3)  = inner.vtcs(:,3);
base.coordMap.innerVtcs = innerSurface.vtcs;
base.coordMap.tris = innerSurface.tris;
base.coordMap.outer = params.outerSurface;
base.coordMap.outerCoords(1,:,1,1)  = outer.vtcs(:,2);
base.coordMap.outerCoords(1,:,1,2)  = outer.vtcs(:,1);
base.coordMap.outerCoords(1,:,1,3)  = outer.vtcs(:,3);
base.coordMap.outerVtcs = outerSurface.vtcs;
base.coordMap.coords = base.coordMap.innerCoords;
base.coordMap.dims = base.hdr.dim([3 2 4])';
base.type = 2;

cd(thispwd);


function val = getBaseField(matFilename,fieldname)
if ~exist(matFilename,'file') % if the anatomy doesn't have these fields set
  val = [];
else % if the volume anatomy has the fields, use them
  load(matFilename);
  eval(sprintf('val = base.%s;',fieldname))
end


  