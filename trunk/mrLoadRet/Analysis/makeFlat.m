% makeFlat.m
%
%       $Id$	
%      usage: makeFlat()
%         by: eli merriam
%       date: 09/27/07
%    purpose: 
%
function retval = makeFlat(view, overlayNum, scan, x, y, s, roi) 


% check arguments
if ~any(nargin == [4 5 6 7])
  help makeFlat
  return
end

% get base info
baseCoordMap = viewGet(view,'baseCoordMap');
baseCoordMapPath = viewGet(view,'baseCoordMapPath');
params.startCoord = viewGet(view,'mouseDownBaseCoords');
baseType = viewGet(view,'baseType');

% parse the parameters
paramsInfo = {};
if ~isempty(baseCoordMap)
  if ~isfield(baseCoordMap, 'innerFileName')
    startPathStr = baseCoordMapPath;
    filterspec = {'*.off','SurfRelax off file';'*WM*.off', 'SurfRelax OFF WM file'; '*.*','All files'};
    title = 'Choose WM OFF file';    
    inner = getPathStrDialog(startPathStr,title,filterspec,'on');
    [params.flatDir baseCoordMap.innerFileName] = fileparts(inner{1});
    % Aborted
    if isempty(inner)
      disp(sprintf('(makeFlat) loading inner (WM) file aborted'));
      return
    end
  end
  % get the inner and outer surfaces. Unfortunately the naming convention
  % is slightly different for flat maps and volumes, so they have to be treated
  % differently
  if baseType == 1
    outerCoordsFileName = [stripext(baseCoordMap.outerFileName) '.off'];
    innerCoordsFileName = [stripext(baseCoordMap.innerFileName) '.off'];
  elseif baseType == 2
    if ~strcmp(baseCoordMap.outerCoordsFileName,'Same as surface')
      outerCoordsFileName = [stripext(baseCoordMap.outerCoordsFileName) '.off'];
    else
      outerCoordsFileName = [stripext(baseCoordMap.outerFileName) '.off'];
    end
    if ~strcmp(baseCoordMap.innerCoordsFileName,'Same as surface')
      innerCoordsFileName = [stripext(baseCoordMap.innerCoordsFileName) '.off'];
    else
      innerCoordsFileName = [stripext(baseCoordMap.innerFileName) '.off'];
    end
  end  

  % we can get all of the file names from the baseCoordMap
  paramsInfo{end+1} = {'flatDir', baseCoordMapPath,'editable=1','Directory from which this flat map was originally created'};
  paramsInfo{end+1} = {'innerCoordsFileName', innerCoordsFileName,'editable=1','Name of inner mesh (aka gray matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'outerCoordsFileName', outerCoordsFileName,'editable=1','Name of outer mesh (aka white matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'curvFileName', baseCoordMap.curvFileName,'editable=1','Name of curvature file from which this flat map was created'};
  paramsInfo{end+1} = {'anatFileName', baseCoordMap.anatFileName,'editable=1','Name of anatomy file from which the xform for this flat map was taken'};
  paramsInfo{end+1} = {'calcCurvFlag', 0, 'type=checkbox'};

else
  % we will have to build a list of file names...
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax OFF file';'*WM*.off','SurfRelax off gray matter file'; '*.*','All files'};
  title = 'Choose WM OFF file';
  innerCoordsFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
  % Aborted
  if isempty(innerCoordsFileName)
    disp(sprintf('(makeFlat) loading inner (WM) file aborted'));
    return
  end
  % make into a cell array
  innerCoordsFileName = cellArray(innerCoordsFileName);

  % only take one inner file for now
  if iscell(innerCoordsFileName)
    innerCoordsFileName = innerCoordsFileName{1};
  end

  % guess which hemisphere
  if regexp(innerCoordsFileName, 'left') 
    params.whichHemi = 'left';
  elseif regexp(innerCoordsFileName, 'right') 
    params.whichHemi = 'right';
  else
    disp('(makeFlat) Cannot guess hemisphere.  The user must choose')
    params.whichHemi = {'left', 'right'};
  end

  % contents of the current directory
  [params.flatDir params.innerCoordsFileName] = fileparts(innerCoordsFileName);
  dirContents = dir(params.flatDir);

  % guess the GM surface
  params.outerCoordsFileName = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'GM')) & (regexp(dirContents(i).name, '.off'))
      params.outerCoordsFileName{end+1} = dirContents(i).name;
    end
  end

  % need to guess the Curv surface
  % this may be empty, in which case, we'll calculate one
  params.curvFileName = {};
  params.calcCurvFlag = 1;
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexpi(dirContents(i).name, 'Curv')) & (regexp(dirContents(i).name, '.vff'))
      params.curvFileName{end+1} = dirContents(i).name;
      params.calcCurvFlag = 0;
    end
  end

  % guess the anatomy file
  % there may be many suitable anatomy files; list them all
  params.anatFileName = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, '.hdr'))
      params.anatFileName{end+1} = dirContents(i).name;
    end
  end

  % parse the parameters
  paramsInfo = {};
  paramsInfo{end+1} = {'flatDir', params.flatDir, 'editable=0', 'directory path for the inner .off surface file'};
  paramsInfo{end+1} = {'whichHemi', params.whichHemi, 'editable=0', 'the hemisphere that the patch comes from -- not really important'};
  paramsInfo{end+1} = {'innerCoordsFileName', sprintf('%s.off', stripext(params.innerCoordsFileName)), 'name of the surface at the white/gray boundary (i.e., the inner surface)'};
  paramsInfo{end+1} = {'outerCoordsFileName', params.outerCoordsFileName, 'name of the surface at the gray/pial boundary (i.e., the outer surface)'};
  if params.calcCurvFlag == 1;
    paramsInfo{end+1} = {'curvFileName', 'will calculate curvature on the fly', 'editable=0', 'name of the curvatue file'};
    paramsInfo{end+1} = {'calcCurvFlag', 1, 'type=checkbox', 'whether or not to recalculate the curvature on the file -- requires surffilt'};
  else
    paramsInfo{end+1} = {'curvFileName', params.curvFileName, 'name of the curvature file'};
    paramsInfo{end+1} = {'calcCurvFlag', 0, 'type=checkbox', 'whether or not to recaculate the curvature on the fly -- requires surffilt'};
  end
  paramsInfo{end+1} = {'anatFileName', params.anatFileName, 'base anatomy file, must be a valid nifti file in register with the session'};
end


paramsInfo{end+1} = {'startCoord', params.startCoord(1:3), 'start flattening from here'};
paramsInfo{end+1} = {'patchRadius', 75, 'round=1', 'incdec=[-5 5]', 'Flat patch radius in mm'};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]', 'the resolution of the flat patch -- a value of 2 doubles the resolution'};
paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'thresholding the surface makes the background two-tone (binary curvature)'};
paramsInfo{end+1} = {'flattenMethod', {'mrFlatMesh','surfRelax'},'use either surfRelax or mrFlatMesh'};

params = mrParamsDialog(paramsInfo, 'makeFlat', []);

if isempty(params)
  return;
end

% check for directory
if ~isdir(params.flatDir)
  params.flatDir = uigetdir(mrGetPref('volumeDirectory'),'Find anatomy directory');
  % user hits cancel
  if params.flatDir == 0
    return
  end
end

% get the nifti header for base anatomy
params.baseHdr = viewGet(view, 'basehdr');

% load the surfaces
[surf, params] = loadSurfHandler(params);

% find the vertex closests to the starting point
params.startVertex = dsearchn(surf.inner.vtcs, params.startCoord);

% set the name of the patch to cut
params.patchFileName = sprintf('%s_Patch_%i_%i_%i_Rad%i.off', ...
                               stripext(params.innerCoordsFileName), ...
                               round(params.startCoord(1)), round(params.startCoord(2)), round(params.startCoord(3)), ...
                               params.patchRadius);

% set the name of the new flat file
params.flatFileName = sprintf('%s_Flat_%i_%i_%i_Rad%i.off', ...
                              stripext(params.innerCoordsFileName), ...
                              round(params.startCoord(1)), round(params.startCoord(2)), round(params.startCoord(3)), ...
                              params.patchRadius);

if strcmp(params.flattenMethod, 'surfRelax')
  disp(sprintf('Flattening using SurfRelax'));
  myCutAndFlatten(params);
elseif strcmp(params.flattenMethod, 'mrFlatMesh')
  disp(sprintf('Flattening using the mrVista mrFlatMesh'));
  [surf, params] = runMrFlatMesh(surf, params);
  writePatchOFF(surf, params);
end

% make it into a MLR4 base anatomy
if isfile(fullfile(params.flatDir, params.flatFileName))

  base = loadFlatOFF(params);
  % get the vol2mag and vol2tal for the base from which the flat patch was defined
  base.vol2mag = params.vol2mag;
  base.vol2tal = params.vol2tal;

  % install it
  disp(sprintf('(makeFlat) installing new flat base anatomy: %s', params.flatFileName));
  if ~isempty(base)
    viewNum = viewGet(view, 'viewNum');
    viewSet(viewNum, 'newbase', base);
    refreshMLRDisplay(viewNum);
  end
end

% remove the 3d patch, b/c it isn't need for anything
if isfile(fullfile(params.flatDir, params.patchFileName))
  system(sprintf('rm -rf %s', fullfile(params.flatDir, params.patchFileName)));
end

return;



function myCutAndFlatten(params)

% set lib path
mylibs = getenv('DYLD_LIBRARY_PATH');
setenv('DYLD_LIBRARY_PATH', sprintf('%s:/Users/eli/src/TFI/sw/lib/', mylibs));

% check for the SurfRelax program called 'surfcut'
[statusCut,result] = system('surfcut');
[statusFlat,result] = system('FlattenSurface.tcl');
if any([statusCut statusFlat]) ~= 0
  disp(sprintf('(makeFlat) Could not run the SurfRelax program surfcut. Make sure that you have SurfRelax correctly installed on your system.'));
  return;
else

  % cut and flatten
  degenFlag = 1; 
  distanceInc = 0;
  
  while degenFlag ~= 0
    disp(sprintf('(makeFlat) Cutting patch with radius of %i mm', params.patchRadius+distanceInc))
    % cut the patch from the 3d mesh
    system(sprintf('surfcut -vertex %i -distance %i %s %s', ...
                   params.startVertex-1, params.patchRadius+distanceInc, ... 
                   fullfile(params.flatDir, params.innerCoordsFileName), ...
                   fullfile(params.flatDir, params.patchFileName)));
    disppercent(inf);
    
    % flatten the patch
    disppercent(-inf, sprintf('(makeFlat) Flattening surface'));
    [degenFlag result] = system(sprintf('FlattenSurface.tcl %s %s %s', ...
                                        fullfile(params.flatDir, params.outerCoordsFileName), ...
                                        fullfile(params.flatDir, params.patchFileName), ...
                                        fullfile(params.flatDir, params.flatFileName)));
    disppercent(inf);
    
    % if FlattenSurface failed, most likely b/c surfcut made a bad patch
    % increase the distance by one and try again.
    if (degenFlag ~= 0 ) && askuser('(makeFlat) Patch is degenerate, should I increase distance and reflatten?')
      distanceInc = distanceInc + 1;
    else
      degenFlag = 0;
      return;
    end

  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surf, params] = loadSurfHandler(params)
% we have already loaded the flat patch
% now load the rest of the surfaces

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = cbiReadNifti(fullfile(params.flatDir, params.anatFileName));
% get vol2tal and vol2mag from the anatomy file
matFileName = [stripext(params.anatFileName) '.mat'];
if(exist([params.flatDir '/' matFileName]))
  load([params.flatDir '/' matFileName]);
  [tf base] = isbase(base);
  params.vol2mag = base.vol2mag;
  params.vol2tal = base.vol2tal;
  clear base
else
  params.vol2mag = [];
  params.vol2tal = [];
end

% load the white matter surface
surf.inner = loadSurfOFF(fullfile(params.flatDir, params.innerCoordsFileName));
surf.inner = xformSurfaceWorld2Array(surf.inner, surf.anat.hdr);

% load the gray matter surface
surf.outer = loadSurfOFF(fullfile(params.flatDir, params.outerCoordsFileName));
surf.outer = xformSurfaceWorld2Array(surf.outer, surf.anat.hdr);

% % read in the curvature file
[surf.curv, hdr] = loadVFF(fullfile(params.flatDir, params.curvFileName));
% needs to be transposed to match the order of the vertices
surf.curv = surf.curv';           

% Extract permutation matrix to keep track of slice orientation.
surf.anat.permutationMatrix = getPermutationMatrix(surf.anat.hdr);

return


function[surf, params] = runMrFlatMesh(surf, params)

% project the surface out to an intermediate cortical deapth
corticalDepth = 0.5;
mesh.vertices = surf.inner.vtcs+corticalDepth*(surf.outer.vtcs-surf.inner.vtcs);

% create the mesh structure that mrFlatMesh expects
mesh.faceIndexList  = surf.inner.tris;
mesh.rgba           = surf.curv;
mesh.normal = surf.inner.vtcs - surf.outer.vtcs;

% run a modified version of the mrFlatMesh code
% this outputs and flattened surface
surf.flat = flattenSurfaceMFM(mesh, params.startCoord, params.patchRadius);

% old method of calculating the surface normals
% figure(999)
% Hp = patch('vertices', mesh.vertices, 'faces', mesh.faceIndexList);
% normals = get(Hp,'vertexnormals');
% close(999);

% we need to figure out whether the flattened patch has been flipped
% during flattening

patch2parent = surf.flat.vertsToUnique(surf.flat.insideNodes);
vIn   = surf.inner.vtcs(patch2parent,:);
vOut  = surf.outer.vtcs(patch2parent,:);
vFlat = surf.flat.locs2d;
f     = surf.flat.uniqueFaceIndexList;

% this is the command to view the patch, in 2D or 3D
%hp = patch('vertices', v, 'faces', f, 'facecolor','none','edgecolor','black');

% loop through all of the faces
disppercent(-inf,'Checking winding direction');
wrapDir = zeros(1,length(f));wrapDirFlat = zeros(1,length(f));
for iFace = 1:length(f);
  % grab a triangle for inner 3D suface
  triIn = vIn(f(iFace,:),:);
  % grab a triangle for outer 3D suface
  triOut = vOut(f(iFace,:),:);
  % calculate the vector normal to the center of the two triangles
  triNorm = (mean(triIn) - mean(triOut)) + mean(triIn);
  % this is a formula that takes the vertices of the triangle and
  % computes the winding direction relative to a fourth point, which
  % is in this case the normal.  i.e. do the vertices of the triangle
  % go in a CW or CCW direction with respect to the normal. If this
  % determinant is positive the direction is CW and if the determinant
  % is negative it is CCW. see wikipedia:
  % http://en.wikipedia.org/wiki/Orientation_(topology)
  wrapDir(iFace) = det([cat(2,triIn, [1 1 1]'); triNorm 1]);
  
  % now do the same for the triangles in the flatpatch
  triFlat = vFlat(f(iFace,:),:);
  % in the flat patch, the z-dimension is always 0
  triFlat(:,3) = 0;
  % we want to compute the winding direction from above the surface
  % (i.e., the direction that we are viewing the surface from)
  triFlatNorm = [0 0 1];
  % same formula as above
  wrapDirFlat(iFace) = det([cat(2,triFlat, [1 1 1]'); triFlatNorm 1]);
  disppercent(iFace/length(f));
end
disppercent(inf)

% now check to see if the winding directions for the flat patch and 3D
% surface are the same or different. Note that because of the
% flattening process, some of the triangles may switch winding
% direction. But whether we should view the patch from above or below
% is decided by which view produces the least number of mismatches in
% winding direction.
match =    sum( sign(wrapDir) == sign(wrapDirFlat) );
misMatch = sum( sign(wrapDir) ~= sign(wrapDirFlat) );

if misMatch > match
  disp(sprintf('(makeFlat) X-Y flipping the flat patch...'));
  surf.flat.locs2d = cat(2, surf.flat.locs2d(:,2), surf.flat.locs2d(:,1));
else
  disp(sprintf('(makeFlat) The patch is oriented properly, not going to flip it...'));
end


return;



% writeOFF.m
%
%      usage: writeOFF()
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writePatchOFF(surf, params)

% check arguments
if ~any(nargin == [ 0 1 2])
  help writeOFF
  return
end

% Vertices
vertices = [surf.flat.locs2d(:,1) surf.flat.locs2d(:,2)];
vertices = cat(1, vertices', zeros(1, length(vertices)));

% triangles(1) is number of vert/triangle: 3
% triangles(2:4) are the vertices of the triangles
% triangles(5) is color: 0
triangles =  surf.flat.uniqueFaceIndexList'-1;
triangles = cat(1, ones(1,length(triangles))*3, triangles, zeros(1,length(triangles),1));

patch2parent = surf.flat.vertsToUnique(surf.flat.insideNodes);

% write the OFF format file 
fid = fopen(fullfile(params.flatDir, params.flatFileName), 'w', 'ieee-be');
fprintf(fid, '#PATCH\n');
fprintf(fid, '#parent_surface=%s\n', fullfile(params.flatDir, params.innerCoordsFileName));
fprintf(fid, '#parent_dimensions=%i %i %i\n', surf.inner.Nvtcs,  surf.inner.Ntris, 1);
fprintf(fid, '#patch_dimensions=%i %i %i\n', length(vertices), length(triangles), 1 );
fprintf(fid, '#parent_vertex_indexes:\n');
for i=1:length(surf.flat.uniqueVertices)
  fprintf(fid, '#%i %i\n', i-1, patch2parent(i)-1);
end

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [size(vertices,2) size(triangles,2) 0], 'int32'); 

% Vertices
fwrite(fid, vertices, 'float32');

% Faces
fwrite(fid, triangles, 'int32');

% Close file
fclose(fid);

return;
