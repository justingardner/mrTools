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
  help mrFlatFromCoord
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
scanXform = viewGet(view,'scanXform',scan);
scanVoxelSize = viewGet(view,'scanVoxelSize',scan);

params.startCoord = xformROIcoords([x;y;s;1], inv(baseSform)*scanXform, baseVoxelSize, baseVoxelSize);


% parse the parameters
paramsInfo = {};
if ~isempty(baseCoordMap)
  if ~isfield(baseCoordMap, 'innerFileName')
    startPathStr = baseCoordMap.flatDir;
    filterspec = {'*.off','SurfRelax off file';'*WM*.off', 'SurfRelax OFF WM file'; '*.*','All files'};
    title = 'Choose WM OFF file';    
    innerFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
    [params.flatDir baseCoordMap.innerFileName] = fileparts(innerFileName{1});
    % Aborted
    if isempty(innerFileName)
      disp(sprintf('(makeFlat) loading inner (WM) file aborted'));
      return
    end
  end
  % we can get all of the file names from the baseCoordMap
  paramsInfo{end+1} = {'flatDir', baseCoordMap.flatDir,'editable=0','Directory from which this flat map was originally created'};
  paramsInfo{end+1} = {'flatFileName', baseCoordMap.flatFileName,'editable=0','Name of original off file from which this flat map was created'};
  paramsInfo{end+1} = {'innerFileName', baseCoordMap.innerFileName,'editable=0','Name of inner mesh (aka gray matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'outerFileName', baseCoordMap.outerFileName,'editable=0','Name of outer mesh (aka white matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'curvFileName', baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this flat map was created'};
  paramsInfo{end+1} = {'anatFileName', baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this flat map was taken'};
  paramsInfo{end+1} = {'calcCurvFlag', 0, 'editable=0'};

else
  % we will have to build a list of file names...
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax OFF file';'*WM*.off','SurfRelax off gray matter file'; '*.*','All files'};
  title = 'Choose flat OFF file';
  innerFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
  % Aborted
  if isempty(innerFileName)
    disp(sprintf('(makeFlat) loading inner (WM) file aborted'));
    return
  end
  % make into a cell array
  innerFileName = cellArray(innerFileName);

  % only take one inner file for now
  if iscell(innerFileName)
    innerFileName = innerFileName{1};
  end

  % guess which hemisphere
  if regexp(innerFileName, 'left') 
    params.whichHemi = 'left';
  elseif regexp(innerFileName, 'right') 
    params.whichHemi = 'right';
  else
    disp('(makeFlat) Cannot guess hemisphere.  The user must choose')
    params.whichHemi = {'left', 'right'};
  end

  % contents of the current directory
  [params.flatDir params.innerFileName] = fileparts(innerFileName);
  dirContents = dir(params.flatDir);

  % guess the GM surface
  params.outerFileName = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'GM')) & (regexp(dirContents(i).name, '.off'))
      params.outerFileName{end+1} = dirContents(i).name;
    end
  end

  % need to guess the Curv surface
  % this may be empty, in which case, we'll calculate one
  params.curvFileName = {};
  params.calcCurvFlag = 1;
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'Curv')) & (regexp(dirContents(i).name, '.vff'))
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
  paramsInfo{end+1} = {'innerFileName', params.innerFileName, 'name of the surface at the white/gray boundary (i.e., the inner surface)'};
  paramsInfo{end+1} = {'outerFileName', params.outerFileName, 'name of the surface at the gray/pial boundary (i.e., the outer surface)'};
  if params.calcCurvFlag == 1;
    paramsInfo{end+1} = {'curvFileName', 'will calculate curvature on the fly', 'editable=0', 'name of the curvatue file'};
    paramsInfo{end+1} = {'calcCurvFlag', 1, 'type=checkbox', 'whether or not to recalculate the curvature on the file -- requires surffilt'};
  else
    paramsInfo{end+1} = {'curvFileName', params.curvFileName, 'name of the curvature file'};
    paramsInfo{end+1} = {'calcCurvFlag', 0, 'type=checkbox', 'whether or not to recaculate the curvature on the fly -- requires surffilt'};
  end
  paramsInfo{end+1} = {'anatFileName', params.anatFileName, 'base anatomy file, must be a valid nifti file in register with the session'};
end


paramsInfo{end+1} = {'startCoord', params.startCoord(1:3)', 'start flattening from here'};
paramsInfo{end+1} = {'patchRadius', 50, 'round=1', 'incdec=[-5 5]', 'flat patch radius'};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]', 'the resolution of the flat patch -- a value of 2 doubles the resolution'};
paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'thresholding the surface makes the background two-tone (binary curvature)'};
paramsInfo{end+1} = {'flattenMethod', {'surfRelax', 'mrFlatMesh'},'use either surfRelax or mrFlatMesh'};

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

% load the surfaces
[surf, params] = loadSurfHandler(params);

% find the vertex closests to the starting point
% need to swap X and Y
params.startCoord = [params.startCoord(2) params.startCoord(1) params.startCoord(3)];
params.startVertex = dsearchn(surf.inner.vtcs, params.startCoord);

% set the name of the patch to cut
params.patchFileName = sprintf('%s_Patch_%i_%i_%i_Rad%i.off', ...
                               stripext(params.innerFileName), ...
                               params.startCoord(1), params.startCoord(2), params.startCoord(3), ...
                               params.patchRadius);

% set the name of the new flat file
params.flatFileName = sprintf('%s_Flat_%i_%i_%i_Rad%i.off', ...
                              stripext(params.innerFileName), ...
                              params.startCoord(1), params.startCoord(2), params.startCoord(3), ...
                              params.patchRadius);

if params.flattenMethod == 'surfRelax'
  disp(sprintf('Flattening using SurfRelax'))
  myCutAndFlatten(params);
elseif params.flattenMethod == 'mrFlatMesh'
  disp(sprintf('Flattening using the mrVista mrFlatMesh'))
  [surf, params] = runMrFlatMesh(surf, params)
  writeOFF(surf, params);
end

% make it into a MLR4 base anatomy
if isfile(fullfile(params.flatDir, params.flatFileName))

  base = loadFlatOFF(params);

  % install it
  disp(sprintf('(makeFlat) installing new flat base anatomy: %s', params.flatFileName));
  if ~isempty(base)
    viewNum = 1;                          % dunno how to figure out the right view num
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
                   params.startVertex, params.patchRadius+distanceInc, ... 
                   fullfile(params.flatDir, params.innerFileName), ...
                   fullfile(params.flatDir, params.patchFileName)));
    disppercent(inf);
    
    % flatten the patch
    disppercent(-inf, sprintf('(makeFlat) Flattening surface'));
    [degenFlag result] = system(sprintf('FlattenSurface.tcl %s %s %s', ...
                                        fullfile(params.flatDir, params.innerFileName), ...
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

% load the white matter surface
surf.inner = loadSurfOFF(fullfile(params.flatDir, params.innerFileName));

% load the gray matter surface
surf.outer = loadSurfOFF(fullfile(params.flatDir, params.outerFileName));

% % read in the curvature file
[surf.curv, hdr] = loadVFF(fullfile(params.flatDir, params.curvFileName));
surf.curv = surf.curv';           % needs to be transposed;

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = cbiReadNifti(fullfile(params.flatDir, params.anatFileName));

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(surf.anat.hdr.qform44(1:3,1:3)));
surf.anat.permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

return


function[surf, params] = runMrFlatMesh(surf, params)

mesh.vertices       = surf.inner.vtcs;
mesh.faceIndexList  = surf.inner.tris;
mesh.rgba           = surf.curv;

figure(999)
Hp = patch('vertices', mesh.vertices, 'faces', mesh.faceIndexList);
mesh.normal = get(Hp,'vertexnormals');
close(999);

surf.flat = myUnfoldMeshFromGUI(mesh, params.startCoord, params.patchRadius);


return;



% writeOFF.m
%
%      usage: writeOFF()
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writeOFF(surf, params)

% check arguments
if ~any(nargin == [ 0 1 2])
  help writeOFF
  return
end

% Vertices
vertices = [surf.flat.locs2d(:,2) surf.flat.locs2d(:,1)] - 2;
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
fprintf(fid, '#parent_surface=%s\n', fullfile(params.flatDir, params.innerFileName));
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
