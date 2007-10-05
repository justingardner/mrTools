% makeFlat.m
%
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
  % we can get all of the file names from the baseCoordMap
  paramsInfo{end+1} = {'surfDir',baseCoordMap.flatDir,'editable=0','Directory from which this flat map was originally created'};
  paramsInfo{end+1} = {'flatFileName',baseCoordMap.flatFileName,'editable=0','Name of original off file from which this flat map was created'};
  paramsInfo{end+1} = {'innerFileName',baseCoordMap.innerFileName,'editable=0','Name of inner mesh (aka gray matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'outerFileName',baseCoordMap.outerFileName,'editable=0','Name of outer mesh (aka white matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'curvFileName',baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this flat map was created'};
  paramsInfo{end+1} = {'anatFileName',baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this flat map was taken'};
else
  % we will have to build a list of file names...
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax off file';'*lat.off','SurfRelax off flat file';'*.*','All files'};
  title = 'Choose flat OFF file';
  innerFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
  % Aborted
  if isempty(innerFileName)
    disp(sprintf('(makeFlat) loading innter file aborted'));
    return
  end
  % make into a cell array
  innerFileName = cellArray(innerFileName);

  % only take one inner file for now
  if iscell(innerFileName)
    innerFileName =   innerFileName{1};
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
  [params.surfDir params.innerFileName] = fileparts(innerFileName);
  dirContents = dir(params.surfDir);

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
paramsInfo{end+1} = {'surfDir', params.surfDir, 'editable=0', 'directory path for the inner .off surface file'};
paramsInfo{end+1} = {'whichHemi', params.whichHemi, 'editable=0', 'the hemisphere that the patch comes from -- not really important'};
paramsInfo{end+1} = {'innerFileName', strcat(params.innerFileName, '.off'), 'name of the surface at the white/gray boundary (i.e., the inner surface'};
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
params = mrParamsDialog(paramsInfo, 'makeFlat', []);

if isempty(params)
  return;
end

% load the surfaces
[surf, params] = loadSurfHandler(params);


% find the vertex closests to the starting point
params.startVertex = dsearchn(surf.inner.vtcs, params.startCoord);

params.patchFileName = sprintf('%s_Patch%i.off', stripext(params.innerFileName), params.startVertex);
params.flatFileName = sprintf('%s_Flat%i.off', stripext(params.innerFileName), params.startVertex);

myCutAndFlatten(params)



return;

function myCutAndFlatten(params)

% set lib path
setenv('DYLD_LIBRARY_PATH', '/Users/eli/src/TFI/sw/lib/');

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
                   fullfile(params.surfDir, params.innerFileName), ...
                   fullfile(params.surfDir, params.patchFileName)));
    disppercent(inf);
    
    % flatten the patch
    disppercent(-inf, sprintf('(makeFlat) Flattening surface'));
    [degenFlag result] = system(sprintf('FlattenSurface.tcl %s %s %s', ...
                                     fullfile(params.surfDir, params.innerFileName), ...
                                     fullfile(params.surfDir, params.patchFileName), ...
                                     fullfile(params.surfDir, params.flatFileName)));
    disppercent(inf);

    % if FlattenSurface failed, most likely b/c surfcut made a bad patch
    % increase the distance by one and try again.
    if degenFlag ~= 0
      disp(sprintf('(makeFlat) Patch is degenerate with distance of %i, increasing distance and recutting', params.patchRadius+distanceInc))
    end
    distanceInc = distanceInc + 1;
  end

  % remove the 3d patch b/c it isn't need for anything
  system(sprintf('rm -rf %s', fullfile(params.surfDir, params.patchFileName)));
  
end


% loadSurfOFF.m
%
%      usage: surf = loadSurfOFF(surffile)
%         by: eli merriam
%       date: 09/25/07
%    purpose: Loading an OFF binary surface into Matlab
%
% This code was ripped directly from Jonas Larsson's tfiReadOFF code
% the main difference is that this code returns a neat structure,
% rather than a bunch of variables
%
% INPUT ARGUMENTS:
% surffile - file name of surface in OFF binary format.
%
% OUTPUT ARGUMENTS:
% vtcs: Nvtcs x 3 matrix of vertex coordinates
% tris: Ntris x 3 matrix of triangle vertex indexes
%       NB!!! Vertex indexes are 1-offset for Matlab compatibility.
%       In the OFF surface vertex indexes are 0-offset.
%
% The surf structure consists of the following variables
% vtcs,tris,Nvtcs,Ntris,Nedges (nParent,nPatch,patch2parent)
%
function surf = loadSurfOFF(surffile)

% check arguments
if ~any(nargin == [1])
  help loadSurfOFF
  return
end

fid = fopen(surffile, 'r', 'ieee-be');
fl = fgetl(fid);
if (~strcmp(fl,'OFF BINARY'))
  if (strcmp(fl,'#PATCH'))

    % get the surface parent name
    fl = fgetl(fid);
    surf.parentSurfaceName = sscanf(fl,'#parent_surface=%s\n');
    
    % parent surface dimensions
    fl = fgetl(fid); 
    surf.nParent = sscanf(fl,'#parent_dimensions=%i %i %i\n');

    % patch surface dimensions
    fl = fgetl(fid); 
    surf.nPatch = sscanf(fl,'#patch_dimensions=%i %i %i\n');

    % parent vertex indexes tag
    fl = fgetl(fid); 
    
    % read patch info
    [a,c] = fscanf(fid,'#%i %i\n',[2 surf.nPatch(1)]);
    if (c ~= surf.nPatch(1)*2)
      disp(c)
      error('error reading file')
    end
    
    % matlab 1-offset
    surf.patch2parent = a'+1; 
    while (~strcmp(fl,'OFF BINARY'))
      fl= fgetl(fid);
    end
  else
    disp('WARNING!! Magic number (OFF BINARY) not found - this file may not be in the correct format!');
  end
end

[Ninfo, count] = fread(fid, 3, 'int32');
if (count~=3) error('error reading file!'); end
surf.Nvtcs  = Ninfo(1);
surf.Ntris  = Ninfo(2);
surf.Nedges = Ninfo(3);

[surf.vtcs, count] = fread(fid, surf.Nvtcs*3, 'float32');
if (count ~= 3*surf.Nvtcs) error('error reading file!'); end

[surf.tris, count] =  fread(fid, surf.Ntris*5, 'int32');
if (count ~= 5*surf.Ntris) error('error reading file!'); end

surf.vtcs = reshape(surf.vtcs, 3, surf.Nvtcs);
surf.tris = reshape(surf.tris, 5, surf.Ntris);

surf.vtcs = surf.vtcs';

%% ATTN ATTN ATTN %%
%% for some reason, I needed to add 2 to the vtcs to make them line up with the anatomy files.
% surf.vtcs = [surf.vtcs(:,2) surf.vtcs(:,1) surf.vtcs(:,3)]; % swaping x and y
surf.vtcs = surf.vtcs + 2;              % adding '2' here, but not sure why -epm



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surf, params] = loadSurfHandler(params)
% we have already loaded the flat patch
% now load the rest of the surfaces

% load the white matter surface
surf.inner = loadSurfOFF(fullfile(params.surfDir, params.innerFileName));

% load the gray matter surface
surf.outer = loadSurfOFF(fullfile(params.surfDir, params.outerFileName));

% % read in the curvature file
% [surf.curv, hdr] = tfiReadVFF(fullfile(params.surfDir, params.curvFileName));
% surf.curv = surf.curv';           % needs to be transposed;

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = cbiReadNifti(fullfile(params.surfDir, params.anatFileName));

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(surf.anat.hdr.qform44(1:3,1:3)));
surf.anat.permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

return

