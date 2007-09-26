% loadFlatOFF.m
%
%        $Id$	
%      usage: loadFlatOFF('flatPatch.off')
%         by: eli merriam
%       date: 09/11/07
%    purpose: load a surfRelax flatpatch, and everything that goes with it
%  
%    Three ways to use this code:
function base = loadFlatOFF(flatFileName)

% check arguments
if ~any(nargin == [0 1 2])
  help loadFlatOFF
  return
end

if ~exist('mrParamsDialog')
  disp(sprintf('(loadFlatOFF) You must have mrTools in your path to run this'));
  return
end

% init arguments
if nargin == 0
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax off flat file';'*.*','All files'};
  title = 'Choose flat OFF file';
  flatFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
  % make into a cell array
  flatFileName = cellArray(flatFileName);
  % Aborted
  if ieNotDefined('flatFileName')
    disp(sprintf('(loadFlatOFF) loading flat patch aborted'));
    return
  end
end

% only take one flat file for now
if iscell(flatFileName)
  flatFileName =   flatFileName{1};
end

% if we are passed in a string, then this is a file name
if isstr(flatFileName)
  basename = sprintf('%s.off', stripext(flatFileName));
  % check for file
  if isfile(flatFileName)
    % read file
    surf.flat = loadSurfOFF(flatFileName);
  else
    disp(sprintf('(loadFlatOFF) %s does not exist', flatFileName));
    return;
  end
  % check that it is indeed a patch
  if ~isfield(surf.flat, 'nPatch')
    disp(sprintf('(loadFlatOFF) %s is not a flat patch file', flatFileName))
    return;
  end
end

% guess which hemisphere
if regexp(flatFileName, 'left') 
  params.whichHemi = 'left';
elseif regexp(flatFileName, 'right') 
  params.whichHemi = 'right';
else
  disp('(loadFlatOFF) Cannot guess hemisphere.  The user must choose')
  params.whichHemi = {'left', 'right'};
end

% contents of the current directory
params.flatDir = fileparts(flatFileName);
dirContents = dir(params.flatDir);

% grab the white matter file name from the parent name
if exist(fullfile(params.flatDir, surf.flat.parentSurfaceName), 'file')
  params.wmFileName{1} = fullfile(params.flatDir, surf.flat.parentSurfaceName);
else
  % if it the parent surface doesn't exist, guess the WM name
  disp(sprintf('(loadFloatPatch) %s does not exist.  Guessing the WM surface name', surf.flat.parentSurfaceName))
  params.wmFileName = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'WM')) & (regexp(dirContents(i).name, '.off'))
      params.wmFileName{end+1} = fullfile(params.flatDir, dirContents(i).name);
    end
  end
end
% ATTN still need to add recursive directory searching


% guess the GM surface
params.gmFileName = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'GM')) & (regexp(dirContents(i).name, '.off'))
    params.gmFileName{end+1} = fullfile(params.flatDir, dirContents(i).name);
  end
end
% ATTN still need to add recursive directory searching


% need to guess the Curv surface
% this may be empty, in which case, we'll calculate one
params.curvFileName = {};
params.calcCurvFlag = 1;
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'Curv')) & (regexp(dirContents(i).name, '.vff'))
    params.curvFileName{end+1} = fullfile(params.flatDir, dirContents(i).name);
    params.calcCurvFlag = 0;
  end
end

% guess the anatomy file
% there may be many suitable anatomy files; list them all
params.anatFileName = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, '.hdr'))
    params.anatFileName{end+1} = fullfile(params.flatDir, dirContents(i).name);
  end
end

% parse the parameters
paramsInfo = {};
paramsInfo{end+1} = {'flatDir', params.flatDir, 'editable=0'};
paramsInfo{end+1} = {'whichHemi', params.whichHemi, 'editable=0'};
paramsInfo{end+1} = {'flatFileName', flatFileName, 'editable=0'};
paramsInfo{end+1} = {'wmFileName', params.wmFileName};
paramsInfo{end+1} = {'gmFileName', params.gmFileName};
if params.calcCurvFlag == 1;
  paramsInfo{end+1} = {'curvFileName', 'will calculate curvature on the fly', 'editable=0'};
  paramsInfo{end+1} = {'calcCurvFlag', 1, 'type=checkbox'};
else
  paramsInfo{end+1} = {'curvFileName', params.curvFileName};
  paramsInfo{end+1} = {'calcCurvFlag', 0, 'type=checkbox'};
end
paramsInfo{end+1} = {'anatFileName', params.anatFileName};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]'};
paramsInfo{end+1} = {'threshold', 0, 'type=checkbox'};
params = mrParamsDialog(paramsInfo, 'loadFlatPatch', []);

% load the rest of the surfaces
[surf, params] = loadSurfHandler(surf, params);

% calculate the base anatomy structure
base = calcFlatBase(surf, params);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surf, params] = loadSurfHandler(surf, params)
% we have already loaded the flat patch
% now load the rest of the surfaces

% load the white matter surface
surf.wm = loadSurfOFF(params.wmFileName);

% load the gray matter surface
surf.gm = loadSurfOFF(params.gmFileName);

% load the curvature
if params.calcCurvFlag==1
  params.baseName = stripext(params.wmFileName);
  params.curvFileName = sprintf('%s_Curv.vff', params.baseName);
  setenv('DYLD_LIBRARY_PATH', '/Users/eli/src/TFI/sw/lib/');
  % check for the SurfRelax program called 'surffilt'
  [status,result] = system('surffilt');
  if status ~= 0
    disp(sprintf('(loadFlatOFF) something is wrong with the surfrelax installation'))
    return;
  else
    % calculate the curvature
    disppercent(-inf, sprintf('(loadFlatOFF) calculating surface curvature'));
    system(sprintf('surffilt -mcurv -iter 1 %s %s', params.wmFileName, params.curvFileName));
    disppercent(inf);
  end
end
% read in the curvature file
[surf.curv, hdr] = tfiReadVFF(params.curvFileName);
surf.curv = surf.curv';           % needs to be transposed;

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = cbiReadNifti(params.anatFileName);

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(surf.anat.hdr.qform44(1:3,1:3)));
surf.anat.permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base] = calcFlatBase(surf, params)
% creates a flat patch base anatomy
% returns a base, which can then be either saved or installed in MLR

flat.whichInx  = surf.flat.patch2parent(:,2);
flat.locsFlat  = surf.flat.vtcs;
flat.curvature = surf.curv(flat.whichInx);
flat.locsInner = surf.wm.vtcs(flat.whichInx,:);
flat.locsOuter = surf.gm.vtcs(flat.whichInx,:);
flat.hdr       = surf.anat.hdr;

% this X-Y swaping only changes the orientation of the image
% isn't a crucial step
flat.locsFlat = [flat.locsFlat(:,2) flat.locsFlat(:,1) flat.locsFlat(:,3)];

% make the lookup table
% first get coordinates
xmin = min(flat.locsFlat(:,1));
xmax = max(flat.locsFlat(:,1));
ymin = min(flat.locsFlat(:,2));
ymax = max(flat.locsFlat(:,2));
x = xmin:(1/params.flatRes):xmax;
y = ymin:(1/params.flatRes):ymax;

% now we need to interp the curvature on to a regular
% grid. we do this by finding the x,y point that is closest
% to the one on the grid.
disppercent(-Inf, 'interpolating surface')
for i = 1:length(x)
  disppercent(i/length(x));
  for j = 1:length(y)
    % find nearest point in curvature
    dist = (flat.locsFlat(:,1)-x(i)).^2+(flat.locsFlat(:,2)-y(j)).^2;
    flat.pos(i,j) = first(find(min(dist)==dist));
    % points that are greater than a distance of 5 away are
    % probably not in the flat patch so mask them out
    if (min(dist) < 5)
      flat.mask(i,j) = 1;
      flat.baseCoordsInner(i,j,:) = flat.locsInner(flat.pos(i,j),:);
      flat.baseCoordsOuter(i,j,:) = flat.locsOuter(flat.pos(i,j),:);
      flat.map(i,j) = flat.curvature(flat.pos(i,j), :);
    else
      flat.mask(i,j) = 0;
      flat.baseCoordsInner(i,j,:) = [0 0 0];
      flat.baseCoordsOuter(i,j,:) = [0 0 0];
      flat.map(i,j) = 0;
    end
  end
end
disppercent(inf)

% now blur image
flat.blurMap(:,:) = blur(flat.map(:,:));
% threshold image
flat.median = median(flat.blurMap(:));
flat.thresholdMap(:,:) = (flat.blurMap(:,:)>median(flat.blurMap(:)))*0.5+0.25;
flat.thresholdMap = blur(flat.thresholdMap);
% mask out points not on map
flat.thresholdMap(~flat.mask(:)) = 0;

% now generate a base structure
clear base;
base.hdr = flat.hdr;
[pathstr, base.name] = fileparts(params.flatFileName);
base.permutationMatrix = surf.anat.permutationMatrix;
% base.params = params;

% load all the flat maps into the base. We
% need to make all the flat images have
% the same width and height.
if params.threshold
  base.data(:,:,1) = flat.thresholdMap;
  base.range = [0 1];
  base.clip = [0 1];
else
  flat.map(flat.map>1) = 1;
  flat.map(flat.map<-1) = -1;
  flat.map = 32*(flat.map+1);
  base.data(:,:,1) = flat.map;
end    

base.coordMap.coords(:,:,1,1) = flat.baseCoordsInner(:,:,2);
base.coordMap.coords(:,:,1,2) = flat.baseCoordsInner(:,:,1);
base.coordMap.coords(:,:,1,3) = flat.baseCoordsInner(:,:,3);

base.coordMap.innerCoords(:,:,1,1) = flat.baseCoordsInner(:,:,2);
base.coordMap.innerCoords(:,:,1,2) = flat.baseCoordsInner(:,:,1);
base.coordMap.innerCoords(:,:,1,3) = flat.baseCoordsInner(:,:,3);

base.coordMap.outerCoords(:,:,1,1) = flat.baseCoordsOuter(:,:,2);
base.coordMap.outerCoords(:,:,1,2) = flat.baseCoordsOuter(:,:,1);
base.coordMap.outerCoords(:,:,1,3) = flat.baseCoordsOuter(:,:,3);

base.coordMap.dims = flat.hdr.dim([2 3 4])';

base.range = [min(min(base.data)) max(max(base.data))];
base.clip = base.range;

% save(savename, 'base');
% v = viewSet(v,'newBase',base);


% ATTN: mrLoadRet ver 3 method
% % grid the data
% imSize = round([max(flat.locsFlat)]);
% base.map = NaN*ones(imSize(1:2));
% x = flat.locsFlat(:,1);
% y = flat.locsFlat(:,2);
% z = flat.curvature;
% yi = [1:imSize(1)]';
% xi = [1:imSize(2)];
% base.map(:,:) = griddata(x,y,z,xi,yi,'linear');
% base.map = mrUpSample(base.map);


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
surf.vtcs = [surf.vtcs(:,2) surf.vtcs(:,1) surf.vtcs(:,3)]; % swaping x and y
surf.vtcs = surf.vtcs + 2;              % adding '2' here, but not sure why -epm


% first entry is # of vertices per face (always 3); last entry is # of colors per face (only 0 allowed)
surf.tris = surf.tris(2:4,:); 
% 1-offset for matlab
surf.tris = surf.tris'+1; 

fclose(fid);

return;
