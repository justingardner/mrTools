% loadFlatOFF.m
%
%        $Id$	
%      usage: loadFlatOFF('flatPatch.off')
%         by: eli merriam
%       date: 09/11/07
%    purpose: load a surfRelax flatpatch, and everything that goes with it
%
function retval = loadFlatOFF(event,viewNum)

% check arguments
if ~any(nargin == [1 2 3])
  help loadFlatPatch
  return
end

if ~exist('mrParamsDialog')
  disp(sprintf('(loadFlatOFF) You must have mrTools in your path to run this'));
  return
end

global gSurf;

% init arguments
if nargin == 1
  % if we are passed in a structure then this is a button callback
  if isstruct(event)
    viewNum = event.viewNum;
    event = event.event;
    retval = [];
    % otherwise it is init event
  elseif isstr(event)
    filename = sprintf('%s.off', stripext(event));
    event = 'init';
    % check for file
    if isfile(filename)
      % read file
      surf.flat = loadOffSurface(filename);
    else
      disp(sprintf('(loadFlatOFF) %s does not exist', filename));
      return;
    end
    % check that it is indeed a patch
    if ~isfield(surf.flat, 'nPatch')
      disp(sprintf('(loadFlatOFF) %s is not a flat patch file', filename))
      return;
    end
  end
end

switch (event)
  case 'init'
    initHandler(filename);
  case 'end'
    endHandler(viewNum);
  case 'dispWM'
    gSurf{viewNum}.whichSurf = 'wm';
    dispSurfHandler(viewNum);
  case 'dispGM'
    gSurf{viewNum}.whichSurf = 'gm';
    dispSurfHandler(viewNum);
  case 'dispPatch'
    gSurf{viewNum}.whichSurf = 'patch';
    dispSurfHandler(viewNum);
  case 'dispAnat'
    dispAnatHandler(viewNum);
  case 'reload'
    loadSurfHandler(viewNum);
  case 'save'
    saveHandler(viewNum);
  case 'close'
    closeHandler(viewNum);
end

return;

function[] = initHandler(filename)
global gSurf;

% get figure handles
viewNum = smartfig('SurfaceViewer');
gSurf{viewNum} = [];
gSurf{viewNum}.fig(1) = viewNum;
gSurf{viewNum}.shutdown = 0;
gSurf{viewNum}.filename = filename;

% set info for callback
gsurf{viewNum}.viewNum = viewNum;

set(gSurf{viewNum}.fig(1),'DeleteFcn',sprintf('loadFlatPatch(''end'',%i)',viewNum));

gSurf{viewNum}.patch = loadOffSurface(filename);

% guess which hemisphere
if regexp(filename, 'left') 
  params.whichHemi = 'left';
elseif regexp(filename, 'right') 
  params.whichHemi = 'right';
else
  disp('(loadFlatOFF) Cannot guess hemisphere.  The user must choose')
  params.whichHemi = {'left', 'right'};
end

% contents of the current directory
dirContents = dir;

% grab the white matter file name from the parent name
if exist(gSurf{viewNum}.patch.parentSurfaceName, 'file')
  params.wmFile{1} = gSurf{viewNum}.patch.parentSurfaceName;
else
  % if it the parent surface doesn't exist, guess the WM name
  disp(sprintf('(loadFloatPatch) %s does not exist.  Guessing the WM surface name', gSurf{viewNum}.patch.parentSurfaceName))
  params.wmFile = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'WM'))
      params.wmFile{end+1} = dirContents(i).name;
    end
  end
end

% need to guess the GM surface
params.gmFile = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'GM'))
    params.gmFile{end+1} = dirContents(i).name;
  end
end

% need to guess the Curv surface
% this may be empty, in which case, we'll calculate one
params.curvFile = {};
gSurf{viewNum}.calcCurvFlag = 1;
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'Curv'))
    params.curvFile{end+1} = dirContents(i).name;
    gSurf{viewNum}.calcCurvFlag = 0;
  end
end

% guess the anatomy file
% there may be many suitable anatomy files; list them all
params.anatFile = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, '.hdr'))
    params.anatFile{end+1} = dirContents(i).name;
  end
end

% parse the parameters
callbackArg.viewNum = viewNum;
paramsInfo = {};
paramsInfo{end+1} = {'thisDir', pwd, 'editable=0'};
paramsInfo{end+1} = {'whichHemi', params.whichHemi, 'editable=0'};

paramsInfo{end+1} = {'wmFile', params.wmFile};
callbackArg.event = 'dispWM';
paramsInfo{end+1} = {'show',0,'type=pushbutton','callback=loadFlatPatch','buttonString=Display WM surface','callbackArg',callbackArg};

paramsInfo{end+1} = {'gmFile', params.gmFile};
callbackArg.event = 'dispGM';
paramsInfo{end+1} = {'show',0,'type=pushbutton','callback=loadFlatPatch','buttonString=Display GM surface','callbackArg',callbackArg};

if isempty(params.curvFile)
  paramsInfo{end+1} = {'curvFile', 'will calculate curvature on the fly', 'editable=0'};
else
  paramsInfo{end+1} = {'curvFile', params.curvFile};
end

paramsInfo{end+1} = {'anatFile', params.anatFile};
callbackArg.event = 'dispAnat';
paramsInfo{end+1} = {'show',0,'type=pushbutton','callback=loadFlatPatch','buttonString=Display anatomy file','callbackArg',callbackArg};


paramsInfo{end+1} = {'patchOverlay', 0, 'type=checkbox'};
callbackArg.event = 'dispPatch';
paramsInfo{end+1} = {'show',0,'type=pushbutton','callback=loadFlatPatch','buttonString=display Patch','callbackArg',callbackArg};

callbackArg.event = 'reload';
paramsInfo{end+1} = {'show',0,'type=pushbutton','callback=loadFlatPatch','buttonString=reload surfaces','callbackArg',callbackArg};

paramsInfo{end+1} = {'azimuth', 90,'Azimuth','round=1','minmax=[-inf inf]','incdec=[-5 5]'};
paramsInfo{end+1} = {'elevation', 0,'Elevation','round=1','minmax=[-inf inf]','incdec=[-5 5]'};
paramsInfo{end+1} = {'slice', 100,'slice number','round=1','minmax=[1 inf]','incdec=[-1 1]'};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]'};
paramsInfo{end+1} = {'threshold', 0, 'type=checkbox'};
callbackArg.event = 'save';
paramsInfo{end+1} = {'save',0,'type=pushbutton','callback=loadFlatPatch','buttonString=Save patch structure','callbackArg',callbackArg,'Save out the patch structure'};

% get the parameters from the user
[gSurf{viewNum}.controlFig gSurf{viewNum}.params] = mrParamsDialog(paramsInfo, 'loadFlatPatch', [], @surfControls, viewNum);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for handling controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surfControls(params, viewNum)
global gSurf;
gSurf{viewNum}.params = params;

figure(gSurf{viewNum}.fig(1));
view([params.azimuth params.elevation]);

% if params.patchOverlay
%   set(gSurf{viewNum}.Hp(2), 'visible', 'on');
% end

if isfield(gSurf{viewNum}, 'oldparams')
  if gSurf{viewNum}.params.slice ~= gSurf{viewNum}.oldparams.slice
    dispAnatHandler(viewNum);
  end
end
gSurf{viewNum}.oldparams = params;
drawnow;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for loading the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadSurfHandler(viewNum)
global gSurf;
params = gSurf{viewNum}.params;

% load the white matter surface
gSurf{viewNum}.wm = loadOffSurface(params.wmFile);

% load the gray matter surface
gSurf{viewNum}.gm = loadOffSurface(params.gmFile);

% load the curvature
if gSurf{1}.calcCurvFlag==1
  params.baseName = stripext(params.wmFile);
  params.curvFile = sprintf('%s_Curv.vff', params.baseName);
  setenv('DYLD_LIBRARY_PATH', '/Users/eli/src/TFI/sw/lib/');
  % check for the SurfRelax program called 'surffilt'
  [status,result] = system('surffilt');
  if status ~= 0
    disp(sprintf('(loadFlatOFF) something is wrong with the surfrelax installation'))
    return;
  else
    % calculate the curvature
    disppercent(-inf, sprintf('(loadFlatOFF) calculating surface curvature'));
    system(sprintf('surffilt -mcurv -iter 1 %s %s', params.wmFile, params.curvFile));
    disppercent(inf);
  end
end
% read in the curvature file
[gSurf{viewNum}.wm.curv, hdr] = tfiReadVFF(params.curvFile);
gSurf{viewNum}.wm.curv = gSurf{viewNum}.wm.curv';           % needs to be transposed;

% read in the anatomy file
[gSurf{viewNum}.anat.data   gSurf{viewNum}.anat.hdr] = cbiReadNifti(params.anatFile);

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(gSurf{viewNum}.anat.hdr.qform44(1:3,1:3)));
gSurf{viewNum}.anat.permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for displaying the anatomy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = dispAnatHandler(viewNum)

global gSurf;
params = gSurf{viewNum}.params;
figure(gSurf{viewNum}.fig(1));
cla;

% check to make sure the anatomy data has been loaded
if ~isfield(gSurf{viewNum}, 'anat')
  loadSurfHandler(viewNum)
end

% display a slice of the anatomy image
imagesc(gSurf{viewNum}.anat.data(:,:,params.slice));
colormap(gray);
axis image;
axis off;

% on display patch vtc if button pressed
if params.patchOverlay
  whichInx = gSurf{viewNum}.patch.patch2parent(:,2);
  wmNodes = gSurf{viewNum}.wm.vtcs(whichInx,:);
  gmNodes = gSurf{viewNum}.gm.vtcs(whichInx,:);
else
  wmNodes = (gSurf{viewNum}.wm.vtcs);
  gmNodes = (gSurf{viewNum}.gm.vtcs);
end

% plot the nodes, displaying both deep and superficial surfaces
wmNodes = wmNodes( find( round(wmNodes(:,3))==params.slice), : );
plot(wmNodes(:,1), wmNodes(:,2), 'b.', 'markersize', 1);

gmNodes = gmNodes( find( round(gmNodes(:,3))==params.slice), : );
plot(gmNodes(:,1), gmNodes(:,2), 'r.', 'markersize', 1);

view([0 90]);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for displaying the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = dispSurfHandler(viewNum)

global gSurf;
params = gSurf{viewNum}.params;

if ~isfield(gSurf{viewNum}, 'wm')
  loadSurfHandler(viewNum);
end

disp(sprintf('Displaying %s surface', gSurf{viewNum}.whichSurf));

vtcs = eval(sprintf('gSurf{viewNum}.%s.vtcs', gSurf{viewNum}.whichSurf));
tris = eval(sprintf('gSurf{viewNum}.%s.tris', gSurf{viewNum}.whichSurf));
if strcmp(gSurf{viewNum}.whichSurf, 'patch')
  curv = gSurf{viewNum}.wm.curv(gSurf{viewNum}.patch.patch2parent(:,2));
else
  curv = gSurf{viewNum}.wm.curv;
end

figure(gSurf{viewNum}.fig(1));
cla;

gSurf{viewNum}.Hp(1) = patch('vertices', vtcs, 'faces', tris, ...
                             'FaceVertexCData', curv, ...
                             'facecolor', 'interp', ...
                             'edgecolor', 'none');
set(gca,'CLim',[-1.2 1.2]);

overlay = NaN(length(curv),3);
overlay(gSurf{viewNum}.patch.patch2parent(:,2),1) = 1;
overlay(gSurf{viewNum}.patch.patch2parent(:,2),2:3) = 0;
gSurf{viewNum}.Hp(2) = patch('vertices', vtcs, 'faces', tris, ...
                             'FaceVertexCData', overlay, 'FaceVertexAlphaData', overlay(:,1)*.1, ...
                             'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',.6);
set(gSurf{viewNum}.Hp(2), 'visible', 'off');

daspect([1,1,1]);
camproj('orthographic');                % 'perspective'

% orient the patch
if strcmp(gSurf{viewNum}.whichSurf, 'patch')
  gSurf{viewNum}.params.elevation = 90;
end
view([gSurf{viewNum}.params.azimuth gSurf{viewNum}.params.elevation]);

colormap(gray);
material dull;
lighting phong;

axis off
rotate3d;
gSurf{viewNum}.display = 1;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for loading a single surface in OFF format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surf = loadOffSurface(surffile);
% [vtcs,tris,Nvtcs,Ntris,Nedges]=tfiReadOFF( surffile );
% surf = tfiReadOFF(surffile)
% [vtcs,tris,Nvtcs,Ntris,Nedges,nParent,nPatch,patch2parent]
%
% PURPOSE: Loading an OFF binary surface into Matlab
%
% INPUT ARGUMENTS:
% surffile - file name of surface in OFF binary format.
%
% OUTPUT ARGUMENTS:
% vtcs: Nvtcs x 3 matrix of vertex coordinates
% tris: Ntris x 3 matrix of triangle vertex indexes
%       NB!!! Vertex indexes are 1-offset for Matlab compatibility.
%       In the OFF surface vertex indexes are 0-offset.

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

surf.vtcs = surf.vtcs + 2;
surf.vtcs = surf.vtcs';
surf.vtcs = [surf.vtcs(:,2) surf.vtcs(:,1) surf.vtcs(:,3)];

% first entry is # of vertices per face (always 3); last entry is # of colors per face (only 0 allowed)
surf.tris = surf.tris(2:4,:); 
% 1-offset for matlab
surf.tris = surf.tris'+1; 

fclose(fid);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrInterrogator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(viewNum)

global gSurf;
% only run this if it is not being called by someonw else
if ~gSurf{viewNum}.shutdown
  gSurf{viewNum}.shutdown = 1;
  if ishandle(gSurf{viewNum}.controlFig)
    close(gSurf{viewNum}.controlFig);
  end
  if ishandle(gSurf{viewNum}.fig(1))
    close(gSurf{viewNum}.fig(1));
  end
  gSurf{viewNum}.shutdown = 0;
end
return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveHandler(viewNum)

global gSurf;
params = gSurf{viewNum}.params;


% now get the filename to save to
thispath = pwd;
if isfield(gSurf{viewNum},'savepath')
  cd(gSurf{viewNum}.savepath);
end

[savename savepath] = uiputfile({'*.mat','matlab files (*.mat)'},'Save as');

cd(thispath);

% check for user cancel
if isequal(savename,0) || isequal(savepath,0)
  return
end

% otherwise set the default savepath
gSurf{viewNum}.savepath = savepath;

% set proper extension
savename = fullfile(savepath,sprintf('%s.mat',stripext(savename)));
disp(sprintf('Saving %s',savename));

% need to save the following fields: curvature, gLocs2d, gLocs3d, pLocs3d, hdr
flat.whichInx = gSurf{viewNum}.patch.patch2parent(:,2);
flat.curvature = gSurf{viewNum}.wm.curv(flat.whichInx);
flat.locsInner = gSurf{viewNum}.wm.vtcs(flat.whichInx,:);
flat.locsOuter = gSurf{viewNum}.gm.vtcs(flat.whichInx,:);
flat.locsFlat  = gSurf{viewNum}.patch.vtcs;
flat.hdr = gSurf{viewNum}.anat.hdr;

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
h = mrWaitBar(0,sprintf('Creating flat image for %s', 'foobar'));
for i = 1:length(x)
  mrWaitBar(i/length(x),h);
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
mrCloseDlg(h);

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
base.name = stripext(gSurf{1}.filename);
base.permutationMatrix = gSurf{viewNum}.anat.permutationMatrix;

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
base.coordMap.dims = flat.hdr.dim([2 3 4])';

base.range = [min(min(base.data)) max(max(base.data))];
base.clip = base.range;

save(savename, 'base');
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
