% mlrLifePlugin
%
%        $Id:$ 
%      usage: mlrLifePlugin(action,<v>)
%         by: justin gardner & franco pestilli
%       date: 09/09/2014
%    purpose: Plugin function for LiFE
%
function retval = mlrLifePlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrLifePlugin) Need a valid view to install plugin'));
  else
    % add the menu item for Diffusion
    mlrAdjustGUI(v,'add','menu','Diffusion','/File/ROI');
    % add one to open Dt6 file
    mlrAdjustGUI(v,'add','menu','Load Dt6','/File/Diffusion/','Callback',@mlrLifeLoadDt6);
    % this menu item allows importation of fascicles as a surface
    mlrAdjustGUI(v,'add','menu','Import fascicles','/File/Base anatomy/Import surface','Callback',@mlrLifeImportFascicles);
    
    % make sure vista is in path
    mlrPath mrTools+vista
    
    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(mlrLifePlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeLoadDt6   %
%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeLoadDt6(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% check system
if ~mlrLifeSystemCheck,return,end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeImportFascicles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeImportFascicles(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% check system
if ~mlrLifeSystemCheck,return,end

% load the anatomy file which has the header
% FIX, FIX, FIX this is just a place-holder - as
% we should figure out a process for getting the right
% xform from the header.
disp(sprintf('(mlrLifePlugin) Choose diffusion weighted image from which the tracts were made'));
h = mlrImageHeaderLoad('.');
if isempty(h),return,end

% cooordinates of tracts are in magnet, so set xform appropriately
h.qform = eye(4);
h.sform = eye(4);
h.hdr.qform_code = 1;
h.hdr.sform_code = 1;

% Put up dialog to choose tractography file
startPathStr = '.';
filterspec = {'*.tck','Tractography file'};
title = 'Choose tractography file to load';
fgFilename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');
if isempty(fgFilename),return,end

% read the file
fg = fgRead(fgFilename);
if isempty(fg), return,end

% extract random subset of fibers
% FOR TESTING ONLY
nFibers = length(fg.fibers);
disp(sprintf('(mlrLifePlugin) Found %i fibers',nFibers));
params = mrParamsDialog({{'nFibersForTesting',5000,'Number of fibers to test on','min',1,'max',nFibers,'incdec',[-1000 1000]}},'Choose number of random fibers to load - this is for testing only');
if isempty(params),return,end
fg = fgExtract(fg, randsample(1:length(fg.fibers), params.nFibersForTesting),'keep');

% Build frames from the fascicles
[X, Y, Z] = mbaBuildFascicleFrame(fg.fibers);

% number of fascicles
nFascicles = length(X);

% initialize bounding box
xMin= inf;yMin = inf;zMin = inf;
xMax = -inf;yMax = -inf;zMax = -inf;

% Build a patch from the frame
nTotalVertices = 0;nTotalTris = 0;
disppercent(-inf,'(mlrLifeImportFascicles) Making fascicle surface');
for i = 1:nFascicles
  fasciclePatches{i} = surf2patch(X{i},Y{i},Z{i},'triangles');
  % compute how many vertices and tris we have all together
  nTotalVertices = size(fasciclePatches{i}.vertices,1) + nTotalVertices;
  nTotalTris = size(fasciclePatches{i}.faces,1) + nTotalTris;
  disppercent(i/nFascicles);
  % calculate bounding box
  xMin = min(xMin,min(fasciclePatches{i}.vertices(:,1)));
  xMax = max(xMax,max(fasciclePatches{i}.vertices(:,1)));
  yMin = min(yMin,min(fasciclePatches{i}.vertices(:,2)));
  yMax = max(yMax,max(fasciclePatches{i}.vertices(:,2)));
  zMin = min(zMin,min(fasciclePatches{i}.vertices(:,3)));
  zMax = max(zMax,max(fasciclePatches{i}.vertices(:,3)));
end
disppercent(inf);

% display bounding box of coordinates
disp(sprintf('(mlrLifePlugin) Bounding box x: %0.1f %0.1f y: %0.1f %0.1f z: %0.1f %0.1f',xMin,xMax,yMin,yMax,zMin,zMax));

% create an MLR base structure
fascicleBase.name = fixBadChars(fg.name);

% set this to be a surface
fascicleBase.type = 2; 

% set nifti header and permutation matrix in base structure
fascicleBase.hdr = mlrImageGetNiftiHeader(h);
fascicleBase.permutationMatrix = getPermutationMatrix(fascicleBase.hdr);
fascicleBase.vol2mag = fascicleBase.hdr.sform44;

% set path and names of files
fascicleBase.coordMap.path = fileparts(fgFilename);
fascicleBase.coordMap.innerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.innerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.curvFileName = '';
fascicleBase.coordMap.anatFileName = fgFilename;

% initalize fields
fascicleBase.data = [];
% set dimensions for coordMap
fascicleBase.coordMap.dims = h.dim;
fascicleBase.coordMap.innerCoords = nan(1,nTotalVertices,1,3);
fascicleBase.coordMap.innerVtcs = nan(nTotalVertices,3);
fascicleBase.coordMap.tris = nan(nTotalTris,3);
nRunningTotalVertices = 0;
nRunningTotalTris = 0;

% now put all fascicles vertices and triangles into one coordMap
disppercent(-inf,sprintf('(mlrLifePlugin) Converting %i fascicles',nFascicles));
for iFascicle = 1:nFascicles
  % number of vertices and triangles
  nVertices = size(fasciclePatches{iFascicle}.vertices,1);
  nTris = size(fasciclePatches{iFascicle}.faces,1);
  % the data which is the grayscale value to color the fascicles with (rand for now)
  fascicleBase.data = [fascicleBase.data rand(1,nVertices)];
  % convert vertices to a coord map which has one x,y,z element for each possible
  % location on the surface (which actually is just a 1xnVerticesx1 image)
  % add these vertices to existing vertices
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,1) = fasciclePatches{iFascicle}.vertices(:,1);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,2) = fasciclePatches{iFascicle}.vertices(:,2);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,3) = fasciclePatches{iFascicle}.vertices(:,3);
  % these are the display vertices which are the same as the coords
  fascicleBase.coordMap.innerVtcs(nRunningTotalVertices+1:nRunningTotalVertices+nVertices,:) = fasciclePatches{iFascicle}.vertices;
  % triangle faces
  fascicleBase.coordMap.tris(nRunningTotalTris+1:nRunningTotalTris+nTris,:) = (fasciclePatches{iFascicle}.faces + nRunningTotalVertices);
  % update runing totals
  nRunningTotalVertices = nRunningTotalVertices + nVertices;
  nRunningTotalTris= nRunningTotalTris + nTris;
  disppercent(iFascicle/nFascicles);
end
disppercent(inf);

% copy the inner to outer since they are all the same for fascicles
fascicleBase.coordMap.outerCoords = fascicleBase.coordMap.innerCoords;
fascicleBase.coordMap.outerVtcs = fascicleBase.coordMap.innerVtcs;

% save the individual fascicles so we can rebuild later
fascicleBase.fascicles.patches = fasciclePatches;
fascicleBase.fascicles.n = length(fasciclePatches);
fascicleBase.fascicles.nTotalVertices = nTotalVertices;
fascicleBase.fascicles.nTotalTris = nTotalTris;

% make it a full base
[tf fascicleBase] = isbase(fascicleBase);

% add it to the view
v = viewSet(v,'newBase',fascicleBase);

% refresh the display to draw it
refreshMLRDisplay(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrLifeSystemCheck   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrLifeSystemCheck 

tf = false;
% check if we have fgRead from vistasoft
if exist('fgRead') ~= 2
  mrWarnDlg(sprintf('(mlrLifePlugin) You need to have vista installed to access function fgRead.\ngit clone https://github.com/vistalab/vistasoft.git vistasoft'));
  return
end

% check if we have mbaBuildFascicleFrame from mba
if exist('mbaBuildFascicleFrame') ~= 2
  mrWarnDlg(sprintf('(mlrLifePlugin) You need to have mba installed to access function mbaBuildFascicleFrame.\ngit clone https://github.com/francopestilli/mba.git mba'));
  return
end
 tf = true;