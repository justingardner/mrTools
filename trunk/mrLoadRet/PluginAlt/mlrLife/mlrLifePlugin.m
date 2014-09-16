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
     disp(sprintf('(DefaultPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    % this installs a new menu item called 'Select Plugins' under /Edit/ROI with the
    % separator turned on above it. It sets the callback to selectPlugins defined below.
    mlrAdjustGUI(v,'add','menu','LiFE','/File/ROI','Callback',@mlrLifeFileOpen);
%    mlrAdjustGUI(v,'add','menu','Import','/File/LiFE/','Callback',@mlrLifeFileOpen);

    mlrAdjustGUI(v,'add','menu','Import fascicles','/File/Base anatomy/Import surface','Callback',@mlrLifeImportFascicles);
    % This is a command that could be used to install some default interrogators
    %mlrAdjustGUI(v,'add','interrogator',{'eventRelatedPlot','glmContrastPlot'});

    % This is a command that could be used to install some default colormaps
    % that will show up when you do /Edit/Overlay
    %mlrAdjustGUI(v,'add','colormap','gray');

    % This is a command that could be used to set a property of an existing menu item
    %mlrAdjustGUI(v,'set','Plots/Mean Time Series','Separator','on');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(DefaultPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeFileOpen    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeFileOpen(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% bring up file dialog to select dt6 file
startPathStr = pwd;
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% now bring up selection dialog
paramsInfo = {{'groupName','DTI'},...
	      {'analysisName','DTI'}};
params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

% convert t1 into a MLR anatomy
% Set required structure fields (additional fields are set to default
% values when viewSet calls isbase).
base.name = stripext(stripext(getLastDir(dt.files.t1)));
base.data = double(t1.img);
% create a header
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;

% set nifti header and permutation matrix in base structure
base.hdr = mlrImageGetNiftiHeader(h);
base.permutationMatrix = getPermutationMatrix(base.hdr);

% make it a full base
[tf base] = isbase(base);

% add it to the view
v = viewSet(v,'newBase',base);

% make group
v = viewSet(v,'newGroup',params.groupName);
v = viewSet(v,'curGroup',params.groupName);

% make nifti header for DTI
h = [];hdr = [];
h.dim = size(dt.dt6);
[tf h] = mlrImageIsHeader(h);
h.pixdim = dt.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = dt.xformToAcpc;
hdr = mlrImageGetNiftiHeader(h);

% save tseries
saveNewTSeries(v,dt.dt6,[],hdr);
scanNum = viewGet(v,'nScans');

% get the overlay
refreshMLRDisplay(v);

% convert RGB into a colormap
% make table
numValues = 8;
values = 0:(256/(numValues-1)):256;
rgbTable = [];
for rValue = 1:length(values)
  for gValue = 1:length(values)
    for bValue = 1:length(values)
      rgbTable(end+1,:) = [values(rValue) values(gValue) values(bValue)];
    end
  end
end
rgbTable = rgbTable(1:2:end,:);

% convert overlay to this table
disppercent(-inf,'Converting vectorRGB to index');
for x = 1:size(o.vectorRGB,1)
  disppercent(x/size(o.vectorRGB,1));
  for y = 1:size(o.vectorRGB,2)
    for z = 1:size(o.vectorRGB,3)
      colorVal = double(squeeze(o.vectorRGB(x,y,z,:)));
      [dummy colorIndex] = min(sum((rgbTable - repmat(colorVal',size(rgbTable,1),1)).^2,2));
      colormapRGB(x,y,z) = colorIndex;
    end
  end
end
disppercent(inf);

% display the overlay
rgbTable = rgbTable/256.0;
% FIX, FIX, FIX Create a proper analysis structure with all the overlays.
mrDispOverlay(colormapRGB,scanNum,params.groupName,v,'overlayName=principleDiffusionDirection','saveName=diffusionAnalysis','cmap',rgbTable);
a = viewGet(v,'curAnalysis');
mrDispOverlay(double(o.b0),scanNum,params.groupName,v,'overlayName=b0','saveName=diffusionAnalysis');

dateString = datestr(now);
b0.name = 'b)';
b0.groupName = params.groupName;
b0.function = '';
b0.reconcileFunction = 'defaultReconcileParams';
b0.data = double(o.b0);
b0.date = dateString;
b0.params = cell(1,viewGet(v,'nScans'));
b0.range = [0 1];
b0.clip = [0 1];
% colormap is made with a little bit less on the dark end
b0.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';
r2.mergeFunction = 'defaultMergeParams';


% create analysis
dti.name = params.analysisName;
dti.type = 'dti';
dti.groupName = params.groupName;
dti.function = '';
dti.reconcileFunction = 'defaultReconcileParams';
dti.mergeFunction = 'defaultMergeParams';
dti.guiFunction = '';
dti.params = params;
dti.overlays = b0;
dti.curOverlay = 1;
dti.date = dateString;
v = viewSet(v,'newAnalysis',dti);

% Save it
saveAnalysis(v,dti.name);


refreshMLRDisplay(v);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeImportFascicles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeImportFascicles(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Load the fascicles from disk
[fg, fgFilename] = fgRead;
if isempty(fg), return,end

% bring up file dialog to select dt6 file
% this is so that we can get the T1 file
% FIX, FIX, FIX - should find some way to automatically
% load the T1 file without having the user have to go find it.
startPathStr = fullfile(fileparts(fgFilename),'..');
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% Get the header from the T1 so that we have the correct sform etc
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;
h.sform = eye(4);
% Build frames from the fascicles
[X, Y, Z] = mbaBuildFascicleFrame(fg.fibers);

% number of fascicles
nFascicles = length(X);

% Build a patch from the frame
for i = 1:nFascicles
  fasciclePatches{i} = surf2patch(X{i},Y{i},Z{i},'triangles');
end

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
fascicleBase.coordMap.anatFileName = dt.files.t1;

% initalize fields
fascicleBase.data = [];
% set dimensions for coordMap
fascicleBase.coordMap.dims = h.dim;
fascicleBase.coordMap.innerVtcs = [];
fascicleBase.coordMap.tris = [];
nTotalVertices = 0;

disppercent(-inf,sprintf('(mlrLifePlugin) Converting %i fascicles',nFascicles));
for iFascicle = 1:nFascicles
  nVertices = size(fasciclePatches{iFascicle}.vertices,1);
  fascicleBase.data = [fascicleBase.data rand(1,nVertices)];
  % convert vertices to a coord map which has one x,y,z element for each possible
  % location on the surface (which actually is just a 1xnVerticesx1 image)
  % add these vertices to existing vertices
  fascicleBase.coordMap.innerCoords(1,nTotalVertices+1:nTotalVertices+nVertices,1,1) = fasciclePatches{iFascicle}.vertices(:,1);
  fascicleBase.coordMap.innerCoords(1,nTotalVertices+1:nTotalVertices+nVertices,1,2) = fasciclePatches{iFascicle}.vertices(:,2);
  fascicleBase.coordMap.innerCoords(1,nTotalVertices+1:nTotalVertices+nVertices,1,3) = fasciclePatches{iFascicle}.vertices(:,3);
  % inner and outer coords are the same
  fascicleBase.coordMap.outerCoords = fascicleBase.coordMap.innerCoords;
  % these are the display vertices which are the same as the coords
  fascicleBase.coordMap.innerVtcs = [fascicleBase.coordMap.innerVtcs ; fasciclePatches{iFascicle}.vertices];
  fascicleBase.coordMap.outerVtcs = fascicleBase.coordMap.innerVtcs;
  fascicleBase.coordMap.tris = [fascicleBase.coordMap.tris ; (fasciclePatches{iFascicle}.faces + nTotalVertices)];
  % update nTotalVertices
  nTotalVertices = nTotalVertices + nVertices;
  disppercent(iFascicle/nFascicles);
end
disppercent(inf);

% make it a full base
[tf fascicleBase] = isbase(fascicleBase);

% add it to the view
v = viewSet(v,'newBase',fascicleBase);

refreshMLRDisplay(v);

keyboard