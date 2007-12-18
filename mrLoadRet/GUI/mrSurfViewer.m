% mrSurfViewer.m
%
%       $Id$	
%      usage: mrSurfViewer(outerSurface,outerCoords,innerSurface,innerCoords,curv,anat)
%         by: justin gardner
%       date: 10/09/07
%    purpose: 
%
function retval = mrSurfViewer(outerSurface,outerCoords,innerSurface,innerCoords,curv,anat)

% check arguments
if ~any(nargin == [1 2 3 4 5 6])
  help mrSurfViewer
  return
end
if nargout == 1
  retval = [];
end

% if passed in a string check to see if
% it needs an extension
if isstr(outerSurface)
  if isfile(sprintf('%s.off',stripext(outerSurface)));
    outerSurface = sprintf('%s.off',stripext(outerSurface));
  end
end

% see how we are being called
if (nargin == 1) && isstr(outerSurface) && ~isfile(outerSurface)
  event = outerSurface;
elseif (nargin == 1) && isstr(outerSurface)
  if ~isempty(strfind(outerSurface,'WM'))
    innerSurface{1} = outerSurface;
    clear outerSurface;
    outerSurface{1} = sprintf('%sGM.off',strtok(stripext(innerSurface{1}),'WM'));
  elseif ~isempty(strfind(outerSurface,'GM'))
    innerSurface{1} = sprintf('%sWM.off',strtok(stripext(outerSurface),'GM'));
    clear outerSurface;
    outerSurface{1} = sprintf('%sGM.off',strtok(stripext(innerSurface{1}),'WM'));
  else
    innerSurface{1} = outerSurface;
    clear outerSurface;
    outerSurface{1} = innerSurface{1};
  end
  event = 'init';
  innerCoords = {};
  outerCoords = {};
  curv = {};
  anat = {};
else
  event = 'init';
  % set defaults
  
  if ieNotDefined('outerCoords'),outerCoords = {};end
  if ieNotDefined('innerSurface'),innerSurface = {};end
  if ieNotDefined('innerCoords'),innerCoords = {};end
  if ieNotDefined('curv'),curv = {};end
  if ieNotDefined('anat'),anat = {};end
  % make everybody a cell array
  outerSurface = cellArray(outerSurface);
  outerCoords = cellArray(outerCoords);
  innerSurface = cellArray(innerCoords);
  innerCoords = cellArray(innerCoords);
  curv = cellArray(curv);
  anat = cellArray(anat);
end

switch (event)
 case 'init'
  retval = initHandler(outerSurface,outerCoords,innerSurface,innerCoords,curv,anat);
 case {'vSlider','hSlider'}
  sliderHandler;
 case {'edit'}
  editHandler;
 otherwise
  disp(sprintf('(mrSurfViewer) Could not find surf file %s',outerSurface));
end

%%%%%%%%%%%%%%%%%%%%%%
%%   init handler   %%
%%%%%%%%%%%%%%%%%%%%%%
function retval = initHandler(outerSurface,outerCoords,innerSurface,innerCoords,curv,anat);

retval = [];
global gSurfViewer;
gSurfViewer.mismatchWarning = 0;
filepath = '';

disppercent(-inf,'(mrSurfaceView) Loading surfaces');

% load the surface
gSurfViewer.outerSurface = loadSurfOFF(sprintf('%s.off',stripext(outerSurface{1})));
if isempty(gSurfViewer.outerSurface)
  disp(sprintf('(mrSurfaceViewer) %s is not a surface file',outer{1}));
  return
end

% look for surfaces with same number of vertices/triangles
surfDir = dir('*.off');
allSurfaces = {};
for i = 1:length(surfDir)
  % load just the header
  surf = myLoadSurface(surfDir(i).name,filepath,1);
  if ~isempty(surf)
    allSurfaces{end+1} = surfDir(i).name;
  end
end
outerCoords = putOnTopOfList('Same as surface',allSurfaces);
outerCoords {end+1} = 'Find file';

% check for an inner surface
if isempty(innerSurface)
  innerSurface = sprintf('%sWM.off',strtok(stripext(outerSurface{1}),'GM'));
  if ismember(innerSurface,allSurfaces)
    innerSurface = putOnTopOfList(innerSurface,allSurfaces);
    gSurfViewer.innerSurface = myLoadSurface(innerSurface{1});
  else
    innerSurface = putOnTopOfList('None',allSurfaces);
  end
elseif length(innerSurface) == 1
  innerSurface = putOnTopOfList(innerSurface{1},allSurfaces);
end
innerCoords = putOnTopOfList('Same as surface',allSurfaces);
innerCoords {end+1} = 'Find file';
innerSurface{end+1} = 'Find file';
gSurfViewer.outerCoords = [];
gSurfViewer.innerCoords = [];

% load any vff file
curv = {};
curvDir = dir('*.vff');
for i = 1:length(curvDir)
  if ~any(strcmp(curvDir(i).name,curv))
    % check length of file matches our patch
    vffhdr = myLoadCurvature(curvDir(i).name, filepath, 1);
    if ~isempty(vffhdr)
      curv{end+1} = curvDir(i).name;
    end
  end
end

% check to see if we have any possible curvatures
if isempty(curv)
  mrWarnDlg(sprintf('(mrSurfViewer) Could not find a matching .vff curvature file. You will need to use the SurfRelax tool surffilt to generate one. It is usually invoked by doing: surffilt -mcurv -iter 2 %s %s_curv.vff',outerSurface{1},stripext(outerSurface{1})));
  return
end

% if we didn't load anything then quit
for i = 1:length(curv)
  gSurfViewer.curv = myLoadCurvature(curv{i},filepath);
  if ~isempty(gSurfViewer.curv),break,end
end
if isempty(gSurfViewer.curv)
  return
end
disppercent(inf);
curv{end+1} = 'Find file';

% guess any nifti file for anatomy
anatDir = dir('*.hdr');
for i = 1:length(anatDir)
  if ~any(strcmp(anatDir(i).name,anat))
    anat{end+1} = anatDir(i).name;
  end
end
% check for 'canonical hdr'
anatCanonicalDir = dir('../*.hdr');
for i= 1:length(anatCanonicalDir)
  if find(strcmp(anatCanonicalDir(i).name,anat))
    anat = putOnTopOfList(anatCanonicalDir(i).name,anat);
  end
end
% if no files found, then ask user to go look for it
if isempty(anat)
  [filename, pathname] = uigetfile({'*.hdr','Nifti file (*.hdr)'},'Find 3D anatomy file');
  if isempty(filename),return,end
  if ~isCurrentPath(pathname)
    filename = fullfile(pathname,filename);
  end
  anat{1} = filename;
end
if isfile(anat{1})
  [gSurfViewer.anat.data gSurfViewer.anat.hdr] = cbiReadNifti(anat{1});
else
  gSurfViewer.anat = [];
end
anat{end+1} = 'Find file';
gSurfViewer.mismatchWarning = 1;

% select the window
gSurfViewer.f = selectGraphWin;
set(gSurfViewer.f,'renderer','OpenGL');

% positions on figure
figLeft = 10;figBottom = 10;
sliderWidth = 20;sliderLength = 200;spacer = 10;
editWidth = 40;editHeight = 20;

% set up horizontal and vertical slider
gSurfViewer.hSliders.v = uicontrol('Style','slider','Position',[figLeft figBottom+sliderWidth sliderWidth sliderLength],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrSurfViewer(''vSlider'')','TooltipString','Rotate around y-axis');
gSurfViewer.hSliders.vText = uicontrol('Style','Edit','Position',[figLeft figBottom+sliderWidth+sliderLength+spacer editWidth editHeight],'Callback','mrSurfViewer(''edit'')','String','0','HorizontalAlignment','Center');
gSurfViewer.hSliders.h = uicontrol('Style','slider','Position',[figLeft+sliderWidth figBottom sliderLength sliderWidth],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrSurfViewer(''hSlider'')','TooltipString','Rotate around z-axis');
gSurfViewer.hSliders.hText = uicontrol('Style','Edit','Position',[figLeft+sliderLength+3*spacer figBottom editWidth editHeight],'Callback','mrSurfViewer(''edit'')','String','0');

% set they we are viewing white matter
gSurfViewer.whichSurface = 1;

% and display surface
dispSurface;
setViewAngle(0,0);

% set up the parameters
paramsInfo = {};
gSurfViewer.guiloc.filenames = 1;
% Now give choice of viewing gray or white
gSurfViewer.whichSurfaceTypes = {'Outer (Gray matter) surface','Outer (Gray matter) coords','Inner (White matter) surface','Inner (White matter) coords','3D Anatomy'};
paramsInfo{end+1} = {'whichSurface',gSurfViewer.whichSurfaceTypes,'callback',@whichSurfaceCallback,'Choose which surface to view'};

paramsInfo{end+1} = {'outerSurface',outerSurface,'The outer (gray matter) surface that will be shown as the outer surface','callback',@switchFile,'callbackArg=outerSurface'};
paramsInfo{end+1} = {'outerCoords',outerCoords,'The outer (gray matter) surface that contains the coordinates for the display surface. Usually this is the same as the outerSurface, but if you have an inflated surface, the display surface and the coords surface will be different','callback',@switchFile,'callbackArg=outerCoords'};
paramsInfo{end+1} = {'innerSurface',innerSurface,'The inner (white matter) display surface','callback',@switchFile,'callbackArg=innerSurface'};
paramsInfo{end+1} = {'innerCoords',innerCoords,'The inner (white matter) coords surface','callback',@switchFile,'callbackArg=innerCoords'};
paramsInfo{end+1} = {'curv',curv,'The curvature file','callback',@switchFile,'callbackArg=curv'};
paramsInfo{end+1} = {'anatomy',anat,'The 3D anatomy file','callback',@switchAnatomy};


% put up dialog
params = mrParamsDialog(paramsInfo,'View surface');

close(gSurfViewer.f);

if isempty(params)
  retval = [];
else
  retval = params;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   whichSurfaceCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichSurfaceCallback(params)

global gSurfViewer;

% get which surface to draw
lastWhichSurface = gSurfViewer.whichSurface;

refreshFlatViewer(find(strcmp(params.whichSurface,gSurfViewer.whichSurfaceTypes)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   refreshFlatViewer   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshFlatViewer(whichSurface,force)

global gSurfViewer;

if ieNotDefined('whichSurface')
  whichSurface = gSurfViewer.whichSurface;
end
if ieNotDefined('force')
  force = 0;
end

% what surface/rois are being displayed now
lastWhichSurface = gSurfViewer.whichSurface;

% get which surface to draw
if force || (whichSurface ~= lastWhichSurface)
  % set which surface and display
  gSurfViewer.whichSurface = whichSurface;
  % 1-4 are surfaces
  if whichSurface <= 4
    % if we are displaying the 3D anatomy, 
    % then switch to the surface view
    if lastWhichSurface > 4
      switchToSurface;
    else
      dispSurface;
    end
  % 5 is the volume
  elseif whichSurface == 5
    % switch to the volume view
    switchToVolume;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function switchToSurface

global gSurfViewer;
set(gSurfViewer.hSliders.v,'Visible','on');
set(gSurfViewer.hSliders.vText,'Visible','on');
set(gSurfViewer.hSliders.h,'Visible','on');
set(gSurfViewer.hSliders.hText,'Visible','on');
set(gSurfViewer.hSliders.h,'SliderStep',[15 45]./360);
set(gSurfViewer.hSliders.h,'Value',0);
set(gSurfViewer.hSliders.h,'Min',-180);
set(gSurfViewer.hSliders.h,'Max',180);
set(gSurfViewer.hSliders.h,'TooltipString','Rotate around z-axis');
set(gSurfViewer.hSliders.v,'Value',0);
set(gSurfViewer.hSliders.vText,'String',0);
set(gSurfViewer.hSliders.hText,'String',0);
dispSurface;
setViewAngle(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToVolume   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function switchToVolume

global gSurfViewer;
initSlice = 127;
set(gSurfViewer.hSliders.h,'Visible','on');
set(gSurfViewer.hSliders.hText,'Visible','on');
set(gSurfViewer.hSliders.v,'Visible','off');
set(gSurfViewer.hSliders.vText,'Visible','off');
set(gSurfViewer.hSliders.h,'SliderStep',[1 16]./256);
set(gSurfViewer.hSliders.h,'Value',initSlice);
set(gSurfViewer.hSliders.h,'Min',1);
set(gSurfViewer.hSliders.h,'Max',256);
set(gSurfViewer.hSliders.h,'TooltipString','Change viewing slice');
set(gSurfViewer.hSliders.hText,'String',num2str(initSlice));
dispVolume(3,initSlice);

%%%%%%%%%%%%%%%%%%%%%%%
%%   sliderHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function sliderHandler

global gSurfViewer;

% get slider position
hPos = round(get(gSurfViewer.hSliders.h,'Value'));
vPos = round(get(gSurfViewer.hSliders.v,'Value'));

% set the edit fields
set(gSurfViewer.hSliders.hText,'String',num2str(hPos));
set(gSurfViewer.hSliders.vText,'String',num2str(vPos));

if gSurfViewer.whichSurface <= 4
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   editHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function editHandler

global gSurfViewer;

hPos = str2num(get(gSurfViewer.hSliders.hText,'String'));
vPos = str2num(get(gSurfViewer.hSliders.vText,'String'));

% make it fit into -180:180
hPos = round(mod(hPos+180,360)-180);
vPos = round(mod(vPos+180,360)-180);

% set slider position
set(gSurfViewer.hSliders.h,'Value',hPos);
set(gSurfViewer.hSliders.v,'Value',vPos);

% set the edit fields
set(gSurfViewer.hSliders.hText,'String',num2str(hPos));
set(gSurfViewer.hSliders.vText,'String',num2str(vPos));

if gSurfViewer.whichSurface <= 4
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   setViewAngle   %%
%%%%%%%%%%%%%%%%%%%%%%
function setViewAngle(hPos,vPos)

% flip the sign to make rotations go in the "right" direction
hPos = -hPos;vPos = -vPos;
% somehow 90 and 180 are a problem for matlab
if abs(vPos) == 90,vPos = sign(vPos)*91;,end
if abs(hPos) == 90,hPos = sign(hPos)*91;,end
if abs(hPos) == 179,hPos = sign(hPos)*179;,end

% set the camera taret to center
camtarget([0 0 0]);

% set the size of the field of view in degrees
% i.e. 90 would be very wide and 1 would be ver
% narrow. 7 seems to fit the whole brain nicely
camva(7);

% set the view angle
view(hPos,vPos);

% change the camera position to avoid the volume
% flipping back and forth, another starnge matlab thing
if (vPos >= 90) || (vPos < -90)
  camup([0 0 -1]);
else
  camup([0 0 1]);
end

%%%%%%%%%%%%%%%%%%%%%
%%   dispSurface   %%
%%%%%%%%%%%%%%%%%%%%%
function dispSurface

global gSurfViewer;
figure(gSurfViewer.f);

% clear the axis
cla;

% get the vertexes/triangles and curvature
if gSurfViewer.whichSurface == 1
  vtcs = gSurfViewer.outerSurface.vtcs;
  tris = gSurfViewer.outerSurface.tris;
  c = gSurfViewer.curv;
elseif gSurfViewer.whichSurface == 2
  if isempty(gSurfViewer.outerCoords),return,end
  vtcs = gSurfViewer.outerCoords.vtcs;
  tris = gSurfViewer.outerCoords.tris;
  c = gSurfViewer.curv;
elseif gSurfViewer.whichSurface == 3
  if isempty(gSurfViewer.innerSurface),return,end
  vtcs = gSurfViewer.innerSurface.vtcs;
  tris = gSurfViewer.innerSurface.tris;
  c = gSurfViewer.curv;
elseif gSurfViewer.whichSurface == 4
  if isempty(gSurfViewer.innerCoords),return,end
  vtcs = gSurfViewer.innerCoords.vtcs;
  tris = gSurfViewer.innerCoords.tris;
  c = gSurfViewer.curv;
end  

% move vertices into center
vtcs(:,1) = vtcs(:,1)-mean(vtcs(:,1));
vtcs(:,2) = vtcs(:,2)-mean(vtcs(:,2));
vtcs(:,3) = vtcs(:,3)-mean(vtcs(:,3));

% not sure why, but this is necessary to set up
% the axis so that right is right...
imagesc(0);

% draw the surface
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', c, ...
      'facecolor', 'interp', ...
      'edgecolor', 'none');

% set axis stuff
axis off;axis equal;colormap(gray);axis tight;
camup('manual');
set(gca,'CLim',[-1.2 1.2]);

if gSurfViewer.whichSurface <= 4
  hPos = round(get(gSurfViewer.hSliders.h,'Value'));
  vPos = round(get(gSurfViewer.hSliders.v,'Value'));
  setViewAngle(hPos,vPos);
end

%%%%%%%%%%%%%%
% dispVolume
%%%%%%%%%%%%%%
function dispVolume(sliceIndex,slice)

global gSurfViewer;
figure(gSurfViewer.f);
cla;

% display a slice of the anatomy image
switch sliceIndex
  case {1}
   img = gSurfViewer.anat.data(slice,:);
  case {2}
   img = gSurfViewer.anat.data(:,slice,:);
  case {3}
   img = gSurfViewer.anat.data(:,:,slice);
end
imagesc(img);
colormap(gray);
axis image;
axis off;
hold on

if min(img(:)) ~= max(img(:))
  set(gca,'CLim',[min(img(:)) max(img(:))]);
end
% display patch and white matter/gray matter
innerNodes = gSurfViewer.innerSurface.vtcs;
outerNodes = gSurfViewer.outerSurface.vtcs;

% Plot the nodes for the gray/white matter surfaces
innerNodes = innerNodes( find( round(innerNodes(:,sliceIndex))==slice), : );
plot(innerNodes(:,1), innerNodes(:,2), '.', 'markersize', 1, 'Color', [0.7 0.7 0.7]);

outerNodes = outerNodes( find( round(outerNodes(:,sliceIndex))==slice), : );
plot(outerNodes(:,1), outerNodes(:,2), 'y.', 'markersize', 1,'Color',[0.7 0.7 0]);

if ~isempty(gSurfViewer.innerCoords)
  innerCoordNodes = gSurfViewer.innerCoords.vtcs;
  innerCoordNodes = innerCoordNodes( find( round(innerCoordNodes(:,sliceIndex))==slice), : );
  plot(innerCoordNodes(:,1), innerCoordNodes(:,2), 'w.', 'markersize', 1);
end
if ~isempty(gSurfViewer.outerCoords)
  outerCoordNodes = gSurfViewer.outerCoords.vtcs;
  outerCoordNodes = outerCoordNodes( find( round(outerCoordNodes(:,sliceIndex))==slice), : );
  plot(outerCoordNodes(:,1), outerCoordNodes(:,2), 'y.', 'markersize', 1);
end

view([0 90]);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for current path
%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isCurrentPath(pathname)

if length(pathname) && (pathname(end)==filesep)
  pathname = pathname(1:end-1);
end
tf = strcmp(pwd,pathname);

%%%%%%%%%%%%%%%%%%%%%%%
%%   switchAnatomy   %%
%%%%%%%%%%%%%%%%%%%%%%%
function switchAnatomy(params)

global gSurfViewer;

% check for find file
if strcmp(params.anatomy,'Find file')
  [filename, pathname] = uigetfile({'*.hdr','Nifti file (*.hdr)'});
  % update the control that displays the choices
  % first check to see if this is a new path
  global gParams;
  if ~isCurrentPath(pathname)
    filename = fullfile(pathname,filename);
  end
  whichControl = gSurfViewer.guiloc.filenames+6;
  currentChoices = get(gParams.ui.varentry{whichControl},'String');
  currentChoices = setdiff(currentChoices,'Find file');
  currentChoices = putOnTopOfList(filename,currentChoices);
  currentChoices{end+1} = 'Find file';
  set(gParams.ui.varentry{whichControl},'String',currentChoices)
  set(gParams.ui.varentry{whichControl},'Value',1)
  params.anatomy = filename;
end

% load the anatomy and view
disppercent(-inf,sprintf('(mrSurfViewer) Load %s',params.anatomy));
[gSurfViewer.anat.data gSurfViewer.anat.hdr] = cbiReadNifti(params.anatomy);
% switch to 3D anatomy view
global gParams
gSurfViewer.whichSurface = 5;
set(gParams.ui.varentry{1},'Value',gSurfViewer.whichSurface)
refreshFlatViewer([],1);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%%   switchFile   %%
%%%%%%%%%%%%%%%%%%%%
function switchFile(whichSurface,params)

global gSurfViewer;

% if the user wants to find a new file
addFilename = 0;
if strcmp(params.(whichSurface),'Find file')
  if strcmp(whichSurface,'curv')
    [filename, pathname] = uigetfile({'*.vff','VFF Curvature files (*.vff)'});
    whichControl = gSurfViewer.guiloc.filenames+5;
  else
    [filename, pathname] = uigetfile({'*.off','OFF Surface files (*.off)'});
    whichControl = gSurfViewer.guiloc.filenames+find(strcmp(whichSurface,{'outerSurface','outerCoords','innerSurface','innerCoords'}));
  end
  % see if file is in a different path
  if ~isCurrentPath(pathname)
    filename = fullfile(pathname,filename);
  end
  addFilename = 1;
else
  filename = params.(whichSurface);
end

% try to load it
disppercent(-inf,sprintf('(mrSurfViewer) Loading %s',filename));
if filename ~= 0
  if strcmp(whichSurface,'curv')
    file = myLoadCurvature(filename);
    whichControl = gSurfViewer.guiloc.filenames+5;
  else
    file = myLoadSurface(filename);
    whichControl = gSurfViewer.guiloc.filenames+find(strcmp(whichSurface,{'outerSurface','outerCoords','innerSurface','innerCoords'}));
  end
else
  file = [];
end
if ~isempty(file)
  if strcmp(whichSurface,'curv')
    gSurfViewer.(whichSurface)=file;
  else
    gSurfViewer.(whichSurface)=file;
    % set the correct one to display
    gSurfViewer.whichSurface = find(strcmp(whichSurface,{'outerSurface','outerCoords','innerSurface','innerCoords'}));
  end    
  % and change the ui control
  global gParams;
  set(gParams.ui.varentry{1},'Value',gSurfViewer.whichSurface)
  % add the filename to the control if necessary
  if addFilename
    currentChoices = get(gParams.ui.varentry{whichControl},'String');
    currentChoices = setdiff(currentChoices,'Find file');
    currentChoices = putOnTopOfList(filename,currentChoices);
    currentChoices{end+1} = 'Find file';
    set(gParams.ui.varentry{whichControl},'String',currentChoices)
    set(gParams.ui.varentry{whichControl},'Value',1)
  end
else
  global gParams;
  % switch back to first on list
  currentChoices = get(gParams.ui.varentry{whichControl},'String');
  set(gParams.ui.varentry{whichControl},'Value',1)
  if ~strcmp(whichSurface,'curv')
    gSurfViewer.(whichSurface)=myLoadSurface(currentChoices{1});
  else
    gSurfViewer.curv=myLoadCurvature(currentChoices{1});
  end
end
refreshFlatViewer([],1);

disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%
function surf = myLoadSurface(filename,filepath,onlyLoadHeader)

% default to loading data
if ieNotDefined('onlyLoadHeader'),onlyLoadHeader = 0;end

% put on path
if ~ieNotDefined('filepath')
  filename = fullfile(filepath,filename);
end

global gSurfViewer;
% load the surface
surf = loadSurfOFF(filename,onlyLoadHeader);
if isempty(surf),return,end

% check that it has the correct number of vertices
if ~isequal([gSurfViewer.outerSurface.Nvtcs gSurfViewer.outerSurface.Ntris],[surf.Nvtcs surf.Ntris])
  % dispaly warning, but only if mismatchWarning is set,
  % this way when we first load surfaces just for checking it
  % won't complain
  if gSurfViewer.mismatchWarning
    mrWarnDlg(sprintf('(mrSurfViewer) Surface %s does not match patch Nvtcs: %i vs %i, Ntris: %i vs %i',filename,gSurfViewer.outerSurface.Nvtcs,surf.Nvtcs,gSurfViewer.outerSurface.Ntris,surf.Ntris));
  end
  surf = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadCurvature   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function curv = myLoadCurvature(filename,filepath,onlyLoadHeader)

% default to loading data
if ieNotDefined('onlyLoadHeader'),onlyLoadHeader = 0;end

% put on path
if ~ieNotDefined('filepath')
  filename = fullfile(filepath,filename);
end

global gSurfViewer;
% load the curvature
[curv curvhdr] = loadVFF(filename,onlyLoadHeader);

if isempty(curvhdr)
  return
end
% % check that it has the correct number of vertices
if ~isequal(gSurfViewer.outerSurface.Nvtcs, curvhdr.size(3))
  % dispaly warning, but only if mismatchWarning is set,
  % this way when we first load surfaces just for checking it
  % won't complain
  if gSurfViewer.mismatchWarning
    mrWarnDlg(sprintf('(mrSurfViewer) Curvature file %s does not match patch Nvtcs: %i vs %i',filename,gSurfViewer.outerSurface.Nvtcs, curvhdr.size(3)));
  end
  curv = [];
  return
end

% return header only
if onlyLoadHeader
  curv = curvhdr;
else 
  curv = curv';
end
