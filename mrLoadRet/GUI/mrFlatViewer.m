% mrFlatViewer.m
%
%       $Id$	
%      usage: mrFlatViewer(flatname,outer,inner,curv,anat,viewNum)
%         by: justin gardner, originally based on surfViewer by eli merriam
%       date: 10/09/07
%    purpose: 
%
function retval = mrFlatViewer(flat,outer,inner,curv,anat,viewNum)

% check arguments
if ~any(nargin == [1 2 3 4 5 6])
  help mrFlatViewer
  return
end
if nargout == 1
  retval = [];
end

% if passed in a string check to see if
% it needs an extension
if isstr(flat)
  if isfile(sprintf('%s.off',stripext(flat)));
    flat = sprintf('%s.off',stripext(flat));
  end
end

% see how we are being called
if (nargin == 1) && isstr(flat) && ~isfile(flat)
  event = flat;
else
  event = 'init';
  % set defaults
  if ieNotDefined('outer'),outer = {};end
  if ieNotDefined('inner'),inner = {};end
  if ieNotDefined('curv'),curv = {};end
  if ieNotDefined('anat'),anat = {};end
  if ieNotDefined('viewNum'),viewNum = [];end
  % make everybody a cell array
  flat = cellArray(flat);
  outer = cellArray(outer);
  inner = cellArray(inner);
  curv = cellArray(curv);
  anat = cellArray(anat);
end

switch (event)
 case 'init'
  retval = initHandler(flat,outer,inner,curv,anat,viewNum);
 case {'vSlider','hSlider'}
  sliderHandler;
 case {'edit'}
  editHandler;
 otherwise
  disp(sprintf('(mrFlatViewer) Could not find flat file %s',event));
end

%%%%%%%%%%%%%%%%%%%%%%
%%   init handler   %%
%%%%%%%%%%%%%%%%%%%%%%
function retval = initHandler(flat,outer,inner,curv,anat,viewNum)

global gFlatViewer;
gFlatViewer = [];

disppercent(-inf,'(mrFlatView) Loading surfaces');
% get more flats
% load the flat
gFlatViewer.flat = loadSurfOFF(sprintf('%s.off',stripext(flat{1})));
if isempty(gFlatViewer.flat) || ~isfield(gFlatViewer.flat,'parentSurfaceName');
  disp(sprintf('(mrFlatViewer) %s is not a flat file',flat{1}));
  return
end
% look for flats with same parent
flatdir = dir('*.off');
for i = 1:length(flatdir)
  if ~isempty(strcmp(lower(flatdir(i).name),'flat')) || ~isempty(strcmp(lower(flatdir(i).name),'patch'))
    flatfile = loadSurfOFF(flatdir(i).name,1);
    if isfield(flatfile,'parentSurfaceName')
      if strcmp(flatfile.parentSurfaceName,gFlatViewer.flat.parentSurfaceName)
	if ~strcmp(flatdir(i).name,flat)
	  flat{end+1} = flatdir(i).name;
	end
      end
    end
  end
end

% remove any paths
gFlatViewer.flat.parentSurfaceName = getLastDir(gFlatViewer.flat.parentSurfaceName);

% load up the surfaces
if isempty(inner)
  % guess the names
  inner{1} = gFlatViewer.flat.parentSurfaceName;
end

% guess anything with the right stem
innerDir = dir(sprintf('%s*.off',strtok(stripext(inner{1}),'WM')));
for i = 1:length(innerDir)
  % don't choose anything we already have or with GM, flat or patch in the title
  if ~any(strcmp(innerDir(i).name,inner)) && isempty(strfind(innerDir(i).name,'GM')) && (isempty(strfind(lower(innerDir(i).name),'flat')) || ~isempty(strfind(lower(innerDir(i).name),'inflate'))) && isempty(strfind(lower(innerDir(i).name),'patch'))
    inner{end+1} = innerDir(i).name;
  end
end

% now try to find the first loadable one
for i = 1:length(inner)
  gFlatViewer.surfaces.inner = myLoadSurface(inner{i});
  if ~isempty(gFlatViewer.surfaces.inner),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.surfaces.inner)
  return
else
  inner = putOnTopOfList(inner{i},inner);
end
inner{end+1} = 'Find file';

% load the outer surface
if isempty(outer)
  % if we weren't passed in anything try to find them
  outer{1} = sprintf('%sGM.off',strtok(stripext(inner{1}),'WM'));
end
% guess anything with the right stem
outerDir = dir(sprintf('%s*.off',strtok(stripext(inner{1}),'WM')));
for i = 1:length(outerDir)
  % don't choose anything we already have or with WM, flat or patch in the title
  if ~any(strcmp(innerDir(i).name,inner)) && isempty(strfind(innerDir(i).name,'WM')) && (isempty(strfind(lower(innerDir(i).name),'flat')) || ~isempty(strfind(lower(innerDir(i).name),'inflate'))) && isempty(strfind(lower(innerDir(i).name),'patch'))
    outer{end+1} = outerDir(i).name;
  end
end
% now try to find the first loadable one
for i = 1:length(outer)
  gFlatViewer.surfaces.outer = myLoadSurface(outer{i});
  if ~isempty(gFlatViewer.surfaces.outer),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.surfaces.outer)
  return
else
  outer = putOnTopOfList(outer{i},outer);
end
outer{end+1} = 'Find file';

% load the curvature
if isempty(curv)
  curvGuess = sprintf('%s_Curv.vff',stripext(inner{1}));
  if isfile(curvGuess)
    curv{1} = curvGuess;
  end
end
% add any vff file
curvDir = dir('*.vff');
for i = 1:length(curvDir)
  if ~any(strcmp(curvDir(i).name,curv))
    curv{end+1} = curvDir(i).name;
  end
end

for i = 1:length(curv)
  gFlatViewer.curv = myLoadCurvature(curv{1});
  if ~isempty(gFlatViewer.curv),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.curv)
  return
else
  curv = putOnTopOfList(curv{i},curv);
end
disppercent(inf);

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
if isfile(anat{1})
  [gFlatViewer.anat.data gFlatViewer.anat.hdr] = cbiReadNifti(anat{1});
else
  gFlatViewer.anat = [];
end

% save the view
gFlatViewer.viewNum = viewNum;

% select the window
gFlatViewer.f = selectGraphWin;

% positions on figure
figLeft = 10;figBottom = 10;
sliderWidth = 20;sliderLength = 200;spacer = 10;
editWidth = 40;editHeight = 20;

% set up horizontal and vertical slider
gFlatViewer.hSliders.v = uicontrol('Style','slider','Position',[figLeft figBottom+sliderWidth sliderWidth sliderLength],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''vSlider'')','TooltipString','Rotate around y-axis');
gFlatViewer.hSliders.vText = uicontrol('Style','Edit','Position',[figLeft figBottom+sliderWidth+sliderLength+spacer editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0','HorizontalAlignment','Center');
gFlatViewer.hSliders.h = uicontrol('Style','slider','Position',[figLeft+sliderWidth figBottom sliderLength sliderWidth],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''hSlider'')','TooltipString','Rotate around z-axis');
gFlatViewer.hSliders.hText = uicontrol('Style','Edit','Position',[figLeft+sliderLength+3*spacer figBottom editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0');


% set they we are viewing white matter
gFlatViewer.whichSurface = 1;
gFlatViewer.patchColoring = 1;
gFlatViewer.displayROIs = 0;
% and display surface
dispSurface;
setViewAngle(0,0);

editable = 0;

% set up the parameters
paramsInfo = {};
if ~editable && (length(flat) == 1)
  paramsInfo{end+1} = {'flatFile',flat{1},'editable=0','The flat patch file'};
else
  paramsInfo{end+1} = {'flatFile',flat,'The flat patch file','callback',@switchFlat};
end
if ~editable && (length(outer) == 1)
  paramsInfo{end+1} = {'outer',outer{1},'editable=0','The outer (gray matter) file'};
else
  paramsInfo{end+1} = {'outer',outer,'The outer (gray matter) file','callback',@switchFile,'callbackArg=outer'};
end
if ~editable && (length(inner) == 1)
  paramsInfo{end+1} = {'inner',inner{1},'editable=0','The inner (white matter) file'};
else
  paramsInfo{end+1} = {'inner',inner,'The inner (white matter) file','callback',@switchFile,'callbackArg=inner'};
end
if ~editable && (length(curv) == 1)
  paramsInfo{end+1} = {'curv',curv{1},'editable=0','The curvature file'};
else
  paramsInfo{end+1} = {'curv',curv,'The curvature file','callback',@switchFile,'callbackArg=curv'};
end
if ~editable && (length(anat) == 1)
  paramsInfo{end+1} = {'anatomy',anat{1},'editable=0','The 3D anatomy file'};
else
  paramsInfo{end+1} = {'anatomy',anat,'The 3D anatomy file','callback',@switchAnatomy};
end
% Now give choice of viewing gray or white
gFlatViewer.whichSurfaceTypes = {'Outer (Gray matter) surface','Inner (White matter) surface','3D Anatomy','Patch'};
paramsInfo{end+1} = {'whichSurface',gFlatViewer.whichSurfaceTypes,'callback',@whichSurfaceCallback,'Choose which surface to view the patch on'};
gFlatViewer.patchColoringTypes = {'Uniform','Rostral in red','Right in red','Dorsal in red','Positive curvature in red','Negative curvature in red','Compressed areas in red','Stretched areas in red','High outer areal distortion in red','High inner areal distortion in red'};
if ~isempty(gFlatViewer.viewNum)
  gFlatViewer.patchColoringTypes{end+1} = 'Current overlay';
end
gFlatViewer.patchColoringTypes{end+1} = 'None';
paramsInfo{end+1} = {'patchColoring',gFlatViewer.patchColoringTypes,'Choose how to color the patch','callback',@patchColoringCallback};
if ~isempty(gFlatViewer.viewNum)
  paramsInfo{end+1} = {'displayROIs',0,'type=checkbox','Display the ROIs','callback',@whichSurfaceCallback};
end

% put up dialog
params = mrParamsDialog(paramsInfo,'View flat patch location on surface');

if isempty(params)
  close(gFlatViewer.f);
  retval = [];
else
  retval = params;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   patchColoringCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function patchColoringCallback(params)

global gFlatViewer
gFlatViewer.patchColoring = find(strcmp(params.patchColoring,gFlatViewer.patchColoringTypes));
if gFlatViewer.whichSurface == 3
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  dispVolume(3,hPos);
else
  dispSurface;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   whichSurfaceCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichSurfaceCallback(params)

global gFlatViewer;

% set the roi drawing
if isfield(params,'displayROIs')
  gFlatViewer.displayROIs = params.displayROIs;
else
  params.displayROIs = gFlatViewer.displayROIs;
end

% get which surface to draw
lastWhichSurface = gFlatViewer.whichSurface;

refreshFlatViewer(find(strcmp(params.whichSurface,gFlatViewer.whichSurfaceTypes)),params.displayROIs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   refreshFlatViewer   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshFlatViewer(whichSurface,displayROIs,force)

global gFlatViewer;

if ieNotDefined('whichSurface')
  whichSurface = gFlatViewer.whichSurface;
end
if ieNotDefined('displayROIs')
  displayROIs = gFlatViewer.displayROIs;
end
if ieNotDefined('force')
  force = 0;
end
% what surface/rois are being displayed now
lastWhichSurface = gFlatViewer.whichSurface;
lastDisplayROIs = gFlatViewer.displayROIs;

% get which surface to draw
if force || (whichSurface ~= lastWhichSurface) || (lastDisplayROIs ~= displayROIs)
  % set which surface and display
  gFlatViewer.whichSurface = whichSurface;
  % 1,2 are surfaces
  if whichSurface <= 2
    % if we are displaying the 3D anatomy, 
    % then switch to the surface view
    if lastWhichSurface > 2
      switchToSurface;
    else
      dispSurface;
    end
  % 3 is the volume
  elseif whichSurface == 3
    % switch to the volume view
    switchToVolume;
  else
    % switch to the volume view
    switchToFlat;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function switchToSurface

global gFlatViewer;
set(gFlatViewer.hSliders.v,'Visible','on');
set(gFlatViewer.hSliders.vText,'Visible','on');
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.h,'SliderStep',[15 45]./360);
set(gFlatViewer.hSliders.h,'Value',0);
set(gFlatViewer.hSliders.h,'Min',-180);
set(gFlatViewer.hSliders.h,'Max',180);
set(gFlatViewer.hSliders.h,'TooltipString','Rotate around z-axis');
set(gFlatViewer.hSliders.v,'Value',0);
set(gFlatViewer.hSliders.vText,'String',0);
set(gFlatViewer.hSliders.hText,'String',0);
dispSurface;
setViewAngle(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToVolume   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function switchToVolume

global gFlatViewer;
initSlice = 127;
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.v,'Visible','off');
set(gFlatViewer.hSliders.vText,'Visible','off');
set(gFlatViewer.hSliders.h,'SliderStep',[1 16]./256);
set(gFlatViewer.hSliders.h,'Value',initSlice);
set(gFlatViewer.hSliders.h,'Min',1);
set(gFlatViewer.hSliders.h,'Max',256);
set(gFlatViewer.hSliders.h,'TooltipString','Change viewing slice');
set(gFlatViewer.hSliders.hText,'String',num2str(initSlice));
dispVolume(3,initSlice);

%%%%%%%%%%%%%%%%%%%%%%
%%   switchToFlat   %%
%%%%%%%%%%%%%%%%%%%%%%
function switchToFlat

global gFlatViewer;
initSlice = 127;
set(gFlatViewer.hSliders.v,'Visible','off');
set(gFlatViewer.hSliders.vText,'Visible','off');
set(gFlatViewer.hSliders.h,'Visible','off');
set(gFlatViewer.hSliders.hText,'Visible','off');
dispSurface;

%%%%%%%%%%%%%%%%%%%%%%%
%%   sliderHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function sliderHandler

global gFlatViewer;

% get slider position
hPos = round(get(gFlatViewer.hSliders.h,'Value'));
vPos = round(get(gFlatViewer.hSliders.v,'Value'));

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if gFlatViewer.whichSurface <= 2
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   editHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function editHandler

global gFlatViewer;

hPos = str2num(get(gFlatViewer.hSliders.hText,'String'));
vPos = str2num(get(gFlatViewer.hSliders.vText,'String'));

% make it fit into -180:180
hPos = round(mod(hPos+180,360)-180);
vPos = round(mod(vPos+180,360)-180);

% set slider position
set(gFlatViewer.hSliders.h,'Value',hPos);
set(gFlatViewer.hSliders.v,'Value',vPos);

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if gFlatViewer.whichSurface <= 2
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

global gFlatViewer;
figure(gFlatViewer.f);

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

% get the vertexes/triangles and curvature
if gFlatViewer.whichSurface == 1
  vtcs = gFlatViewer.surfaces.outer.vtcs;
  tris = gFlatViewer.surfaces.outer.tris;
  c = gFlatViewer.curv;
elseif gFlatViewer.whichSurface == 2
  vtcs = gFlatViewer.surfaces.inner.vtcs;
  tris = gFlatViewer.surfaces.inner.tris;
  c = gFlatViewer.curv;
else
  vtcs = gFlatViewer.flat.vtcs;
  tris = gFlatViewer.flat.tris;
  c = gFlatViewer.curv(patchVtcs);
%  c = (c-min(c))./((max(c)-min(c)))>0.5;
  view([0 90]);
end  

% clear the axis
cla;

% not sure why, but this is necessary to set up
% the axis so that right is right...
imagesc(0);

% get the colors that we want to show for that patch
[co alpha] = getPatchColoring;

% now set the overlay
if gFlatViewer.whichSurface <= 2
  overlay = NaN(length(c),3);
  overlay(patchVtcs,:) = co;
else
  overlay(:,:) = co;
end

% get the roi overlay
if isfield(gFlatViewer,'viewNum') && gFlatViewer.displayROIs
  % recompute roiOverlay
  if ~isfield(gFlatViewer,'roiOverlays') || ...
	length(gFlatViewer.roiOverlays) < gFlatViewer.whichSurface || ...
	isempty(gFlatViewer.roiOverlays{gFlatViewer.whichSurface})
    % get the vertices for which to calculate the roi overlay
    if gFlatViewer.whichSurface <= 2
      baseCoords = round(vtcs);
    else
      baseCoords = round(gFlatViewer.surfaces.inner.vtcs(patchVtcs,:));
    end
    % and compute them
    gFlatViewer.roiOverlays{gFlatViewer.whichSurface} = computeROIOverlay(baseCoords);
  end
  % get the overlay from the global
  roiOverlay = gFlatViewer.roiOverlays{gFlatViewer.whichSurface};
else
  roiOverlay = [];
end

% move vertices into center
vtcs(:,1) = vtcs(:,1)-mean(vtcs(:,1));
vtcs(:,2) = vtcs(:,2)-mean(vtcs(:,2));
vtcs(:,3) = vtcs(:,3)-mean(vtcs(:,3));

% draw the surface and the overlay
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', c, ...
      'facecolor', 'interp', ...
      'edgecolor', 'none');
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', overlay, ...
      'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',alpha);
if ~isempty(roiOverlay)
  patch('vertices', vtcs, 'faces', tris, ...
	'FaceVertexCData', roiOverlay, ...
	'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',.4);
end

% set axis stuff
axis off;axis equal;colormap(gray);axis tight;
camup('manual');
set(gca,'CLim',[-1.2 1.2]);

if gFlatViewer.whichSurface <= 2
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  vPos = round(get(gFlatViewer.hSliders.v,'Value'));
  setViewAngle(hPos,vPos);
end
%%%%%%%%%%%%%%
% dispVolume
%%%%%%%%%%%%%%
function dispVolume(sliceIndex,slice)

global gFlatViewer;
figure(gFlatViewer.f);
cla;

% display a slice of the anatomy image
switch sliceIndex
  case {1}
   img = gFlatViewer.anat.data(slice,:);
  case {2}
   img = gFlatViewer.anat.data(:,slice,:);
  case {3}
   img = gFlatViewer.anat.data(:,:,slice);
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
whichInx = gFlatViewer.flat.patch2parent(:,2);
wmPatchNodes = gFlatViewer.surfaces.inner.vtcs(whichInx,:);
gmPatchNodes = gFlatViewer.surfaces.outer.vtcs(whichInx,:);

% get full white matter/gray matter nodes
wmNodes = gFlatViewer.surfaces.inner.vtcs;
gmNodes = gFlatViewer.surfaces.outer.vtcs;

% Plot the nodes for the gray/white matter surfaces
wmNodes = wmNodes( find( round(wmNodes(:,sliceIndex))==slice), : );
plot(wmNodes(:,1), wmNodes(:,2), 'w.', 'markersize', 1);

gmNodes = gmNodes( find( round(gmNodes(:,sliceIndex))==slice), : );
plot(gmNodes(:,1), gmNodes(:,2), 'y.', 'markersize', 1);


% plot the patch nodes, displaying both deep and superficial surfaces
co = getPatchColoring;
if ~isnan(co(1))
  wmco = co(find( round(wmPatchNodes(:,sliceIndex))==slice),:);
  % make into magenta vs blue
  i = wmco(:,1);
  wmco(:,1) = i;
  wmco(:,2) = 0;
  wmco(:,3) = max(i,1-i);
else
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),1)=1;
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),2)=1;
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),3)=1;
end

wmPatchNodes = wmPatchNodes( find( round(wmPatchNodes(:,sliceIndex))==slice), : );

if any(gFlatViewer.patchColoring == [1 length(gFlatViewer.patchColoringTypes)])
  % draw all the points in the same color (if there are any)
  if ~ieNotDefined('wmco')
    plot(wmPatchNodes(:,1), wmPatchNodes(:,2), '.', 'markersize', 1,'Color',wmco(1,:));
  end
  % otherwise each pixel has to be set
else
  for i = 1:length(wmPatchNodes(:,1))
    plot(wmPatchNodes(i,1), wmPatchNodes(i,2), '.', 'markersize', 1,'Color',wmco(i,:)');
  end
end

if ~isnan(co(1))
  gmco = co(find( round(gmPatchNodes(:,sliceIndex))==slice),:);
else
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),1)=1;
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),2)=1;
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),3)=0;
end
gmPatchNodes = gmPatchNodes( find( round(gmPatchNodes(:,sliceIndex))==slice), : );
% uniform patch coloring
if any(gFlatViewer.patchColoring == [1 length(gFlatViewer.patchColoringTypes)])
  % draw all the points in the same color (if there are any)
  if ~ieNotDefined('gmco')
    plot(gmPatchNodes(:,1), gmPatchNodes(:,2), '.', 'markersize', 1,'Color',gmco(1,:));
  end
  % otherwise each pixel has to be set
else
  for i = 1:length(gmPatchNodes(:,1))
    plot(gmPatchNodes(i,1), gmPatchNodes(i,2), '.', 'markersize', 1,'Color',gmco(i,:));
  end
end

view([0 90]);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getPatchColoring   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [co alpha] = getPatchColoring

global gFlatViewer;

alpha = 0.6;

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

% one is uniform, 2-4 are red/blue
switch gFlatViewer.patchColoring
 % uniform
 case 1
  % make everybody red
  co = ones(1,gFlatViewer.flat.Nvtcs);
  % anatomical directions
 case {2,3,4}
  co = gFlatViewer.surfaces.outer.vtcs(patchVtcs,gFlatViewer.patchColoring-1)';
  co = (co-min(co))./(max(co)-min(co));
  % curvature
 case {5,6}
  % get curvature
  curv = gFlatViewer.curv(patchVtcs)';
  curv = (curv-min(curv))./(max(curv)-min(curv));
  if gFlatViewer.patchColoring == 6
    curv = 1-curv;
  end
  % flatten distribution
  co = flattenDistribution(curv);
 % Areal distortion
 case {7,8,9,10}
  tris = gFlatViewer.flat.tris;
  % get area in volume
  if gFlatViewer.patchColoring == 9
    % use outer surface
    volumeArea = getTriangleArea(gFlatViewer.surfaces.outer.vtcs(patchVtcs,:),tris);
  else
    % use inner surface
    volumeArea = getTriangleArea(gFlatViewer.surfaces.inner.vtcs(patchVtcs,:),tris);
  end
  patchArea = getTriangleArea(gFlatViewer.flat.vtcs,tris);
  % get the distortion as the ratio of the area in the
  % patch to the volume
  trisDistortion = patchArea./volumeArea;
  % now convert this into a color for each vertex
  distortion = ones(1,length(patchVtcs));
  % note that this is a shortcut, it just
  % sets each vertex to one value of the distortion
  % even though the vertex may belong to many triangles
  distortion(tris(:,1)) = trisDistortion;
  distortion(tris(:,2)) = trisDistortion;
  distortion(tris(:,3)) = trisDistortion;
  % normalize
  if gFlatViewer.patchColoring < 9
    co = (distortion-min(distortion))./(max(distortion)-min(distortion));
    co = flattenDistribution(co);
    % invert colors
    distortion = log10(distortion);
    if gFlatViewer.patchColoring == 8
      co = 1-co;
      title(sprintf('Stretch: max=%0.2fx median=%0.2fx mean=%0.2fx',10^max(distortion),10^median(distortion(distortion>0)),10^mean(distortion(distortion>0))));
    else
      title(sprintf('Compression: max=%0.2fx median=%0.2fx mean=%0.2fx',10^abs(min(distortion)),10^abs(median(distortion(distortion<0))),10^abs(mean(distortion(distortion<0)))));
    end
  else
    % set the colors to the absolute value
    % of the log of the distortion. This scales
    % so that doubling or halving the area from
    % the volume to the patch give you the same
    % number. Anything above 10x distortion is 1
    distortion = abs(log10(distortion));
    co = distortion;
    co(co>1) = 1;
    % print out some statistics
    distortion = 10.^distortion;
    title(sprintf('Distortion: max=%0.2fx median=%0.2fx mean=%0.2fx',max(distortion),median(distortion),mean(distortion)));
  end
 %current overlay
 case {11}
  baseXform = gFlatViewer.anat.hdr.sform44;
  whichInx = gFlatViewer.flat.patch2parent(:,2);
  % get the coordinates from the right surface
  if gFlatViewer.whichSurface == 1
    baseCoords = gFlatViewer.surfaces.outer.vtcs(whichInx,:);
  else
    baseCoords = gFlatViewer.surfaces.inner.vtcs(whichInx,:);    
  end
  % swap x/y
  temp = baseCoords(:,2);
  baseCoords(:,2) = baseCoords(:,1);
  baseCoords(:,1) = temp;
  baseCoords(:,4) = 1;
  baseCoords = baseCoords';
  baseDims = [size(baseCoords,2) 1];
  % get the view
  v = viewGet([],'view',gFlatViewer.viewNum);
  if ~isempty(viewGet(v,'currentOverlay'))
    overlay = computeOverlay(v,baseXform,baseCoords,baseDims);
    overlay.RGB(overlay.alphaMap==0) = nan;
    co = squeeze(overlay.RGB);
  else
    co = zeros(size(baseCoords,2),3);
    co(:) = nan;
  end
  alpha = 1;
  return
 % no coloring
 case {length(gFlatViewer.patchColoringTypes)}
  co = nan;
  return
end

% intermediate values turn to gray, so avoid them
co((co>0.3)&(co<0.5)) = 0.3;
co((co>0.5)&(co<0.7)) = 0.7;

% make into RGB
co(2:3,:) = [1-co;1-co];
co = co';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getTriangleArea   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function area = getTriangleArea(vtcs,tris)

% get length of a each side of triangles
a = sqrt(sum((vtcs(tris(:,1),:)-vtcs(tris(:,2),:))'.^2));
b = sqrt(sum((vtcs(tris(:,2),:)-vtcs(tris(:,3),:))'.^2));
c = sqrt(sum((vtcs(tris(:,3),:)-vtcs(tris(:,1),:))'.^2));

% get semiperimeter (i.e. 1/2 perimeter)
p = (a+b+c)/2;

% use Heron's formula for size of triangle given side lengths
area = sqrt(p.*(p-a).*(p-b).*(p-c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   flattenDistribution   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function co = flattenDistribution(co,nbins)

if ~exist('nbins','var'),nbins = 10;end

% sort the values
[co sortIndex] = sort(co);
% get binsize
binSize = floor(length(co)/nbins);
% add even number of values back into each bin
for i = 1:nbins
  co(sortIndex((i-1)*binSize+1:min(length(co),i*binSize))) = i/nbins;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   computeROIOverlay   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiOverlay = computeROIOverlay(baseCoords);

global gFlatViewer;
roiOverlay = [];
disppercent(-inf,'(mrFlatViewer) Computing ROI Overlay');

% get view information
v = viewGet([],'view',gFlatViewer.viewNum);
numROIs = viewGet(v,'numROIs');
baseXform = viewGet(v,'baseXform');
baseVoxelSize = [1 1 1];
% swap x/y
tmp = baseCoords(:,2);
baseCoords(:,2) = baseCoords(:,1);
baseCoords(:,1) = tmp;
  
% init the overlay
roiOverlay = zeros(size(baseCoords,1),3);
roiOverlay(:) = nan;

% deal with selected ROI color
selectedROI = viewGet(v,'currentroi');
selectedROIColor = mrGetPref('selectedROIColor');

% get which ROIs to do
showROIs = viewGet(v,'showROIs');
if strcmp(showROIs,'none')
  return
end
if strfind(showROIs,'selected')
  rois = selectedROI;
else
  rois = 1:numROIs;
end

for roinum = rois
  % get ROI info
  roiCoords = viewGet(v,'roiCoords',roinum);
  roiXform = viewGet(v,'roiXform',roinum);
  roiVoxelSize = viewGet(v,'roiVoxelSize',roinum);
  if roinum ~= selectedROI
    roiColorRGB = viewGet(v,'roiColorRGB',roinum);
  else
    if strcmp(selectedROIColor,'none')
      roiColorRGB = viewGet(v,'roiColorRGB',roinum);
    else
      roiColorRGB = color2RGB(selectedROIColor);
    end
  end
  % get the base coord that match the roi
  roiBaseCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));
  roiBaseCoords = roiBaseCoords(1:3,:)';
  roiVertices = find(ismember(baseCoords,roiBaseCoords,'rows'));
  % and set them to the roi color
  roiOverlay(roiVertices,1) = roiColorRGB(1);
  roiOverlay(roiVertices,2) = roiColorRGB(2);
  roiOverlay(roiVertices,3) = roiColorRGB(3);
  disppercent(roinum/numROIs);
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%%   switchAnatomy   %%
%%%%%%%%%%%%%%%%%%%%%%%
function switchAnatomy(params)

global gFlatViewer;

% load the anatomy and view
disppercent(-inf,sprintf('(mrFlatViewer) Load %s',params.anatomy));
[gFlatViewer.anat.data gFlatViewer.anat.hdr] = cbiReadNifti(params.anatomy);
% switch to 3D anatomy view
global gParams
gFlatViewer.whichSurface = 3;
set(gParams.ui.varentry{6},'Value',gFlatViewer.whichSurface)
refreshFlatViewer([],[],1);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%%   switchFlat   %%
%%%%%%%%%%%%%%%%%%%%
function switchFlat(params)

global gFlatViewer;

% load the anatomy and view
disppercent(-inf,sprintf('(mrFlatViewer) Load %s',params.flatFile));
gFlatViewer.flat = loadSurfOFF(params.flatFile);
% switch to flat view
global gParams
refreshFlatViewer([],[],1);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%%   switchFile   %%
%%%%%%%%%%%%%%%%%%%%
function switchFile(whichSurface,params)

global gFlatViewer;

% if the user wants to find a new file
addFilename = 0;
if strcmp(params.(whichSurface),'Find file')
  if strcmp(whichSurface,'curv')
    [filename, pathname] = uigetfile({'*.vff','VFF Curvature files (*.vff)'});
    whichControl = 1+find(strcmp(whichSurface,{'outer','inner'}));
  else
    [filename, pathname] = uigetfile({'*.off','OFF Surface files (*.off)'});
    whichControl = 4;
  end
  % uigetfile seems to return a path with filesep at end
  if length(pathname) && (pathname(end)==filesep)
    pathname = pathname(1:end-1);
  end
  % see if file is in a different path
  if ~strcmp(pwd,pathname)
    filename = fullfile(pathname,filename);
  end
  addFilename = 1;
else
  filename = params.(whichSurface);
end

% try to load it
disppercent(-inf,sprintf('(mrFlatViewer) Loading %s',filename));
if filename ~= 0
  if strcmp(whichSurface,'curv')
    file = myLoadCurvature(filename);
    whichControl = 4;
  else
    file = myLoadSurface(filename);
    whichControl = 1+find(strcmp(whichSurface,{'outer','inner'}));
  end
else
  file = [];
end
if ~isempty(file)
  if strcmp(whichSurface,'curv')
    gFlatViewer.(whichSurface)=file;
  else
    gFlatViewer.surfaces.(whichSurface)=file;
    % set the correct one to display
    gFlatViewer.whichSurface = find(strcmp(whichSurface,{'outer','inner'}));
  end    
  % and change the ui control
  global gParams;
  set(gParams.ui.varentry{6},'Value',gFlatViewer.whichSurface)
  % add the filename to the control if necessary
  if addFilename
    currentChoices = get(gParams.ui.varentry{gFlatViewer.whichSurface+1},'String');
    currentChoices = setdiff(currentChoices,'Find file');
    currentChoices = putOnTopOfList(filename,currentChoices);
    currentChoices{end+1} = 'Find file';
    set(gParams.ui.varentry{gFlatViewer.whichSurface+1},'String',currentChoices)
    set(gParams.ui.varentry{gFlatViewer.whichSurface+1},'Value',1)
  end
else
  global gParams;
  % switch back to first on list
  currentChoices = get(gParams.ui.varentry{whichControl},'String');
  set(gParams.ui.varentry{whichControl},'Value',1)
  if ~strcmp(whichSurface,'curv')
    gFlatViewer.surfaces.(whichSurface)=myLoadSurface(currentChoices{1});
  else
    gFlatViewer.curv=myLoadCurvature(currentChoices{1});
  end
end
refreshFlatViewer([],[],1);

disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%
function surf = myLoadSurface(filename)

global gFlatViewer;
% load the surface
surf = loadSurfOFF(filename);
if isempty(surf),return,end

% check that it has the correct number of vertices
if ~isequal(gFlatViewer.flat.nParent(1:2),[surf.Nvtcs surf.Ntris]')
  disp(sprintf('(mrFlatViewer) Surface %s does not match patch Nvtcs: %i vs %i, Ntris: %i vs %i',filename,gFlatViewer.flat.nParent(1),surf.Nvtcs,gFlatViewer.flat.nParent(2),surf.Ntris));
  surf = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadCurvature   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function curv = myLoadCurvature(filename)

global gFlatViewer;
% load the curvature
curv = loadVFF(filename)';
if isempty(curv),return,end

% check that it has the correct number of vertices
if ~isequal(gFlatViewer.flat.nParent(1),size(curv,1))
  disp(sprintf('(mrFlatViewer) Curvature file %s does not match patch Nvtcs: %i vs %i',filename,gFlatViewer.flat.nParent(1),size(curv,1)));
  curv = [];
  return
end

