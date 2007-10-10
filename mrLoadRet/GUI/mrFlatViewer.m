% mrFlatViewer.m
%
%      usage: mrFlatViewer(flatname)
%         by: modified by jg from surfViewer by eli merriam
%       date: 10/09/07
%    purpose: 
%
function retval = mrFlatViewer(flat,GM,WM,curv)

% check arguments
if ~any(nargin == [1 4])
  help mrFlatViewer
  return
end

% see how we are being called
if (nargin == 1) && isstr(flat) && ~isfile(flat)
  event = flat;
else
  event = 'init';
  % set defaults
  if ieNotDefined('GM'),GM = {};end
  if ieNotDefined('WM'),WM = {};end
  if ieNotDefined('curv'),curv = {};end
  % make everybody a cell array
  flat = cellArray(flat);
  GM = cellArray(GM);
  WM = cellArray(WM);
  curv = cellArray(curv);
end

switch (event)
 case 'init'
  initHandler(flat,GM,WM,curv);
 case {'vSlider','hSlider'}
  sliderHandler;
 case {'edit'}
  editHandler;
end

%%%%%%%%%%%%%%%%%%%%%%
%%   init handler   %%
%%%%%%%%%%%%%%%%%%%%%%
function initHandler(flat,GM,WM,curv)

global gFlatViewer;

disppercent(-inf,'(mrFlatView) Loading surfaces');
% load the flat
gFlatViewer.flat = loadSurfOFF(flat{1});
if isempty(gFlatViewer.flat) || ~isfield(gFlatViewer.flat,'parentSurfaceName');
  disp(sprintf('(mrFlatViewer) %s is not a flat file',flat{1}));
  return
end

% load up the surfaces
if isempty(WM)
  WM{1} = gFlatViewer.flat.parentSurfaceName;
end
gFlatViewer.surfaces.WM = loadSurfOFF(WM{1});

% load the gray
if isempty(GM)
  GM{1} = sprintf('%sGM.off',strtok(stripext(WM{1}),'WM'));
end
gFlatViewer.surfaces.GM = loadSurfOFF(GM{1});

% load the curvature
if isempty(curv)
  curv{1} = sprintf('%s_Curv.vff',stripext(WM{1}));
end
gFlatViewer.curv = loadVFF(curv{1})';
disppercent(inf);

gFlatViewer.f = selectGraphWin;
%gFlatViewer.f = gcf;

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
gFlatViewer.whichSurface= 1;

% and display surface
dispSurface;

axis off;axis equal;colormap(gray);axis tight;
camup('manual');
set(gca,'CLim',[-1.2 1.2]);
setViewAngle(0,0);
editable = 0;

% set up the parameters
paramsInfo = {};
if ~editable && (length(flat) == 1)
  paramsInfo{end+1} = {'flatFile',flat{1},'editable=0','The flat patch file'};
else
  paramsInfo{end+1} = {'flatFile',flat,'The flat patch file'};
end
if ~editable && (length(GM) == 1)
  paramsInfo{end+1} = {'GMSurface',GM{1},'editable=0','The gray matter file'};
else
  paramsInfo{end+1} = {'GMSurface',GM,'The gray matter file'};
end
if ~editable && (length(WM) == 1)
  paramsInfo{end+1} = {'WMSurface',WM{1},'editable=0','The white matter file'};
else
  paramsInfo{end+1} = {'WMSurface',WM,'The white matter file'};
end
if ~editable && (length(curv) == 1)
  paramsInfo{end+1} = {'curvature',curv{1},'editable=0','The curvature file'};
else
  paramsInfo{end+1} = {'curvature',curv,'The curvature file'};
end
% Now give choice of viewing gray or white
paramsInfo{end+1} = {'whichSurface',{'Gray matter','White matter'},'callback',@whichSurfaceCallback,'Choose which surface to view the patch on'};
mrParamsDialog(paramsInfo,'View flat patch location on surface');


close(gFlatViewer.f);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   whichSurfaceCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichSurfaceCallback(params)

global gFlatViewer;

% get which surface to draw
whichSurface = find(strcmp(params.whichSurface,{'Gray matter','White matter'}));
if whichSurface ~= gFlatViewer.whichSurface
  gFlatViewer.whichSurface = whichSurface;
  dispSurface;
end


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

setViewAngle(hPos,vPos);

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

setViewAngle(hPos,vPos);

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
% get the vertexes/triangles and curvature
if gFlatViewer.whichSurface == 1
  vtcs = gFlatViewer.surfaces.WM.vtcs;
  tris = gFlatViewer.surfaces.WM.tris;
else
  vtcs = gFlatViewer.surfaces.GM.vtcs;
  tris = gFlatViewer.surfaces.GM.tris;
end
  
c = gFlatViewer.curv;

% make vertices into center
vtcs(:,1) = vtcs(:,1)-mean(vtcs(:,1));
vtcs(:,2) = vtcs(:,2)-mean(vtcs(:,2));
vtcs(:,3) = vtcs(:,3)-mean(vtcs(:,3));

% clear the axis
cla;

% and draw the surface
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', c, ...
      'facecolor', 'interp', ...
      'edgecolor', 'none');

% and draw the patch in red
overlay = NaN(length(c),3);
overlay(gFlatViewer.flat.patch2parent(:,2),1) = 1;
overlay(gFlatViewer.flat.patch2parent(:,2),2:3) = 0;
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', overlay, 'FaceVertexAlphaData', overlay(:,1)*.1, ...
      'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',.6);
