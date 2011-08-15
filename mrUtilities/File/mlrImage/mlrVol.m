% mlrVol.m
%
%        $Id:$ 
%      usage: mlrVol(filename)
%         by: justin gardner
%       date: 08/15/11
%    purpose: Volume viewer to display volume multi-dimensional data (usually 3D volumes or 4D epi sequences)
%
function retval = mlrVol(filename,varargin)

% check arguments
if nargin < 1
  help mlrVol
  return
end

if isstr(filename)
  % will create two globals:
  % gVols contains each volumes data and header that we are currently displaying
  % gSystem which scontains system parameters
  initSystem(varargin);

  % load the volume
  loadVolume(filename);

  % display the volume
  refreshDisplay;
else
  % command argument
  switch (filename)
   case 1
    mouseDownHandler(varargin{1});
   case 2
    mouseMoveHandler(varargin{1});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mouseMoveHadnler    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseMoveHandler(sysNum)

global gSystem;

coord = getMouseCoord;
if isempty(coord)
  set(gcf,'Pointer','arrow');
else
  set(gcf,'Pointer','fullcrosshair');
  set(gSystem.hTextboxX,'String',coord(1));
  set(gSystem.hTextboxY,'String',coord(2));
  set(gSystem.hTextboxZ,'String',coord(3));
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mouseDownHandler    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(sysNum)

global gVols;
global gSystem;

% FIX: get current vol, this should be a function of which axis was clicked
vol = gVols(end);

% get the mouse pos
coord = getMouseCoord;
if isempty(coord),return,end
% and set the volume coordinate
setVolCoord(gSystem.n,coord);

% refresh the display
refreshDisplay;

%%%%%%%%%%%%%%%%%%%%%%%
%    getMouseCoord    %
%%%%%%%%%%%%%%%%%%%%%%%
function coord = getMouseCoord

coord = [];viewNum = [];

global gVols;
global gSystem;

vol = gVols(end);

% figure out which axis we are on
pointerLoc = get(gcf,'CurrentPoint');
pos = get(gcf,'Position');
pos = pointerLoc./pos(3:4);
a = subplot(1,3,ceil(pos(1)*3));
viewNum = find(gSystem.a == a);

% get pointer loc
pointerLoc = get(a,'CurrentPoint');
pointerLoc = round(pointerLoc(1,2:-1:1));
% check image boundaries
if (pointerLoc(1) < 1) | (pointerLoc(1) > vol.h.dim(vol.xDim(viewNum))),return,end
if (pointerLoc(2) < 1) | (pointerLoc(2) > vol.h.dim(vol.yDim(viewNum))),return,end

% now get the point on the image
coord(vol.xDim(viewNum)) = pointerLoc(1);
coord(vol.yDim(viewNum)) = pointerLoc(2);
coord(vol.viewDim(viewNum)) = vol.coord(vol.viewDim(viewNum));


%%%%%%%%%%%%%%%%%%%%%%%%
%    refreshDisplay    %
%%%%%%%%%%%%%%%%%%%%%%%%
function refreshDisplay

global gVols;
global gSystem;

f = gcf;a = gca;
for iVol = 1:gSystem.n
  % display the coordinate
  dispVolume(gVols(iVol))
  % set that we have displayed this coordinate
  gVols(iVol).curCoord = gVols(iVol).coord;
end

drawnow;
figure(f);%axes(a);
%%%%%%%%%%%%%%%%%%%%
%    dispVolume    %
%%%%%%%%%%%%%%%%%%%%
function dispVolume(vol)

global gSystem;

for iView = 1:3
  % see if we need to redisplay
  if ~isequal(vol.curCoord([vol.xDim(iView) vol.yDim(iView)]),vol.coord([vol.xDim(iView) vol.yDim(iView)]))
    % get the slice
    dispSlice = squeeze(vol.data(vol.viewIndexes{iView,1},vol.viewIndexes{iView,2},vol.viewIndexes{iView,3}));
  
    % clip
    %  dispSlice = clipImage(dispSlice);
  
    % get the correct axis to draw into
    figure(vol.fig(iView));
    a = subplot(1,3,vol.subplotNum(iView));

    % and display the image
    cla(a);imagesc(dispSlice,'Parent',a);
    axis equal
    axis off
  
    % set title
    titleStr = sprintf('%s (%i %i %i)',vol.h.filename,vol.coord(1),vol.coord(2),vol.coord(3));
    h = title(titleStr,'Interpreter','none');
  end
end

%%%%%%%%%%%%%%%%%%%
%    clipImage    %
%%%%%%%%%%%%%%%%%%%
function img = clipImage(img)

% Choose a sensible clipping value
histThresh = length(img(:))/1000;
[cnt, val] = hist(img(:),100);
goodVals = find(cnt>histThresh);
if isempty(goodVals)
  clipMin = 0;clipMax = 0;
else
  clipMin = val(min(goodVals));
  clipMax = val(max(goodVals));
end

% and convert the image
img(img<clipMin) = clipMin;
img(img>clipMax) = clipMax;
if (clipMax-clipMin) > 0
  img = 255*(img-clipMin)./(img-clipMin);
end

%%%%%%%%%%%%%%%%%%%%
%    loadVolume    %
%%%%%%%%%%%%%%%%%%%%
function loadVolume(filename)

global gVols;
global gSystem;

[d h] = mlrLoadImage(filename);
if isempty(d),return,end

% updated number of volumes
gSystem.n = gSystem.n+1;

% update the fields in gVols
gVols(gSystem.n).data = d;
gVols(gSystem.n).h = h;

% set which figure numbers this volume will display into
gVols(gSystem.n).fig(1:3) = gSystem.fig(1:3);
%gVols(gSystem.n).subplotNum(1:3) = gSystem.n;
gVols(gSystem.n).subplotNum(1:3) = 1:3;

% FIX: these should be set with reference to permutation matrix
gVols(gSystem.n).viewDim(1:3) = 1:3;

% for convenience set which dimension is x and y for each of the views
% given the viewDim info set above
for i = 1:3
  otherDims = setdiff(1:3,gVols(gSystem.n).viewDim(i));
  gVols(gSystem.n).xDim(i) = otherDims(1);
  gVols(gSystem.n).yDim(i) = otherDims(2);
end

% update the coordinate (start displaying in middle of volume)
setVolCoord(gSystem.n,round(h.dim(1:3)/2));

% now this sets the indexes from the volume for which the
% image will be displayed
for iView = 1:3
  for jAxis = 1:3
    if jAxis ~= gVols(gSystem.n).viewDim(iView)
      gVols(gSystem.n).viewIndexes{iView,jAxis} = 1:gVols(gSystem.n).h.dim(jAxis);
    else
      gVols(gSystem.n).viewIndexes{iView,jAxis} = gVols(gSystem.n).coord(jAxis);
    end
  end
end

gVols(gSystem.n).curCoord = [nan nan nan];

%%%%%%%%%%%%%%%%%%%%%
%    setVolCoord    %
%%%%%%%%%%%%%%%%%%%%%
function setVolCoord(iVol,coord)

global gVols;

% set the current x,y,z coordinate (start in middle)
gVols(iVol).coord = coord;

% now this sets the indexes from the volume for which the
% image will be displayed
for iView = 1:3
  for jAxis = 1:3
    if jAxis ~= gVols(iVol).viewDim(iView)
      gVols(iVol).viewIndexes{iView,jAxis} = 1:gVols(iVol).h.dim(jAxis);
    else
      gVols(iVol).viewIndexes{iView,jAxis} = gVols(iVol).coord(jAxis);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%    initSystem    %
%%%%%%%%%%%%%%%%%%%%
function initSystem(args)

global gVols;
global gSystem;

% number of loaded volumes
gSystem.n = 0;

% empty vols
vols = [];

% put up the display windows
if 0
gSystem.fig(1) = smartfig('mlrVolCor','reuse');
clf;colormap(gray);axis equal;
gSystem.fig(2) = smartfig('mlrVolSag','reuse');
clf;colormap(gray);axis equal;
gSystem.fig(3) = smartfig('mlrVolAxial','reuse');
clf;colormap(gray);axis equal;

% set the mouse move functions
set(gSystem.fig(1),'WindowButtonMotionFcn','mlrVol(1,1)');
set(gSystem.fig(2),'WindowButtonMotionFcn','mlrVol(1,2)');
set(gSystem.fig(3),'WindowButtonMotionFcn','mlrVol(1,3)');
end

gSystem.buttonWidth = 120;
gSystem.buttonMargin = 10;
gSystem.buttonHeight = 20;
gSystem.buttonLeftMargin = 10;
gSystem.buttonBottomMargin = 10;

gSystem.fig(1) = smartfig('mlrVol','reuse');
gSystem.fig(2) = gSystem.fig(1);
gSystem.fig(3) = gSystem.fig(1);
clf;colormap(gray);
set(gSystem.fig(1),'WindowButtonDownFcn','mlrVol(1,nan)');
set(gSystem.fig(1),'WindowButtonMotionFcn','mlrVol(2,nan)');
gSystem.a(1) = subplot(1,3,1);cla;
gSystem.a(2) = subplot(1,3,2);cla;
gSystem.a(3) = subplot(1,3,3);cla;
set(gSystem.fig(1),'Pointer','fullcrosshair');
gSystem.hTextboxX = makeTextbox(1,'',1,1,1);
gSystem.hTextboxY = makeTextbox(1,'',1,2,1);
gSystem.hTextboxZ = makeTextbox(1,'',1,3,1);

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(viewNum,displayString,rownum,colnum,uisize)

h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(viewNum,rownum,colnum,uisize),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(viewNum,rownum,colnum,uisize)

global gSystem;

% get figure position
figpos = get(gSystem.fig(viewNum),'Position');

% set this buttons width
thisButtonWidth = gSystem.buttonWidth*uisize+(uisize-1)*gSystem.buttonMargin;

% set the position for the button
pos(1) = (gSystem.buttonWidth+gSystem.buttonMargin)*(colnum-1)+gSystem.buttonLeftMargin;
pos(2) = gSystem.buttonBottomMargin + (gSystem.buttonHeight+gSystem.buttonMargin)*(rownum-1)+gSystem.buttonHeight;
pos(3) = thisButtonWidth;
pos(4) = gSystem.buttonHeight;
