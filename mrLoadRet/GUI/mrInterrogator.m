% mrInterrogator.m
%
%      usage: mrInterrogator()
%         by: justin gardner
%       date: 03/14/07
%    purpose: this functions sets up the figure to have an interrogator
%             start by calling
%             mrInterrogator('init',fignum);
%             turn off
%             mrInterrogator('end');
% e.g.
% fignum = figure;
% imagesc(rand(100,150));
% mrInterrogator('init',fignum);
function retval = mrInterrogator(event,fignum,view)

% check arguments
if ~any(nargin == [1 2 3])
  help mrInterrogator
  return
end

% some basic info about location of controls
global gMrInterrogator;
gMrInterrogator.rightMargin = 5;
gMrInterrogator.topMargin = 5;
gMrInterrogator.buttonWidth = 50;
gMrInterrogator.buttonHeight = 20;
gMrInterrogator.margin = 5;
gMrInterrogator.fontsize = 10;
gMrInterrogator.fontname = 'Helvetica';

switch (event)
 case 'init'
   initHandler(fignum,view);
 case 'end'
   endHandler;
 case 'mouseMove'
   mouseMoveHandler;
 case 'mouseUp'
   mouseUpHandler;
 case 'mouseDown'
   mouseDownHandler;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test whether mouse is in image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mouseInImage(xpos,ypos)

global gMrInterrogator;

if isnan(xpos)
  retval = 0;
else
  retval = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousemove
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseMoveHandler(fignum)

global gMrInterrogator;

% get pointer
[x y s] = getMouseCoords;

% check location in bounds on image
if mouseInImage(x,y)
  % set pointer to crosshairs
  set(gMrInterrogator.fignum,'pointer','fullcrosshair');
  % convert to overlay coordinats
  % set the xpos/ypos textbox 
  set(gMrInterrogator.hX,'String',sprintf('x=%i',x));
  set(gMrInterrogator.hY,'String',sprintf('y=%i',y));
else
  % set pointer to arrow
  set(gMrInterrogator.fignum,'pointer','arrow');
  % set strings to empty
  set(gMrInterrogator.hX,'String','');
  set(gMrInterrogator.hY,'String','');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s] = getMouseCoords

global gMrInterrogator;

pointerLoc = get(gMrInterrogator.axesnum,'CurrentPoint');
xpos = round(pointerLoc(1,1));ypos=round(pointerLoc(1,2));

% get the coordinate mapping
global MLR;
view = MLR.views{end};
coords = viewGet(view,'curSliceBaseCoords');

if (xpos >= 1) && (xpos <= size(coords,2)) && (ypos >= 1) && (ypos <= size(coords,1))
  x = coords(ypos,xpos,1);
  y = coords(ypos,xpos,2);
  s = viewGet(view,'currentSlice');
else
  x = nan;
  y = nan;
  s = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseUpHandler(fignum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(fignum)

global gMrInterrogator;

% get pointer
[x y s] = getMouseCoords;

if mouseInImage(x,y)
  % Draw graph
  global MLR;
  view = MLR.views{end};
  overlayNum = viewGet(view,'currentOverlay');
  analysisNum = viewGet(view,'currentAnalysis');
  interrogator = viewGet(view,'interrogator',overlayNum,analysisNum);
  scanNum = viewGet(view,'currentScan');
  feval(interrogator,view,overlayNum,scanNum,x,y,s);
  disp(sprintf('Mouse at %i %i %i',x,y,s));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrInterrogator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler

global gMrInterrogator;

% set the callbacks back to their originals
set(gMrInterrogator.fignum,'WindowButtonMotionFcn',gMrInterrogator.windowButtonMotionFcn);
set(gMrInterrogator.fignum,'WindowButtonDownFcn',gMrInterrogator.windowButtonDownFcn);
set(gMrInterrogator.fignum,'WindowButtonUpFcn',gMrInterrogator.windowButtonUpFcn);

% set the pointer back
set(gMrInterrogator.fignum,'pointer',gMrInterrogator.pointer);

% turn off the text boxes
set(gMrInterrogator.hX,'visible','off');
set(gMrInterrogator.hY,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(fignum,view)

global gMrInterrogator;
gMrInterrogator.fignum = fignum;
gMrInterrogator.guide = guidata(fignum);
figure(fignum);gMrInterrogator.axesnum = gMrInterrogator.guide.axis;

% remember old callbacks
gMrInterrogator.windowButtonMotionFcn = get(fignum,'WindowButtonMotionFcn');
gMrInterrogator.windowButtonDownFcn = get(fignum,'WindowButtonDownFcn');
gMrInterrogator.windowButtonUpFcn = get(fignum,'WindowButtonUpFcn');

% set the callbacks appropriately
set(fignum,'WindowButtonMotionFcn',sprintf('mrInterrogator(''mouseMove'')'));
set(fignum,'WindowButtonDownFcn',sprintf('mrInterrogator(''mouseDown'')'));
set(fignum,'WindowButtonUpFcn',sprintf('mrInterrogator(''mouseUp'')'));

% set pointer to crosshairs
gMrInterrogator.pointer = get(fignum,'pointer');

% set the x and y textbox
gMrInterrogator.hX = makeTextbox('',1,2,1);
gMrInterrogator.hY = makeTextbox('',1,1,1);

% set the x/y min/max
a = axis(gMrInterrogator.axesnum);
gMrInterrogator.xmin = a(1);
gMrInterrogator.xmax = a(2);
gMrInterrogator.ymin = a(3);
gMrInterrogator.ymax = a(4);

% set info for callback
gMrInterrogator.view = view;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(displayString,rownum,colnum,uisize)

global gMrInterrogator;
h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gMrInterrogator.fontsize,'FontName',gMrInterrogator.fontname,'HorizontalAlignment','Center');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(rownum,colnum,uisize)

% get global parameters
global gMrInterrogator;

% get figure position
figpos = get(gMrInterrogator.fignum,'Position');

% set this buttons width
thisButtonWidth = gMrInterrogator.buttonWidth*uisize+(uisize-1)*gMrInterrogator.margin;

% set the position for the button
pos(1) = figpos(3)-gMrInterrogator.margin - (gMrInterrogator.buttonWidth+gMrInterrogator.margin)*(colnum-1)-gMrInterrogator.rightMargin-gMrInterrogator.buttonWidth;
pos(2) = figpos(4)-gMrInterrogator.buttonHeight-gMrInterrogator.topMargin - (gMrInterrogator.buttonHeight+gMrInterrogator.margin)*(rownum-1);
pos(3) = thisButtonWidth;
pos(4) = gMrInterrogator.buttonHeight;
