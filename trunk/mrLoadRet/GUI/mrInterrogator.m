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
function retval = mrInterrogator(event,fignum,viewNum)

% check arguments
if ~any(nargin == [1 2 3])
  help mrInterrogator
  return
end

% some basic info about location of controls
global gMrInterrogator;
gMrInterrogator.leftMargin = 5;
gMrInterrogator.rightMargin = 5;
gMrInterrogator.topMargin = 5;
gMrInterrogator.bottomMargin = 5;
gMrInterrogator.buttonWidth = 50;
gMrInterrogator.buttonHeight = 20;
gMrInterrogator.margin = 5;
gMrInterrogator.fontsize = 10;
gMrInterrogator.fontname = 'Helvetica';

switch (event)
 case 'init'
   initHandler(fignum,viewNum);
 case 'end'
   endHandler;
 case 'mouseMove'
   mouseMoveHandler;
 case 'mouseUp'
   mouseUpHandler;
 case 'mouseDown'
   mouseDownHandler;
 case 'interrogator'
   interrogatorHandler;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change in interrogator field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interrogatorHandler

global gMrInterrogator;

% get new string
interrogator = get(gMrInterrogator.hInterrogator,'String');

% if not a valid function, go back to old one
if exist(interrogator)~=2
  set(gMrInterrogator.hInterrogator,'String',gMrInterrogator.interrogator);
else
  gMrInterrogator.interrogator = interrogator;
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
  set(gMrInterrogator.hPos,'String',sprintf('[%i %i %i]',x,y,s));
else
  % set pointer to arrow
  set(gMrInterrogator.fignum,'pointer','arrow');
  % set strings to empty
  set(gMrInterrogator.hPos,'String','');
end

% eval the old handler
eval(gMrInterrogator.windowButtonMotionFcn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s] = getMouseCoords

global gMrInterrogator;

pointerLoc = get(gMrInterrogator.axesnum,'CurrentPoint');
xpos = round(pointerLoc(1,1));ypos=round(pointerLoc(1,2));

% get the coordinate mapping
global MLR;
view = MLR.views{gMrInterrogator.viewNum};
coords = viewGet(view,'curSliceOverlayCoords');
if isempty(coords)
  coords = viewGet(view,'curSliceBaseCoords');
  % FIX FIX FIX this assumes 1/2 coords, there 
  % should always be a coords that knows about 
  % the epi, not just the overlay
  if ~isempty(coords)
    coords(:,:,1) = round(coords(:,:,1)/2);
    coords(:,:,2) = round(coords(:,:,2)/2);
  end
end

if ~isempty(coords) && (xpos >= 1) && (xpos <= size(coords,2)) && (ypos >= 1) && (ypos <= size(coords,1))
  x = round(coords(ypos,xpos,1));
  y = round(coords(ypos,xpos,2));
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

global gMrInterrogator;

% eval the old handler
eval(gMrInterrogator.windowButtonUpFcn);

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
  view = MLR.views{gMrInterrogator.viewNum};
  overlayNum = viewGet(view,'currentOverlay');
  analysisNum = viewGet(view,'currentAnalysis');
  scanNum = viewGet(view,'currentScan');
  feval(gMrInterrogator.interrogator,view,overlayNum,scanNum,x,y,s);
end

% eval the old handler
eval(gMrInterrogator.windowButtonDownFcn);

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
set(gMrInterrogator.hPos,'visible','off');
set(gMrInterrogator.hInterrogator,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(fignum,viewNum)

global gMrInterrogator;

% see if this is a restart
restart = 0;
if isfield(gMrInterrogator,'fignum') && isequal(gMrInterrogator.fignum,fignum)
  disp('(mrInterrogator) Restarting');
  restart = 1;
end

% get figure handles
gMrInterrogator.fignum = fignum;
gMrInterrogator.guide = guidata(fignum);
figure(fignum);gMrInterrogator.axesnum = gMrInterrogator.guide.axis;

if ~restart
  % remember old callbacks
  gMrInterrogator.windowButtonMotionFcn = get(fignum,'WindowButtonMotionFcn');
  gMrInterrogator.windowButtonDownFcn = get(fignum,'WindowButtonDownFcn');
  gMrInterrogator.windowButtonUpFcn = get(fignum,'WindowButtonUpFcn');
end

% set the callbacks appropriately
set(fignum,'WindowButtonMotionFcn',sprintf('mrInterrogator(''mouseMove'')'));
set(fignum,'WindowButtonDownFcn',sprintf('mrInterrogator(''mouseDown'')'));
set(fignum,'WindowButtonUpFcn',sprintf('mrInterrogator(''mouseUp'')'));

% set pointer to crosshairs
gMrInterrogator.pointer = get(fignum,'pointer');

if ~restart
  % set the x and y textbox
  gMrInterrogator.hPos = makeTextbox('',1,4,2);
  gMrInterrogator.hInterrogator = makeTextentry('test','interrogator',1,1,3);
else
  set(gMrInterrogator.hPos,'visible','on');
  set(gMrInterrogator.hInterrogator,'visible','on');
end

% set the x/y min/max
a = axis(gMrInterrogator.axesnum);
gMrInterrogator.xmin = a(1);
gMrInterrogator.xmax = a(2);
gMrInterrogator.ymin = a(3);
gMrInterrogator.ymax = a(4);

% set info for callback
gMrInterrogator.viewNum = viewNum;

% set interrogator field
global MLR;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
gMrInterrogator.interrogator = viewGet(view,'interrogator',overlayNum,analysisNum);
set(gMrInterrogator.hInterrogator,'String',gMrInterrogator.interrogator);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(displayString,rownum,colnum,uisize)

global gMrInterrogator;
h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gMrInterrogator.fontsize,'FontName',gMrInterrogator.fontname,'HorizontalAlignment','Center');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextentry(displayString,callback,rownum,colnum,uisize)

% make callback string
if isnumeric(callback)
  callback = sprintf('mrInterrogator(%f)',callback);
else
  callback = sprintf('mrInterrogator(''%s'')',callback);
end  

global gMrInterrogator;

h = uicontrol('Style','edit','Callback',callback,'String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gMrInterrogator.fontsize,'FontName',gMrInterrogator.fontname);

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
%pos(1) = figpos(3)-gMrInterrogator.margin - (gMrInterrogator.buttonWidth+gMrInterrogator.margin)*(colnum-1)-gMrInterrogator.rightMargin-gMrInterrogator.buttonWidth;
%pos(2) = figpos(4)-gMrInterrogator.buttonHeight-gMrInterrogator.topMargin - (gMrInterrogator.buttonHeight+gMrInterrogator.margin)*(rownum-1);
pos(1) = (gMrInterrogator.buttonWidth+gMrInterrogator.margin)*(colnum-1)+gMrInterrogator.leftMargin;
pos(2) = gMrInterrogator.bottomMargin + (gMrInterrogator.buttonHeight+gMrInterrogator.margin)*(rownum-1)+gMrInterrogator.buttonHeight;
pos(3) = thisButtonWidth;
pos(4) = gMrInterrogator.buttonHeight;
