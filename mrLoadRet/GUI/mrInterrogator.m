% mrInterrogator.m
%
%      usage: mrInterrogator()
%         by: justin gardner
%       date: 03/14/07
%    purpose: this functions sets up the figure to have an interrogator
%             start by calling
%             mrInterrogator('init',viewNum);
%             turn off
%             mrInterrogator('end',viewNum);
function retval = mrInterrogator(event,viewNum,val)

% check arguments
if ~any(nargin == [1 2 3])
  help mrInterrogator
  return
end

% some basic info about location of controls
mrGlobals;
MLR.interrogator{viewNum}.leftMargin = 5;
MLR.interrogator{viewNum}.rightMargin = 5;
MLR.interrogator{viewNum}.topMargin = 5;
MLR.interrogator{viewNum}.bottomMargin = 5;
MLR.interrogator{viewNum}.buttonWidth = 50;
MLR.interrogator{viewNum}.buttonHeight = 20;
MLR.interrogator{viewNum}.margin = 5;
MLR.interrogator{viewNum}.fontsize = 10;
MLR.interrogator{viewNum}.fontname = 'Helvetica';

switch (event)
 case 'init'
   initHandler(viewNum);
 case 'end'
   endHandler(viewNum);
 case 'mouseMove'
   mouseMoveHandler(viewNum);
 case 'mouseUp'
   mouseUpHandler(viewNum);
 case 'mouseDown'
   mouseDownHandler(viewNum);
 case 'interrogator'
   interrogatorHandler(viewNum);
 case 'updateInterrogator'
   updateInterrogatorHandler(viewNum,val);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the interrogator function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateInterrogatorHandler(viewNum,interrogator)

mrGlobals;

% if not a valid function, go back to old one
if exist(interrogator)==2
  set(MLR.interrogator{viewNum}.hInterrogator,'String',interrogator);
  MLR.interrogator{viewNum}.interrogator = interrogator;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change in interrogator field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interrogatorHandler(viewNum)

mrGlobals;

% get new string
interrogator = get(MLR.interrogator{viewNum}.hInterrogator,'String');

% if not a valid function, go back to old one
if exist(interrogator)~=2
  set(MLR.interrogator{viewNum}.hInterrogator,'String',MLR.interrogator{viewNum}.interrogator);
else
  MLR.interrogator{viewNum}.interrogator = interrogator;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test whether mouse is in image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mouseInImage(xpos,ypos)

mrGlobals;

if isnan(xpos)
  retval = 0;
else
  retval = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousemove
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseMoveHandler(viewNum)

mrGlobals;

% get pointer
[x y s] = getMouseCoords(viewNum);

% check location in bounds on image
if mouseInImage(x,y)
  % set pointer to crosshairs
  set(MLR.interrogator{viewNum}.fignum,'pointer','fullcrosshair');
  % convert to overlay coordinats
  % set the xpos/ypos textbox 
  set(MLR.interrogator{viewNum}.hPos,'String',sprintf('[%i %i %i]',x,y,s));
else
  % set pointer to arrow
  set(MLR.interrogator{viewNum}.fignum,'pointer','arrow');
  % set strings to empty
  set(MLR.interrogator{viewNum}.hPos,'String','');
end

% eval the old handler
eval(MLR.interrogator{viewNum}.windowButtonMotionFcn);

% this snippet of code gets the current default interrogator
% function. It shouldn't go in mousemove because it makes mousemvoe
% slow--but more importantly if the users sets the interrogator
% themselves, then they don't want to use the default. But what
% happens if you change views or analyses? Should the interrogator change?
if 0
% check the interrogator
global MLR;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
interrogator = viewGet(view,'interrogator',overlayNum,analysisNum);
% if it is different from current one, then reset it
if ~strcmp(MLR.interrogator{viewNum}.interrogator,interrogator)
  MLR.interrogator{viewNum}.interrogator = interrogator;
  set(MLR.interrogator{viewNum}.hInterrogator,'String',MLR.interrogator{viewNum}.interrogator);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s] = getMouseCoords(viewNum)

mrGlobals;

pointerLoc = get(MLR.interrogator{viewNum}.axesnum,'CurrentPoint');
xpos = round(pointerLoc(1,1));ypos=round(pointerLoc(1,2));

% get the coordinate mapping
view = MLR.views{viewNum};
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
function mouseUpHandler(viewNum)

mrGlobals;

% eval the old handler
eval(MLR.interrogator{viewNum}.windowButtonUpFcn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(viewNum)

mrGlobals;

% get pointer
[x y s] = getMouseCoords(viewNum);

if mouseInImage(x,y)
  global MLR;
  view = MLR.views{viewNum};
  % make a waiting cursor
  %set(MLR.interrogator{viewNum}.fignum,'Pointer','watch');
  % Draw graph
  overlayNum = viewGet(view,'currentOverlay');
  analysisNum = viewGet(view,'currentAnalysis');
  scanNum = viewGet(view,'currentScan');
  feval(MLR.interrogator{viewNum}.interrogator,view,overlayNum,scanNum,x,y,s);
  % reset to full crosshair
  %set(MLR.interrogator{viewNum}.fignum,'Pointer','fullcrosshair');
end

% eval the old handler
eval(MLR.interrogator{viewNum}.windowButtonDownFcn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrInterrogator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(viewNum)

mrGlobals;

% set the callbacks back to their originals
set(MLR.interrogator{viewNum}.fignum,'WindowButtonMotionFcn',MLR.interrogator{viewNum}.windowButtonMotionFcn);
set(MLR.interrogator{viewNum}.fignum,'WindowButtonDownFcn',MLR.interrogator{viewNum}.windowButtonDownFcn);
set(MLR.interrogator{viewNum}.fignum,'WindowButtonUpFcn',MLR.interrogator{viewNum}.windowButtonUpFcn);

% set the pointer back
set(MLR.interrogator{viewNum}.fignum,'pointer',MLR.interrogator{viewNum}.pointer);

% turn off the text boxes
set(MLR.interrogator{viewNum}.hPos,'visible','off');
set(MLR.interrogator{viewNum}.hInterrogator,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(viewNum)

mrGlobals;

fignum = viewGet(MLR.views{viewNum},'figNum');

% see if this is a restart
restart = 0;
if isfield(MLR.interrogator{viewNum},'fignum') && isequal(MLR.interrogator{viewNum}.fignum,fignum)
  disp('(mrInterrogator) Restarting');
  restart = 1;
end

% get figure handles
MLR.interrogator{viewNum}.fignum = fignum;
MLR.interrogator{viewNum}.guide = guidata(fignum);
figure(fignum);MLR.interrogator{viewNum}.axesnum = MLR.interrogator{viewNum}.guide.axis;

if ~restart
  % remember old callbacks
  MLR.interrogator{viewNum}.windowButtonMotionFcn = get(fignum,'WindowButtonMotionFcn');
  MLR.interrogator{viewNum}.windowButtonDownFcn = get(fignum,'WindowButtonDownFcn');
  MLR.interrogator{viewNum}.windowButtonUpFcn = get(fignum,'WindowButtonUpFcn');
end

% set the callbacks appropriately
set(fignum,'WindowButtonMotionFcn',sprintf('mrInterrogator(''mouseMove'',%i)',viewNum));
set(fignum,'WindowButtonDownFcn',sprintf('mrInterrogator(''mouseDown'',%i)',viewNum));
set(fignum,'WindowButtonUpFcn',sprintf('mrInterrogator(''mouseUp'',%i)',viewNum));

% set pointer to crosshairs
MLR.interrogator{viewNum}.pointer = get(fignum,'pointer');

if ~restart
  % set the x and y textbox
  MLR.interrogator{viewNum}.hPos = makeTextbox(viewNum,'',1,4,2);
  MLR.interrogator{viewNum}.hInterrogator = makeTextentry(viewNum,'test','interrogator',1,1,3);
else
  set(MLR.interrogator{viewNum}.hPos,'visible','on');
  set(MLR.interrogator{viewNum}.hInterrogator,'visible','on');
end

% set the x/y min/max
a = axis(MLR.interrogator{viewNum}.axesnum);
MLR.interrogator{viewNum}.xmin = a(1);
MLR.interrogator{viewNum}.xmax = a(2);
MLR.interrogator{viewNum}.ymin = a(3);
MLR.interrogator{viewNum}.ymax = a(4);

% set info for callback
MLR.interrogator{viewNum}.viewNum = viewNum;

% set interrogator field
global MLR;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
MLR.interrogator{viewNum}.interrogator = viewGet(view,'interrogator',overlayNum,analysisNum);
set(MLR.interrogator{viewNum}.hInterrogator,'String',MLR.interrogator{viewNum}.interrogator);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(viewNum,displayString,rownum,colnum,uisize)

mrGlobals;
h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(viewNum,rownum,colnum,uisize),'FontSize',MLR.interrogator{viewNum}.fontsize,'FontName',MLR.interrogator{viewNum}.fontname,'HorizontalAlignment','Center');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextentry(viewNum,displayString,callback,rownum,colnum,uisize)

% make callback string
if isnumeric(callback)
  callback = sprintf('mrInterrogator(%f,%i)',callback,viewNum);
else
  callback = sprintf('mrInterrogator(''%s'',%i)',callback,viewNum);
end  

mrGlobals;

h = uicontrol('Style','edit','Callback',callback,'String',displayString,'Position',getUIControlPos(viewNum,rownum,colnum,uisize),'FontSize',MLR.interrogator{viewNum}.fontsize,'FontName',MLR.interrogator{viewNum}.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(viewNum,rownum,colnum,uisize)

% get global parameters
mrGlobals;

% get figure position
figpos = get(MLR.interrogator{viewNum}.fignum,'Position');

% set this buttons width
thisButtonWidth = MLR.interrogator{viewNum}.buttonWidth*uisize+(uisize-1)*MLR.interrogator{viewNum}.margin;

% set the position for the button
%pos(1) = figpos(3)-MLR.interrogator{viewNum}.margin - (MLR.interrogator{viewNum}.buttonWidth+MLR.interrogator{viewNum}.margin)*(colnum-1)-MLR.interrogator{viewNum}.rightMargin-MLR.interrogator{viewNum}.buttonWidth;
%pos(2) = figpos(4)-MLR.interrogator{viewNum}.buttonHeight-MLR.interrogator{viewNum}.topMargin - (MLR.interrogator{viewNum}.buttonHeight+MLR.interrogator{viewNum}.margin)*(rownum-1);
pos(1) = (MLR.interrogator{viewNum}.buttonWidth+MLR.interrogator{viewNum}.margin)*(colnum-1)+MLR.interrogator{viewNum}.leftMargin;
pos(2) = MLR.interrogator{viewNum}.bottomMargin + (MLR.interrogator{viewNum}.buttonHeight+MLR.interrogator{viewNum}.margin)*(rownum-1)+MLR.interrogator{viewNum}.buttonHeight;
pos(3) = thisButtonWidth;
pos(4) = MLR.interrogator{viewNum}.buttonHeight;
