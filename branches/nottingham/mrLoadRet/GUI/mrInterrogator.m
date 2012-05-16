% mrInterrogator.m
%
%      usage: mrInterrogator()
%         by: justin gardner
%       date: 03/14/07
%        $Id$
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
MLR.interrogator{viewNum}.leftMargin = .015;
MLR.interrogator{viewNum}.rightMargin = .02;
MLR.interrogator{viewNum}.topMargin = .005;
MLR.interrogator{viewNum}.bottomMargin = .01;
MLR.interrogator{viewNum}.buttonWidth = .065;
MLR.interrogator{viewNum}.buttonHeight = .035;
MLR.interrogator{viewNum}.margin = 0;
MLR.interrogator{viewNum}.fontsize = 10;
MLR.interrogator{viewNum}.fontname = 'Helvetica';

switch (event)
    case 'init'
        initHandler(viewNum);
    case 'end'
        endHandler(viewNum);
    case 'isactive'
        retval = isActiveHandler(viewNum);
    case 'mouseMove'
        mouseMoveHandler(viewNum);
    case 'mouseUp'
        mouseUpHandler(viewNum);
    case 'mouseDown'
        mouseDownHandler(viewNum);
    case 'interrogator'
        interrogatorHandler(viewNum);
    case 'defaultInterrogators'
        defaultInterrogatorsHandler(viewNum);
    case 'updateInterrogator'
        updateInterrogatorHandler(viewNum,val);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the interrogator function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateInterrogatorHandler(viewNum,interrogator)

if isActiveHandler(viewNum)
  mrGlobals;
  v = MLR.views{viewNum};

  % set list of interrogators
  interrogatorList = getDefaultInterrogators(v);
  set(MLR.interrogator{viewNum}.hInterrogatorLabel,'String',interrogatorList);

  % if not a valid function, go back to old one
  if isempty(which(interrogator))
    set(MLR.interrogator{viewNum}.hInterrogator,'String',interrogator);
    MLR.interrogator{viewNum}.interrogator = interrogator;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change in interrogator field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interrogatorHandler(viewNum)

mrGlobals;

% get new string
interrogator = get(MLR.interrogator{viewNum}.hInterrogator,'String');

% if not a valid function, go back to old one
if isfield(MLR.interrogator{viewNum},'hInterrogator')
  if isempty(which(interrogator))
    set(MLR.interrogator{viewNum}.hInterrogator,'String',MLR.interrogator{viewNum}.interrogator);
  else
    MLR.interrogator{viewNum}.interrogator = interrogator;
    % and set the interrogator for the overlay
    global MLR;
    v= MLR.views{viewNum};
    v = viewSet(v,'interrogator',interrogator);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change in interrogator field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defaultInterrogatorsHandler(viewNum)

mrGlobals;

% get which item was chosen
value = get(MLR.interrogator{viewNum}.hInterrogatorLabel,'Value');
if value > 1
  % get the interrogator string
  interrogatorList = get(MLR.interrogator{viewNum}.hInterrogatorLabel,'String');
  newInterrogator = interrogatorList{value};
  % set the popupmenu back to the top value
  set(MLR.interrogator{viewNum}.hInterrogatorLabel,'Value',1);
  % see if the interrogator exists
  if ~isempty(which(newInterrogator))
    % set the string
    set(MLR.interrogator{viewNum}.hInterrogator,'String',newInterrogator);
    MLR.interrogator{viewNum}.interrogator = newInterrogator;
  end
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
[x y s xBase yBase sBase xTal yTal zTal] = getMouseCoords(viewNum);

% check location in bounds on image
if mouseInImage(x,y)
    % set pointer to crosshairs
    set(MLR.interrogator{viewNum}.fignum,'pointer','fullcrosshair');
    % set the xpos/ypos textbox
    set(MLR.interrogator{viewNum}.hPos,'String',sprintf('[%i %i %i]',x,y,s));
    % set the overlay value
    overlayVal = viewGet(MLR.interrogator{viewNum}.viewNum,'overlayDataVal',x,y,s);
    set(MLR.interrogator{viewNum}.hOverlayLabel,'String',viewGet(MLR.interrogator{viewNum}.viewNum,'overlayName'));
    if ~isempty(overlayVal)
      set(MLR.interrogator{viewNum}.hOverlay,'String',sprintf('[%f]',overlayVal));
    else
      set(MLR.interrogator{viewNum}.hOverlay,'String','');
    end
else
    % set pointer to arrow
    set(MLR.interrogator{viewNum}.fignum,'pointer','arrow');
    % set strings to empty
    set(MLR.interrogator{viewNum}.hPos,'String','');
    set(MLR.interrogator{viewNum}.hOverlay,'String','');
end
if mouseInImage(xBase,yBase)
  set(MLR.interrogator{viewNum}.hPosBase,'String',sprintf('[%0.4g %0.4g %0.4g]',xBase,yBase,sBase));
else
  set(MLR.interrogator{viewNum}.hPosBase,'String','');
end

if mouseInImage(xTal, yTal)
  set(MLR.interrogator{viewNum}.hPosTalLabel,'visible','on');
  set(MLR.interrogator{viewNum}.hPosTalLabel,'String','Tal');
  set(MLR.interrogator{viewNum}.hPosTal,'visible','on');
  set(MLR.interrogator{viewNum}.hPosTal,'String',sprintf('[%0.4g %0.4g %0.4g]',xTal,yTal,zTal));
% % else
% %   set(MLR.interrogator{viewNum}.hPosTalLabel,'visible','off');
% %   set(MLR.interrogator{viewNum}.hPosTal,'visible','off');
% % end
% % 
% % if ~mouseInImage(xTal,yTal) && mouseInImage(xBase,yBase) % display magnet coordinates as well
elseif mouseInImage(xBase,yBase) % display magnet coordinates as well
  % get the *original* base2mag 
  base2mag = viewGet(MLR.interrogator{viewNum}.viewNum,'baseQform');
  base2mag = base2mag * shiftOriginXform;
  magCoords = base2mag*[xBase yBase sBase 1]';
  set(MLR.interrogator{viewNum}.hPosTalLabel,'visible','on');
  set(MLR.interrogator{viewNum}.hPosTalLabel,'String','Magnet');
  set(MLR.interrogator{viewNum}.hPosTal,'visible','on');
  set(MLR.interrogator{viewNum}.hPosTal,'String',sprintf('[%0.1f %0.1f %0.1f]',magCoords(1),magCoords(2),magCoords(3)));
else
	set(MLR.interrogator{viewNum}.hPosTal,'visible','off');
end


% eval the old handler
%eval(MLR.interrogator{viewNum}.windowButtonMotionFcn);

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
function [x y s xBase yBase sBase xTal yTal zTal] = getMouseCoords(viewNum)

% get the view
mrGlobals;
view = viewGet([],'view',viewNum);

% no bases
if viewGet(view,'numBase') == 0
  x = nan;y= nan;s = nan;xBase = nan;yBase = nan;sBase = nan;xTal = nan;yTal = nan;zTal = nan;
  return
end
baseType = viewGet(view,'baseType');

% get location of pointer
pointerLoc = get(MLR.interrogator{viewNum}.axesnum,'CurrentPoint');

if baseType <= 1
  mouseY = round(pointerLoc(1,1));
  mouseX = round(pointerLoc(1,2));

  % get base coordinates
  baseCoords = viewGet(view,'cursliceBaseCoords');
  % convert mouse to baseCoords
  if (mouseX>0) && (mouseX<=size(baseCoords,1)) && (mouseY>0) && (mouseY<=size(baseCoords,2))
    xBase = baseCoords(mouseX,mouseY,1);
    yBase = baseCoords(mouseX,mouseY,2);
    sBase = baseCoords(mouseX,mouseY,3);
  else
    x = nan;y = nan; s = nan;
    xBase = nan;yBase = nan; sBase = nan;
    xTal = nan; yTal = nan; zTal = nan;
    return
  end
else
  % handle getting coordinates for surface
  baseSurface = viewGet(view,'baseSurface');
  baseDims = viewGet(view,'baseSurfaceDims');
  pos = [];xBase = nan; yBase = nan; sBase = nan;
  % check mouse bounding box coords against baseDims
  % for a quick check to see if we are in the volume
  if all(pointerLoc(1,:) >= 0)
    % then use select3d which is slooow, but accurate
    hobj = get(MLR.interrogator{viewNum}.axesnum,'Children');
    % make sure we are using the correct object (should be the 3D
    % brain). Use end here because with searchForVoxel we plot a
    % point on the image and the brain object always seems to be
    % last. But if this is not always the case, then we may need
    % to do a little more work here to find the correct object
    [pos v vi] = select3d(hobj(end));
    % convert the index to the coordinates
    if ~isempty(pos)
      baseCoordMap = viewGet(view,'baseCoordMap');
      %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
      corticalSlice = max(1,ceil(mean(viewGet(view,'corticalDepth'))*size(baseCoordMap.coords,5)));
      pos = round(squeeze(baseCoordMap.coords(1,vi,1,:,corticalSlice)));
      %pos = round(squeeze(baseCoordMap.coords(1,vi,1,:)));
      xBase = pos(1);yBase = pos(2);sBase = pos(3);
    end
  end      
end

% transform from base coordinates to talairach coordinates, if the base has a talairach transform defined
base2tal = viewGet(view,'base2tal'); % keyboard
if(~isempty(base2tal))
  talCoords = round(base2tal * [xBase yBase sBase 1]');
  xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
else
  xTal = nan; yTal = nan; zTal = nan;
end


% transform from base coordinates into scan coordinates
base2scan = viewGet(view,'base2scan');
if isempty(base2scan), x = nan; y=nan; s=nan; return,  end
transformed = round(base2scan*[xBase yBase sBase 1]');

x = transformed(1);
y = transformed(2);
s = transformed(3);

% get the scan dims to make sure we haven't jumped off end
% of scan
scanDims = viewGet(view,'scanDims',viewGet(view,'curScan'));
if ((x < 1) || (x > scanDims(1)) || ...
        (y < 1) || (y > scanDims(2)) || ...
        (s < 1) || (s > scanDims(3)))
    x = nan;y = nan;s = nan;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseUpHandler(viewNum)

mrGlobals;

% eval the old handler
%eval(MLR.interrogator{viewNum}.windowButtonUpFcn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousedown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(viewNum)

mrGlobals;

% get pointer
[x y s xBase yBase sBase xTal yTal zTal] = getMouseCoords(viewNum);

% set the points in the global, so that they
% are accessible by doing a viewGet
MLR.interrogator{viewNum}.mouseDownBaseCoords = [xBase yBase sBase];
MLR.interrogator{viewNum}.mouseDownScanCoords = [x y s];

if mouseInImage(x,y)
    global MLR;
    if isempty(which(MLR.interrogator{viewNum}.interrogator))
      disp(sprintf('(mrInterrogator) Cannot find interrogator function'));
      return
    end
    view = MLR.views{viewNum};
    % make a waiting cursor
    %set(MLR.interrogator{viewNum}.fignum,'Pointer','watch');
    % find all rois that the user clicked on
    roi = {};roinums = [];
    switch lower(viewGet(view,'showROIs'))
     case {'hide'}
      roinums = [];
     case {'selected','selected perimeter'}
      roinums = viewGet(view,'currentROI');
     case {'all','all perimeter'}
      roinums = 1:viewGet(view,'nROIs');
     case {'group','group perimeter'}
      roinums = viewGet(view,'roiGroup');
    end
    for roinum = roinums
      roicoords = getROICoordinates(view,roinum,0);
      % see if this is a matching roi
      if ismember(round([xBase yBase sBase]),roicoords','rows')
	% get the roi
	roi{end+1} = viewGet(view,'roi',roinum);
      end
    end
    % get some info
    overlayNum = viewGet(view,'currentOverlay');
    analysisNum = viewGet(view,'currentAnalysis');
    scanNum = viewGet(view,'currentScan');
    % call interrogator
    feval(MLR.interrogator{viewNum}.interrogator,view,overlayNum,scanNum,x,y,s,roi);
    % reset to full crosshair
    %set(MLR.interrogator{viewNum}.fignum,'Pointer','fullcrosshair');
end

% eval the old handler
%eval(MLR.interrogator{viewNum}.windowButtonDownFcn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrInterrogator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(viewNum)

mrGlobals;
v = MLR.views{viewNum};
fignum = viewGet(v,'figNum');

% check to see if it is running
if ~isfield(MLR.interrogator{viewNum},'fignum') || ~isequal(MLR.interrogator{viewNum}.fignum,fignum)
  return
end

% set the callbacks back to their originals
set(MLR.interrogator{viewNum}.fignum,'WindowButtonMotionFcn',MLR.interrogator{viewNum}.windowButtonMotionFcn);
set(MLR.interrogator{viewNum}.fignum,'WindowButtonDownFcn',MLR.interrogator{viewNum}.windowButtonDownFcn);
set(MLR.interrogator{viewNum}.fignum,'WindowButtonUpFcn',MLR.interrogator{viewNum}.windowButtonUpFcn);

% set the pointer back
set(MLR.interrogator{viewNum}.fignum,'pointer',MLR.interrogator{viewNum}.pointer);

% turn off the text boxes
set(MLR.interrogator{viewNum}.hPos,'visible','off');
set(MLR.interrogator{viewNum}.hPosLabel,'visible','off');
set(MLR.interrogator{viewNum}.hPosBase,'visible','off');
set(MLR.interrogator{viewNum}.hPosBaseLabel,'visible','off');
set(MLR.interrogator{viewNum}.hPosTal,'visible','off');
set(MLR.interrogator{viewNum}.hPosTalLabel,'visible','off');
set(MLR.interrogator{viewNum}.hOverlay,'visible','off');
set(MLR.interrogator{viewNum}.hOverlayLabel,'visible','off');
set(MLR.interrogator{viewNum}.hInterrogator,'visible','off');
set(MLR.interrogator{viewNum}.hInterrogatorLabel,'visible','off');

% turn on free rotation
if viewGet(v,'baseType') == 2
  mlrSetRotate3d(v,1);
end

% make inactive
MLR.interrogator{viewNum}.isActive = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(viewNum)

mrGlobals;
v = MLR.views{viewNum};
fignum = viewGet(MLR.views{viewNum},'figNum');

% see if this is a restart
restart = 0;
if isfield(MLR.interrogator{viewNum},'fignum') && isequal(MLR.interrogator{viewNum}.fignum,fignum)
    disp('(mrInterrogator) Restarting');
    restart = 1;
end

% turn off free rotation
mlrSetRotate3d(v,0);

% get figure handles
MLR.interrogator{viewNum}.fignum = fignum;
MLR.interrogator{viewNum}.guide = guidata(fignum);
figure(fignum);MLR.interrogator{viewNum}.axesnum = MLR.interrogator{viewNum}.guide.axis;
MLR.interrogator{viewNum}.isActive = 1;

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

% get default interrogators
interrogatorList = getDefaultInterrogators(v);

if ~restart
    % set the x and y textbox
    MLR.interrogator{viewNum}.hInterrogatorLabel = makePopupmenu(viewNum,interrogatorList,'defaultInterrogators',1,1,2);
    MLR.interrogator{viewNum}.hInterrogator = makeTextentry(viewNum,'test','interrogator',1,3,3);
    MLR.interrogator{viewNum}.hPos = makeTextbox(viewNum,'',1,6,2);
    MLR.interrogator{viewNum}.hPosBase = makeTextbox(viewNum,'',1,8,2);
    MLR.interrogator{viewNum}.hPosTal = makeTextbox(viewNum,'',1,10,2);
    MLR.interrogator{viewNum}.hOverlay = makeTextbox(viewNum,'',1,12,4);
    MLR.interrogator{viewNum}.hPosLabel = makeTextbox(viewNum,'Scan',2,6,2);
    MLR.interrogator{viewNum}.hPosBaseLabel = makeTextbox(viewNum,'Base',2,8,2);
    MLR.interrogator{viewNum}.hPosTalLabel = makeTextbox(viewNum,'Talairach',2,10,2);
    MLR.interrogator{viewNum}.hOverlayLabel = makeTextbox(viewNum,'Overlay',2,12,4);
else
    set(MLR.interrogator{viewNum}.hPos,'visible','on');
    set(MLR.interrogator{viewNum}.hPosBase,'visible','on');
    set(MLR.interrogator{viewNum}.hPosTal,'visible','on');
    set(MLR.interrogator{viewNum}.hOverlay,'visible','on');
    set(MLR.interrogator{viewNum}.hPosLabel,'visible','on');
    set(MLR.interrogator{viewNum}.hPosBaseLabel,'visible','on');
    set(MLR.interrogator{viewNum}.hPosTalLabel,'visible','on');
    set(MLR.interrogator{viewNum}.hOverlayLabel,'visible','on');
    set(MLR.interrogator{viewNum}.hInterrogator,'visible','on');
    set(MLR.interrogator{viewNum}.hInterrogatorLabel,'visible','on');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = isActiveHandler(viewNum)

global MLR
if isfield(MLR,'interrogator') && (length(MLR.interrogator) >= viewNum) && isfield(MLR.interrogator{viewNum},'isActive')
  retval =  MLR.interrogator{viewNum}.isActive;
else
  retval = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(viewNum,displayString,rownum,colnum,uisize)

mrGlobals;
h = uicontrol('Style','text','Units','normalized','String',displayString,'Position',getUIControlPos(viewNum,rownum,colnum,uisize),...
  'FontSize',MLR.interrogator{viewNum}.fontsize,'FontName',MLR.interrogator{viewNum}.fontname,'HorizontalAlignment','Center');

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

h = uicontrol('Style','edit','Units','normalized','Callback',callback,'String',displayString,...
              'Position',getUIControlPos(viewNum,rownum,colnum,uisize)+[0 .005 0 0],...
              'FontSize',MLR.interrogator{viewNum}.fontsize,'FontName',MLR.interrogator{viewNum}.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makePopupmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makePopupmenu(viewNum,displayString,callback,rownum,colnum,uisize)

% make callback string
if isnumeric(callback)
  callback = sprintf('mrInterrogator(%f,%i)',callback,viewNum);
else
  callback = sprintf('mrInterrogator(''%s'',%i)',callback,viewNum);
end

if ~iscell(displayString)
  choices{1} = displayString;
else
  if iscell(displayString{1})
    choices = displayString{1};
  else
    choices = displayString;
  end
end

mrGlobals;
h = uicontrol('Style','Popupmenu','Units','normalized','Callback',callback,'Max',length(choices),'Min',1,'String',choices,'Value',1,'Position',getUIControlPos(viewNum,rownum,colnum,uisize),'FontSize',MLR.interrogator{viewNum}.fontsize,'FontName',MLR.interrogator{viewNum}.fontname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(viewNum,rownum,colnum,uisize)

% get global parameters
mrGlobals;

% set this buttons width
thisButtonWidth = MLR.interrogator{viewNum}.buttonWidth*uisize+(uisize-1)*MLR.interrogator{viewNum}.margin;

% set the position for the button
%pos(1) = figpos(3)-MLR.interrogator{viewNum}.margin - (MLR.interrogator{viewNum}.buttonWidth+MLR.interrogator{viewNum}.margin)*(colnum-1)-MLR.interrogator{viewNum}.rightMargin-MLR.interrogator{viewNum}.buttonWidth;
%pos(2) = figpos(4)-MLR.interrogator{viewNum}.buttonHeight-MLR.interrogator{viewNum}.topMargin - (MLR.interrogator{viewNum}.buttonHeight+MLR.interrogator{viewNum}.margin)*(rownum-1);
pos(1) = (MLR.interrogator{viewNum}.buttonWidth+MLR.interrogator{viewNum}.margin)*(colnum-1)+MLR.interrogator{viewNum}.leftMargin;
pos(2) = MLR.interrogator{viewNum}.bottomMargin + (MLR.interrogator{viewNum}.buttonHeight+MLR.interrogator{viewNum}.margin)*(rownum-1);%+MLR.interrogator{viewNum}.buttonHeight;
pos(3) = thisButtonWidth;
pos(4) = MLR.interrogator{viewNum}.buttonHeight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get default interrogators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interrogatorList = getDefaultInterrogators(v);

interrogatorList = mrGetPref('defaultInterrogators');
if isstr(interrogatorList)
  interrogatorList = commaDelimitedToCell(interrogatorList);
end

%get names of interrogators in interrogators directory
interrogatorsDirectory = which('mrLoadRetGUI');
interrogatorsDirectory = [interrogatorsDirectory(1:strfind(interrogatorsDirectory,'GUI/mrLoadRetGUI.m')-1) 'Plugins/Interrogators/'];
interrogatorFiles =  dir([interrogatorsDirectory '*.m']);
for iFile=1:length(interrogatorFiles)
   interrogatorList{end+1} = stripext(interrogatorFiles(iFile).name);
end

% put the interrogator associated with this overlay on the list
overlayInterrogator = viewGet(v,'interrogator');
if ~isempty(overlayInterrogator)
  interrogatorList{end+1} = overlayInterrogator;
end
interrogatorList = putOnTopOfList('Interrogator',interrogatorList);

% add default interrogators from viewGet. These are per-session
% ones that are added by mlrAdjustGUI
defaultInterrogators = viewGet(v,'defaultInterrogators');
if ~isempty(defaultInterrogators)
  interrogatorList = {interrogatorList{:} defaultInterrogators{:}};
end