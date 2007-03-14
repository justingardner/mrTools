% mrDefaultParamsGUI.m
%
%      usage: mrDefaultParamsGUI()
%         by: justin gardner
%       date: 03/13/07
%    purpose: creates a dialog for selection of parameters
%             see wiki for details
%
function params = mrDefaultParamsGUI(varargin)

% check arguments
if ~any(nargin == [1 2])
  help mrDefaultParamsGUI
  return
end

% if this is a cell array, it means to open up the figure
% using the variable name, default value pairs given
if (nargin == 1) && iscell(varargin{1})
  params = initFigure(varargin{1});
% otherwise it is a callback
else
  handleCallbacks(varargin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up figure in first palce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = initFigure(vars)

global gParams;

% parse the input parameter string
[vars varinfo] = mrDefaultParamsParse(vars);
gParams.varinfo = varinfo;

% get variable names
gParams.vars = vars;

% some basic info about location of controls
gParams.leftMargin = 10;
gParams.topMargin = 10;
gParams.buttonWidth = 100;
gParams.buttonHeight = 20;
gParams.margin = 5;
gParams.fontsize = 12;
gParams.fontname = 'Helvetica';

% get the figure 
if ~isfield(gParams,'fignum') || (gParams.fignum == -1);
  % open figure, and turn off menu
  gParams.fignum = figure;
  set(gParams.fignum,'MenuBar','none');
  set(gParams.fignum,'NumberTitle','off');
  set(gParams.fignum,'Name','Set parameters');
else
  figure(gParams.fignum);
end

numrows = length(gParams.varinfo)+1;
numcols = 4;
% set height of figure according to how many rows we have
if ~isfield(gParams,'figpos')
  figpos = get(gParams.fignum,'Position');
else
  figpos = gParams.figpos;
end
figpos(4) = 2*gParams.topMargin+numrows*gParams.buttonHeight+(numrows-1)*gParams.margin;
figpos(3) = 2*gParams.leftMargin+numcols*gParams.buttonWidth+(numcols-1)*gParams.margin;
set(gParams.fignum,'Position',figpos);

% make entry buttons
for i = 1:length(gParams.varinfo)
  gParams.ui.varname(i) = makeTextbox(gParams.varinfo{i}.name,i,1,1);
  if isfield(gParams.varinfo{i},'incdec')
    gParams.ui.varentry(i) = makeTextentryWithIncdec(gParams.varinfo{i}.value,i,i,2,3);
  elseif strcmp(gParams.varinfo{i}.type,'checkbox')
    gParams.ui.varentry(i) = makeCheckbox(gParams.varinfo{i}.value,i,2,.25);
  elseif strcmp(gParams.varinfo{i}.type,'popupmenu') || iscell(gParams.varinfo{i}.value)
    gParams.ui.varentry(i) = makePopupmenu(gParams.varinfo{i}.value,i,2,3);
  else
    gParams.ui.varentry(i) = makeTextentry(gParams.varinfo{i}.value,i,i,2,3);
  end
end

% make ok and cancel buttons
makeButton('Help','help',numrows,1,1);
makeButton('OK','ok',numrows,numcols,1);
makeButton('Cancel','cancel',numrows,numcols-1,1);

% wait for user to hit ok or cancel (which sets uiresume)
uiwait;

% check return value
if gParams.ok
  % return the var entries
  for i = 1:length(gParams.varinfo)
    if strcmp(gParams.varinfo{i}.type,'checkbox')
      params.(gParams.varinfo{i}.name) = num2str(get(gParams.ui.varentry(i),'Value'));
    elseif strcmp(gParams.varinfo{i}.type,'popupmenu')
      choices = get(gParams.ui.varentry(i),'String');
      params.(gParams.varinfo{i}.name) = choices{get(gParams.ui.varentry(i),'Value')};
    else
      params.(gParams.varinfo{i}.name) = get(gParams.ui.varentry(i),'String');
    end
    % if non numeric then convert back to a number
    if ~any(strcmp(gParams.varinfo{i}.type,{'string' 'popupmenu'}))
      params.(gParams.varinfo{i}.name) = str2num(params.(gParams.varinfo{i}.name));
    end
  end
  params.paramInfo = vars;
else
  % otherwise return empty
  params = [];
end

% close figure
figpos = get(gParams.fignum,'Position');
close(gParams.fignum);
% close help
helpcloseHandler;
helpFigpos = gParams.helpFigpos;

clear global gParams;
global gParams;
gParams.figpos = figpos;
gParams.helpFigpos = helpFigpos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handleCallbacks(args)

event = args{1};

% if it is a number than an entry field has been updated
if ~isstr(event)
  if length(args) == 1
    buttonHandler(event);
  else
    buttonHandler(event,args{2});
  end
else
  switch lower(event)
   case {'ok'}
    okHandler;
   case {'help'}
    helpHandler;
   case {'helpclose'}
    helpcloseHandler;
   case {'cancel'}
    cancelHandler;
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for button handler
%%%%%%%%%%%%%%%%%%%%
function buttonHandler(event,incdec)

global gParams;

% get the value of the text field
val = get(gParams.ui.varentry(event),'string');

% if this is supposed to be a number, then make sure it is.
if ~strcmp(gParams.varinfo{event}.type,'string')
  % convert to number
  val = str2num(val);
  % check for incdec (this is for buttons that increment or decrement
  % the values, if one of these was passed, we will have by how much
  % we need to increment or decrement the value
  if isfield(gParams.varinfo{event},'incdec') && exist('incdec','var')
    val = val+incdec;
  end
  % check for minmax violations
  if isfield(gParams.varinfo{event},'minmax')
    if (val < gParams.varinfo{event}.minmax(1))
      disp(sprintf('(mrDefaultParamsGUI) Value %f lower than minimum %f',val,gParams.varinfo{event}.minmax(1)));
      val = [];
    elseif (val > gParams.varinfo{event}.minmax(2))
      disp(sprintf('(mrDefaultParamsGUI) Value %f greater than maximum %f',val,gParams.varinfo{event}.minmax(2)));
      val = [];
    end
  end
  % if the variable is empty, then set it back to default, this happens
  % if we got an invalid number entry
  if isempty(val)
    set(gParams.ui.varentry(event),'string',gParams.varinfo{event}.value);
  % otherwise remember this string as the default
  else
    gParams.varinfo{event}.value = num2str(val);
    set(gParams.ui.varentry(event),'string',gParams.varinfo{event}.value);
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for ok
%%%%%%%%%%%%%%%%%%%%
function helpHandler

global gParams;
if isfield(gParams,'helpFignum') && (gParams.helpFignum ~= -1)
  figure(gParams.helpFignum);
else
  gParams.helpFignum = figure;
end

% turn off menu/title etc.
set(gParams.helpFignum,'MenuBar','none');
set(gParams.helpFignum,'NumberTitle','off');
set(gParams.helpFignum,'Name','Parameter help');

% figure out how many rows
numrows = length(gParams.varinfo)+1;
numcols = 8;

% set the position and size
if ~isfield(gParams,'helpFigpos') || isempty(gParams.helpFigpos)
  figpos = get(gParams.helpFignum,'Position');
else
  figpos = gParams.helpFigpos;
end
figpos(4) = 2*gParams.topMargin+numrows*gParams.buttonHeight+(numrows-1)*gParams.margin;
figpos(3) = 2*gParams.leftMargin+numcols*gParams.buttonWidth+(numcols-1)*gParams.margin;
set(gParams.helpFignum,'Position',figpos);

% put up the info
for i = 1:length(gParams.varinfo)
  makeTextbox(gParams.varinfo{i}.name,i,1,1);
  set(makeTextbox(gParams.varinfo{i}.description,i,2,numcols-1),'HorizontalAlignment','Left');
end

% make close button
makeButton('Close','helpclose',numrows,numcols,1);

%%%%%%%%%%%%%%%%%%%%
% callback for helpclose
%%%%%%%%%%%%%%%%%%%%
function helpcloseHandler

global gParams;
if isfield(gParams,'helpFignum') && (gParams.helpFignum ~= -1)
  gParams.helpFigpos = get(gParams.helpFignum,'Position');
  close(gParams.helpFignum);
  gParams.helpFignum = -1;
else
  if ~isfield(gParams,'helpFigpos')
    gParams.helpFigpos = [];
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for ok
%%%%%%%%%%%%%%%%%%%%
function okHandler

global gParams;
gParams.ok = 1;
uiresume;

%%%%%%%%%%%%%%%%%%%%
% callback for cancel
%%%%%%%%%%%%%%%%%%%%
function cancelHandler

global gParams;
gParams.ok = 0;
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeButton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeButton(displayString,callback,rownum,colnum,uisize)

% make callback string
if isnumeric(callback)
  callback = sprintf('mrDefaultParamsGUI(%f)',callback);
else
  callback = sprintf('mrDefaultParamsGUI(''%s'')',callback);
end  

global gParams;

h = uicontrol('Style','pushbutton','Callback',callback,'String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(displayString,rownum,colnum,uisize)

global gParams;
h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname,'HorizontalAlignment','Right');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextentry(displayString,callback,rownum,colnum,uisize)

% make callback string
if isnumeric(callback)
  callback = sprintf('mrDefaultParamsGUI(%f)',callback);
else
  callback = sprintf('mrDefaultParamsGUI(''%s'')',callback);
end  

global gParams;

h = uicontrol('Style','edit','Callback',callback,'String',displayString,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makePopupmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makePopupmenu(displayString,rownum,colnum,uisize)

if ~iscell(displayString)
  choices{1} = displayString;
else
  choices = displayString;
end

global gParams;
h = uicontrol('Style','Popupmenu','Max',length(choices),'Min',1,'String',choices,'Value',1,'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeCheckbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeCheckbox(displayString,rownum,colnum,uisize)

global gParams;

h = uicontrol('Style','checkbox','Value',str2num(displayString),'Position',getUIControlPos(rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry w/inc dec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextentryWithIncdec(displayString,callback,rownum,colnum,uisize)

global gParams;

% make callback string
deccallback = sprintf('mrDefaultParamsGUI(%f,%f)',callback,gParams.varinfo{callback}.incdec(1));
inccallback = sprintf('mrDefaultParamsGUI(%f,%f)',callback,gParams.varinfo{callback}.incdec(2));

callback = sprintf('mrDefaultParamsGUI(%f)',callback);

% make inc and dec buttons
h = uicontrol('Style','pushbutton','Callback',deccallback,'String','<','Position',getUIControlPos(rownum,colnum,1),'FontSize',gParams.fontsize,'FontName',gParams.fontname);
h = uicontrol('Style','pushbutton','Callback',inccallback,'String','>','Position',getUIControlPos(rownum,colnum+2,1),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

% make text control
h = uicontrol('Style','edit','Callback',callback,'String',displayString,'Position',getUIControlPos(rownum,colnum+1,uisize-2),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(rownum,colnum,uisize)

% get global parameters
global gParams;

% get figure position
figpos = get(gParams.fignum,'Position');

% set this buttons width
thisButtonWidth = gParams.buttonWidth*uisize+(uisize-1)*gParams.margin;

% set the position for the button
pos(1) = gParams.margin + (gParams.buttonWidth+gParams.margin)*(colnum-1) + gParams.leftMargin;
pos(2) = figpos(4)-gParams.buttonHeight-gParams.topMargin - (gParams.buttonHeight+gParams.margin)*(rownum-1);
pos(3) = thisButtonWidth;
pos(4) = gParams.buttonHeight;
