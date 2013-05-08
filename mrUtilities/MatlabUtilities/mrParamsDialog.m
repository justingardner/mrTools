% mrParamsDialog.m
%
%      usage: mrParamsDialog(paramsInfo,<titleString>,<buttonWidth>,<callback>,<callbackArg>,<okCallback>,<cancelCallback>)
%  alt usage: You can also set variables explicitly, with the following syntax:
%             mrParamsDialog(paramsInfo,<titleString>,varargin)
%       e.g.: mrParamsDialog(paramsInfo,'This is the title','buttonWidth=1.5');
%             valid variable names are (buttonWidth,callback,callbackArg,okCallback,cancelCallback)
%             and also ignoreKeys (which keeps mrParamsDialog from allowing ESC to close it)
%         by: justin gardner
%       date: 03/13/07
%    purpose: creates a dialog for selection of parameters
%             see wiki for details
%
function [params params2] = mrParamsDialog(varargin)

% check arguments
if nargin < 1
  help mrParamsDialog
  return
end

% if this is a cell array, it means to open up the figure
% using the variable name, default value pairs given
if iscell(varargin{1})
  % if empty paramsInfo just return
  if isempty(varargin{1}),params = [];return,end
  % otherwise init the dialog
  [params params2] = initFigure(varargin{1},varargin);
  % otherwise it is a callback
else
  handleCallbacks(varargin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up figure in first palce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params params2] = initFigure(vars,otherParams)

global gParams;

% close any existing one
if isfield(gParams,'fignum')
  if ishandle(gParams.fignum)
    if isfield(gParams,'okCallback')
      feval(gParams.okCallback);
    end
    closeHandler;
  end
  global gParams;
end

global mrDEFAULTS;

% parse the input parameter string
[gParams.vars gParams.varinfo numrows numcols] = mrParamsParse(vars);

% get maximum length of var name
maxChars = 0;
for i = 1:length(gParams.vars)
  maxChars = max(length(gParams.vars{i}{1}),maxChars);
end  

% some basic info about location of controls
gParams.leftMargin = 10;
gParams.topMargin = 10;
gParams.buttonWidth = min(max(100,maxChars*7),200);
mver = matlabVersionNumber;
if strcmp(computer,'MACI') || (mver > 7.4)
  gParams.buttonHeight = 26;
else
  gParams.buttonHeight = 20;
end  
gParams.margin = 5;
gParams.fontsize = 12;
gParams.fontname = 'Helvetica';

% parse the otherParams. The first otherParams is always the title
if length(otherParams) > 1
  titleStr = otherParams{2};
  gParams.figlocstr = sprintf('mrParamsDialog_%s',fixBadChars(titleStr));
else
  titleStr = 'Set parameters';
  gParams.figlocstr = 'mrParamsDialog';
end

% now if there is a second otherParams and it is a string, then
% it means that we have been passed in a "getArgs" type. Otherwise
% it is the old calling convention in which the order determines what
% variable is being set
buttonWidth = [];callback = [];callbackArg = [];okCallback = [];cancelCallback = [];modal=[];
if (length(otherParams) > 2)
  if isstr(otherParams{3})
    getArgs(otherParams(3:end));
  else
    % get the arguments the old way, by order
    buttonWidth = otherParams{3};
    if length(otherParams) > 3,callback = otherParams{4}; end
    if length(otherParams) > 4,callbackArg = otherParams{5}; end
    if length(otherParams) > 5,okCallback = otherParams{6}; end
    if length(otherParams) > 6,cancelCallback = otherParams{7}; end
  end
end

% default to using a modal dialog if there is no callback
if isempty(modal)
  if ~isempty(callback)
    modal = 0;
  else
    modal = 1;
  end
end
% set the buttonWidth      
if ~isempty(buttonWidth)
  gParams.buttonWidth = gParams.buttonWidth*buttonWidth;
end

% get the figure
if ~isfield(gParams,'fignum') || (gParams.fignum == -1);
  % open figure, and turn off menu
  gParams.fignum = figure;
else
  figure(gParams.fignum);
end
set(gParams.fignum,'MenuBar','none');
set(gParams.fignum,'NumberTitle','off');
set(gParams.fignum,'Name',titleStr);
set(gParams.fignum,'closeRequestFcn','mrParamsDialog(''close'')');

% set height of figure according to how many rows we have
figpos = mrGetFigLoc(fixBadChars(gParams.figlocstr));
if isempty(figpos)
  figpos = get(gParams.fignum,'Position');
end
% if we have more than 25 rows then split into multiple columns
% but at most we make 6 multi columns
figMultiCols = min(ceil(numrows/25),6);
figrows = ceil(numrows/figMultiCols);
figcols = numcols*figMultiCols;
% for really big ones, reduce the button size
if (numcols > 2) && (figMultiCols > 3)
  gParams.buttonWidth = round(gParams.buttonWidth/2);
end
% set them in gParams
gParams.figrows = figrows;
gParams.figMultiCols = figMultiCols;
gParams.numcols = numcols;
gParams.numrows = numrows;
% set the figure position
figpos(4) = 2*gParams.topMargin+figrows*gParams.buttonHeight+(figrows-1)*gParams.margin;
figpos(3) = 2*gParams.leftMargin+figcols*gParams.buttonWidth+(figcols-1)*gParams.margin;
set(gParams.fignum,'Position',figpos);

% make entry buttons
rownum = 1;
for i = 1:length(gParams.varinfo)
  % make ui for varname
  gParams.ui.varname(i) = makeTextbox(gParams.fignum,gParams.varinfo{i}.name,rownum,1,1);
  % make ui entry dependent on what type we have
  if isfield(gParams.varinfo{i},'incdec')
    [gParams.ui.varentry{i} gParams.ui.incdec{i}(1) gParams.ui.incdec{i}(2)] =...
      makeTextentryWithIncdec(gParams.fignum,gParams.varinfo{i}.value,i,rownum,2,3);
    enableArrows(mrStr2num(gParams.varinfo{i}.value),i);
  elseif strcmp(gParams.varinfo{i}.type,'string')
    gParams.ui.varentry{i} = makeTextentry(gParams.fignum,gParams.varinfo{i}.value,i,rownum,2,3,gParams.varinfo{i}.editable);
  elseif strcmp(gParams.varinfo{i}.type,'checkbox')
    gParams.ui.varentry{i} = makeCheckbox(gParams.fignum,num2str(gParams.varinfo{i}.value),i,rownum,2,.5);
  elseif strcmp(gParams.varinfo{i}.type,'pushbutton')
    if isfield(gParams.varinfo{i},'buttonString')
      gParams.ui.varentry{i} = makeButton(gParams.fignum,gParams.varinfo{i}.buttonString,i,rownum,2,3);
    else
      gParams.ui.varentry{i} = makeButton(gParams.fignum,'',i,rownum,2,3);
    end
  elseif strcmp(gParams.varinfo{i}.type,'popupmenu') || iscell(gParams.varinfo{i}.value)
    gParams.ui.varentry{i} = makePopupmenu(gParams.fignum,gParams.varinfo{i}.value,i,rownum,2,3);
  elseif strcmp(gParams.varinfo{i}.type,'statictext')
    gParams.ui.varentry{i} = makeTextentry(gParams.fignum,gParams.varinfo{i}.value,i,rownum,2,3,0);
  elseif strcmp(gParams.varinfo{i}.type,'array')
    gParams.ui.varentry{i} = makeArrayentry(gParams.fignum,gParams.varinfo{i}.value,i,rownum,numcols,gParams.varinfo{i}.editable);
    rownum = rownum+size(gParams.varinfo{i}.value,1)-1;
  else
    gParams.ui.varentry{i} = makeTextentry(gParams.fignum,gParams.varinfo{i}.value,i,rownum,2,3,gParams.varinfo{i}.editable);
  end
  % check to see if we have to disable the entry field
  if isfield(gParams.varinfo{i},'enable') && isequal(gParams.varinfo{i}.enable,0)
    for j = 1:length(gParams.ui.varentry{i})
      set(gParams.ui.varentry{i}(j),'enable','off');
    end
  end
  rownum = rownum+1;
end

% set ok and cancel callback
if ~isempty(okCallback)
  gParams.okCallback = okCallback;
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
end
if ~isempty(cancelCallback)
  gParams.cancelCallback = cancelCallback;
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
end

% for each value that controls another one, call the buttonHandler to
% set up the correct dependency
for i = 1:length(gParams.varinfo)
  if isfield(gParams.varinfo{i},'controls')
    buttonHandler(i);
  end
end

gParams.callback = [];
% see if this has a callback, in which case we don't
% need to make ok/cancel buttons
if ~isempty(callback)
  gParams.callback = callback;
  % if another argument is specified that should
  % be sent as an argument to the callback function
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
  params = gParams.fignum;
  params2 = mrParamsGet(vars);
  % if another argument is specified than put up 
  % an ok button with the callback
  if ~isempty(okCallback)
    makeButton(gParams.fignum,'OK','ok',numrows,numcols,1);
  end
  % if a final argument is specified than put up 
  % an ok button with the callback
  if ~isempty(cancelCallback)
    makeButton(gParams.fignum,'Cancel','cancel',numrows,numcols-1,1);
  end
  makeButton(gParams.fignum,'Help','help',numrows,1,1);
  if ~modal,return,end
else
  gParams.callback = [];
end
% make ok and cancel buttons
if gParams.numcols > 2
  makeButton(gParams.fignum,'OK','ok',numrows,numcols,1);
  makeButton(gParams.fignum,'Cancel','cancel',numrows,numcols-1,1);
  makeButton(gParams.fignum,'Help','help',numrows,1,1);
else
  makeButton(gParams.fignum,'OK','ok',numrows,numcols+0.5,0.5);
  makeButton(gParams.fignum,'Cancel','cancel',numrows,numcols-0.1,0.5);
  makeButton(gParams.fignum,'Help','help',numrows,numcols-1,0.5);
end  

% set the input control to the first field that is editable
focusSet = 0;
if isfield(gParams,'ui') && isfield(gParams.ui,'varentry')
  % set the first editable field to have the keyboard focus
  for i = 1:length(gParams.varinfo)
    if gParams.varinfo{i}.editable
      % check to see if it is an array of handles
      if isequal(size(gParams.ui.varentry{i}),[1 1])
	H = gParams.ui.varentry{i};
      else
	H = gParams.ui.varentry{i}(1);
      end
      % confirm that we have a handle
      if ishandle(H)
	uicontrol(H);
	focusSet = 1;
	break;
      end
    end
  end
end

% if focus has not been set, then set the focus to the figure
% so the keyboard handler will let you Esc to cancel / enter to ok
if ~focusSet
  figure(gParams.fignum);
end

% set keyboard function
if ieNotDefined('ignoreKeys')
  set(gParams.fignum,'KeyPressFcn',@mrParamsKeyPressFcn);
end

% wait for user to hit ok or cancel (which sets uiresume)
uiwait;

% this can happen if we have opened up (and closed a second
% mrParamsDialog--this should be removed when this function
% is fixed to be able to run multiple simultaneous mrParamsDialogs;
if ieNotDefined('gParams'),params=[];params2=[];return,end

% check return value
if gParams.ok
  params = mrParamsGet(vars);
else
  % otherwise return empty
  params = [];
end
params2 = [];

closeHandler;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrParamsKeyPressFcn(figHandle,keyEvent)

switch (keyEvent.Key)
  case {'return'}
   okHandler;
  case {'escape'}
   closeHandler;
  case {'f1'}
   helpHandler;
end

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
    case {'close'}
      closeHandler;
    case {'helpclose'}
      helpcloseHandler;
    case {'cancel'}
      cancelHandler;
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for button handler
%%%%%%%%%%%%%%%%%%%%
function buttonHandler(varnum,incdec)

global gParams;

% if this is a push button then call it's callback
if strcmp(gParams.varinfo{varnum}.type,'pushbutton')
  if isfield(gParams.varinfo{varnum},'callback')
    args = {};getVars = 0;
    % if it wants optional arguments, pass that
    if isfield(gParams.varinfo{varnum},'callbackArg')
      args{end+1} = gParams.varinfo{varnum}.callbackArg;
    end
    % if the function wants the current parameter settings, pass that
    if isfield(gParams.varinfo{varnum},'passParams') && (gParams.varinfo{varnum}.passParams == 1)
      args{end+1} = mrParamsGet(gParams.vars);
    end
    % create the string to call the function
    funcall = 'gParams.varinfo{varnum}.value = feval(gParams.varinfo{varnum}.callback';
    for i = 1:length(args)
      funcall = sprintf('%s,args{%i}',funcall,i);
    end
    funcall = sprintf('%s);',funcall);
    % and call it
    eval(funcall);
  else
    disp(sprintf('(mrParamsDialog) Pushbutton %s does not have a callback',gParams.varinfo{varnum}.name));
  end
  return
end

% if this is supposed to be a number, then make sure it is.
if ~any(strcmp(gParams.varinfo{varnum}.type,{'string','array'}))
  if strcmp(gParams.varinfo{varnum}.type,'checkbox')
    val = get(gParams.ui.varentry{varnum},'Value');
  elseif strcmp(gParams.varinfo{varnum}.type,'popupmenu')
    val = [];
    if isfield(gParams.varinfo{varnum},'controls')
      % get the value from the list of values
      val = get(gParams.ui.varentry{varnum},'Value');
      val = gParams.varinfo{varnum}.value{val};
      if isstr(val),val=mrStr2num(val);end
    end
  else
    % get the value of the text field
    val = get(gParams.ui.varentry{varnum},'string');
    % convert to number
    val = mrStr2num(val);
  end
  % check for incdec (this is for buttons that increment or decrement
  % the values, if one of these was passed, we will have by how much
  % we need to increment or decrement the value
  if isfield(gParams.varinfo{varnum},'incdec') && exist('incdec','var')
    val = val+incdec;
  end
  % check if number needs to be round
  if isfield(gParams.varinfo{varnum},'round') && gParams.varinfo{varnum}.round
    val = round(val);
  end
  % check for minmax violations
  if isfield(gParams.varinfo{varnum},'minmax')
    if (val < gParams.varinfo{varnum}.minmax(1))
      disp(sprintf('(mrParamsDialog) Value %f lower than minimum %f',val,gParams.varinfo{varnum}.minmax(1)));
      val = [];
    elseif (val > gParams.varinfo{varnum}.minmax(2))
      disp(sprintf('(mrParamsDialog) Value %f greater than maximum %f',val,gParams.varinfo{varnum}.minmax(2)));
      val = [];
    end
  end
  % if the variable is empty, then set it back to default, this happens
  % if we got an invalid number entry
  if isempty(val)
    if ~strcmp(gParams.varinfo{varnum}.type,'popupmenu')
      set(gParams.ui.varentry{varnum},'string',gParams.varinfo{varnum}.value);
    end
    % otherwise remember this string as the default
  else
    if ~any(strcmp(gParams.varinfo{varnum}.type,{'popupmenu','checkbox'}))
      gParams.varinfo{varnum}.value = num2str(val);
      set(gParams.ui.varentry{varnum},'string',gParams.varinfo{varnum}.value);
    end
    % now check to see if this variable controls another one
    if isfield(gParams.varinfo{varnum},'controls')
      % go through all the controlled values
      for i = gParams.varinfo{varnum}.controls
        % enable or disable all the controlled fields
        if (val == 0)
          set(gParams.ui.varentry{i},'Enable','off');
        else
          set(gParams.ui.varentry{i},'Enable','on');
        end
        % now store the value currently being displayed
        if gParams.varinfo{i}.oldControlVal
          % get the current value
          currentValue = get(gParams.ui.varentry{i},'String');
          if strcmp(gParams.varinfo{i}.type,'popupmenu')
	    valueNum = get(gParams.ui.varentry{i},'Value');
	    % valueNum should not be 0, but sometimes it can get in that
	    % way if you click funny on the popupmenu and it doesn't select
	    % anything, this is a fix from getting lots of debug messages
	    % when that happens
	    if valueNum > 0
	      currentValue = putOnTopOfList(currentValue{valueNum},currentValue);
	    else
	      currentValue = putOnTopOfList(currentValue{1},currentValue);
	    end
          elseif strcmp(gParams.varinfo{i}.type,'checkbox')
            currentValue = get(gParams.ui.varentry{i},'Value');
          end
          % and save it
	  if strcmp(gParams.varinfo{i}.type,'array')
	    for k = 1:length(currentValue)
	      currentNumericValue(k) = str2num(currentValue{k});
	    end
	    currentNumericValue = reshape(currentNumericValue,size(currentValue));
	    gParams.varinfo{i}.allValues{gParams.varinfo{i}.oldControlVal}=currentNumericValue;
	  else
	    gParams.varinfo{i}.allValues{gParams.varinfo{i}.oldControlVal}=currentValue;
	  end
        end
        % switch to new value
        if (val >=1) && (val <= length(gParams.varinfo{i}.allValues))
          gParams.varinfo{i}.value = gParams.varinfo{i}.allValues{val};
	  % if this is an array, we have to set each individual array item
          if strcmp(gParams.varinfo{i}.type,'array')
	    for k = 1:length(gParams.varinfo{i}.allValues{val})
	      set(gParams.ui.varentry{i}(k),'String',gParams.varinfo{i}.allValues{val}(k));
	    end
	  else
	    if ~strcmp(gParams.varinfo{i}.type,'checkbox')
	      set(gParams.ui.varentry{i},'String',gParams.varinfo{i}.allValues{val});
	    end
	  end
	  % so more things to set for these types
	  if strcmp(gParams.varinfo{i}.type,'popupmenu')
            set(gParams.ui.varentry{i},'Value',1);
          elseif strcmp(gParams.varinfo{i}.type,'checkbox')
	    if isstr(gParams.varinfo{i}.value)
	      set(gParams.ui.varentry{i},'Value',str2num(gParams.varinfo{i}.value));
	    else
	      set(gParams.ui.varentry{i},'Value',gParams.varinfo{i}.value);
	    end
          end
          gParams.varinfo{i}.oldControlVal = val;
        end
      end
    end
  end
  % if the field has incdec, see how they should be grayed or not
  if isfield(gParams.varinfo{varnum},'incdec') && isfield(gParams.varinfo{varnum},'minmax')
    enableArrows(val,varnum)
  end
end
% update params
if isfield(gParams, 'callback')
  if ~isempty(gParams.callback)
    gParams.params = mrParamsGet(gParams.vars);
    if isfield(gParams,'callbackArg')
      feval(gParams.callback,gParams.params,gParams.callbackArg);
    else
      feval(gParams.callback,gParams.params);
    end
  end
end

% handle callbacks for non-push buttons
if ~ieNotDefined('gParams')
  if isfield(gParams.varinfo{varnum},'callback')
    if isfield(gParams.varinfo{varnum},'callbackArg')
      % create the string to call the function
      feval(gParams.varinfo{varnum}.callback,gParams.varinfo{varnum}.callbackArg,mrParamsGet(gParams.vars));
    else
      % create the string to call the function
      feval(gParams.varinfo{varnum}.callback,mrParamsGet(gParams.vars));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turn on or off incdec arrows depending on minmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableArrows(val,varnum)

global gParams;

if isfield(gParams.varinfo{varnum},'incdec') && isfield(gParams.varinfo{varnum},'minmax')
  if isnumeric(val)
    % turn on or off dec arrow
    if (val+gParams.varinfo{varnum}.incdec(1)) < gParams.varinfo{varnum}.minmax(1)
      set(gParams.ui.incdec{varnum}(1),'Enable','off');
    else
      set(gParams.ui.incdec{varnum}(1),'Enable','on');
    end
    % turn on or off inc arrow
    if (val+gParams.varinfo{varnum}.incdec(2)) > gParams.varinfo{varnum}.minmax(2)
      set(gParams.ui.incdec{varnum}(2),'Enable','off');
    else
      set(gParams.ui.incdec{varnum}(2),'Enable','on');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for help
%%%%%%%%%%%%%%%%%%%%
function helpHandler

global gParams;
global mrDEFAULTS;

% some settings that control how many rows/cols will be displayed
charsPerRow = 90;
numRowsPerCol = 30;

% figure out how many rows
numrows = 1;
% add number of rows each line needs
for i = 1:length(gParams.varinfo)
  numrows = numrows+max(1,ceil(length(gParams.varinfo{i}.description)/charsPerRow));
end
numcols = 5;

% see if we have to split into more columns
figMultiCols = ceil(numrows/numRowsPerCol);
figcols = numcols*figMultiCols;

% keep the settings
gParams.help.numcols = numcols;
gParams.help.numrows = numrows;
gParams.help.figrows = min(numrows,numRowsPerCol);
gParams.help.figMultiCols = min(figMultiCols,2);

% loop variable - controls which parameter we are displaying
startParam = 1;

iFig = 0;
for iCol = 1:2:figMultiCols

  % this figure num
  iFig = iFig + 1;

  % open the figure
  if isfield(gParams,'helpFignum') && (length(gParams.helpFignum)>=iFig) && (gParams.helpFignum(iFig) ~= -1)
    figure(gParams.helpFignum(iFig));
  else
    gParams.helpFignum(iFig) = figure;
  end

  % turn off menu/title etc.
  set(gParams.helpFignum(iFig),'MenuBar','none');
  set(gParams.helpFignum(iFig),'NumberTitle','off');
  set(gParams.helpFignum(iFig),'Name','Parameter help');

  % set close handler
  set(gParams.helpFignum(iFig),'DeleteFcn',@helpcloseHandler);

  % get number of columns for this figure
  thisMultiCols = min(figMultiCols-iCol+1,2);

  % set the position and size
  figpos = mrGetFigLoc('mrParamsDialogHelp');
  if isempty(figpos)
    figpos = get(gParams.helpFignum(iFig),'Position');
  end
  figpos(4) = 2*gParams.topMargin+gParams.help.figrows*gParams.buttonHeight+(gParams.help.figrows-1)*gParams.margin;
  figpos(3) = 2*gParams.leftMargin+thisMultiCols*numcols*gParams.buttonWidth+(thisMultiCols*numcols-1)*gParams.margin;
  figpos(2) = figpos(2) - 40*(iFig-1)/2;
  figpos(1) = figpos(1) + 40*(iFig-1)/2;
  set(gParams.helpFignum(iFig),'Position',figpos);
  
  % put up the info
  rownum = 1;
  for i = startParam:length(gParams.varinfo)
    % how many lines are needed for this line of text
    numLines = max(1,ceil(length(gParams.varinfo{i}.description)/charsPerRow));
    if (rownum+numLines) >= numRowsPerCol*2,break,end
    % make the textbox for the parameter name
    makeTextbox(gParams.helpFignum(iFig),gParams.varinfo{i}.name,rownum,1,1,numLines,1);
    % make the textbox for the help info
    set(makeTextbox(gParams.helpFignum(iFig),gParams.varinfo{i}.description,rownum,2,numcols-1,numLines,1),'HorizontalAlignment','Left');
    rownum = rownum+numLines;
  end
  startParam = i;
  
  % make close button
  makeButton(gParams.helpFignum(iFig),'Close','helpclose',min(figMultiCols-iCol+1,2)*gParams.help.figrows,numcols,1,1);
end

%%%%%%%%%%%%%%%%%%%%
% callback for close
%%%%%%%%%%%%%%%%%%%%
function closeHandler

global gParams;
if isempty(gParams),return,end
% close figure
mrSetFigLoc(fixBadChars(gParams.figlocstr),get(gParams.fignum,'Position'));
delete(gParams.fignum);

% close help
helpcloseHandler;

clear global gParams;
drawnow
% save figure locations .mrDefaults
saveMrDefaults;

%%%%%%%%%%%%%%%%%%%%
% callback for helpclose
%%%%%%%%%%%%%%%%%%%%
function helpcloseHandler(varargin)

global gParams;

if ~isfield(gParams,'helpFigpos')
  gParams.helpFigpos = [];
end

% close all help figures
if isfield(gParams,'helpFignum') 
  for iFig = length(gParams.helpFignum):-1:1
    if gParams.helpFignum(iFig) ~= -1
      mrSetFigLoc('mrParamsDialogHelp',get(gParams.helpFignum(iFig),'Position'));
      delete(gParams.helpFignum(iFig));
      gParams.helpFignum(iFig) = -1;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for ok
%%%%%%%%%%%%%%%%%%%%
function okHandler

global gParams;
gParams.ok = 1;
if isfield(gParams,'okCallback')
  if isfield(gParams,'callbackArg')
    feval(gParams.okCallback,gParams.callbackArg);
  else
    feval(gParams.okCallback);
  end
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%
% callback for cancel
%%%%%%%%%%%%%%%%%%%%
function cancelHandler

global gParams;
gParams.ok = 0;
if isfield(gParams,'cancelCallback')
  if isfield(gParams,'callbackArg')
    feval(gParams.cancelCallback,gParams.callbackArg);
  else
    feval(gParams.cancelCallback);
  end
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeButton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeButton(fignum,displayString,callback,rownum,colnum,uisize,isHelpDialog)

if ieNotDefined('isHelpDialog'),isHelpDialog=0;end
% make callback string
if isnumeric(callback)
  callback = sprintf('mrParamsDialog(%f)',callback);
else
  callback = sprintf('mrParamsDialog(''%s'')',callback);
end

global gParams;

h = uicontrol(fignum,'Style','pushbutton','Callback',callback,'String',displayString,'Position',getUIControlPos(fignum,rownum,colnum,uisize,[],isHelpDialog),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextbox makes an uneditable text box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(fignum,displayString,rownum,colnum,uisize,uisizev,isHelpDialog)

if ieNotDefined('isHelpDialog'),isHelpDialog=0;end
if ieNotDefined('uisizev'),uisizev=1;,end
global gParams;
h = uicontrol(fignum,'Style','text','String',displayString,'Position',getUIControlPos(fignum,rownum,colnum,uisize,uisizev,isHelpDialog),'FontSize',gParams.fontsize,'FontName',gParams.fontname,'HorizontalAlignment','Right');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextentry(fignum,displayString,callback,rownum,colnum,uisize,editable)

if ieNotDefined('editable'),editable=1;end

if editable
  style = 'edit';
else
  style = 'text';
end

% make callback string
if isnumeric(callback)
  callback = sprintf('mrParamsDialog(%f)',callback);
else
  callback = sprintf('mrParamsDialog(''%s'')',callback);
end

global gParams;

h = uicontrol(fignum,'Style',style,'Callback',callback,'String',displayString,'Position',getUIControlPos(fignum,rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeArrayentry makes a uicontrol to handle array entry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeArrayentry(fignum,array,callback,rownum,numcols,editable)

if ieNotDefined('editable'),editable=1;end

if editable
  style = 'edit';
else
  style = 'text';
end

% make callback string
if isnumeric(callback)
  callback = sprintf('mrParamsDialog(%f)',callback);
else
  callback = sprintf('mrParamsDialog(''%s'')',callback);
end

global gParams;

for i = 1:size(array,1)
  for j = 1:size(array,2)
    h(i,j) = uicontrol(fignum,'Style',style,'Callback',callback,'String',array(i,j),'Position',getUIControlPos(fignum,rownum+i-1,2+j-1,1),'FontSize',gParams.fontsize,'FontName',gParams.fontname);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makePopupmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makePopupmenu(fignum,displayString,callback,rownum,colnum,uisize)

callback = sprintf('mrParamsDialog(%f)',callback);

if ~iscell(displayString)
  choices{1} = displayString;
else
  if iscell(displayString{1})
    choices = displayString{1};
  else
    choices = displayString;
  end
end

global gParams;
h = uicontrol(fignum,'Style','Popupmenu','Callback',callback,'Max',length(choices),'Min',1,'String',choices,'Value',1,'Position',getUIControlPos(fignum,rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeCheckbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeCheckbox(fignum,displayString,callback,rownum,colnum,uisize)

global gParams;

% make callback string
callback = sprintf('mrParamsDialog(%f)',callback);

h = uicontrol(fignum,'Style','checkbox','Value',mrStr2num(displayString),'Callback',callback,'Position',getUIControlPos(fignum,rownum,colnum,uisize),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTextentry makes a uicontrol to handle text entry w/inc dec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h hl hr] = makeTextentryWithIncdec(fignum,displayString,callback,rownum,colnum,uisize)

global gParams;

% make callback string
deccallback = sprintf('mrParamsDialog(%f,%f)',callback,gParams.varinfo{callback}.incdec(1));
inccallback = sprintf('mrParamsDialog(%f,%f)',callback,gParams.varinfo{callback}.incdec(2));

callback = sprintf('mrParamsDialog(%f)',callback);

% make inc and dec buttons
hl = uicontrol(fignum,'Style','pushbutton','Callback',deccallback,'String','<','Position',getUIControlPos(fignum,rownum,colnum,1),'FontSize',gParams.fontsize,'FontName',gParams.fontname);
hr = uicontrol(fignum,'Style','pushbutton','Callback',inccallback,'String','>','Position',getUIControlPos(fignum,rownum,colnum+2,1),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

% make text control
h = uicontrol(fignum,'Style','edit','Callback',callback,'String',displayString,'Position',getUIControlPos(fignum,rownum,colnum+1,uisize-2),'FontSize',gParams.fontsize,'FontName',gParams.fontname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(fignum,rownum,colnum,uisize,uisizev,isHelpDialog)

if ieNotDefined('isHelpDialog'),isHelpDialog = 0;end

% get global parameters
global gParams;

% if we have too many parameters, we make them into two columns
% help and regular dialogs may have different numbers of rows/cols
if ~isHelpDialog
  figrows = gParams.figrows;
  numrows = gParams.numrows;
  numcols = gParams.numcols;
else
  figrows = gParams.help.figrows;
  numrows = gParams.help.numrows;
  numcols = gParams.help.numcols;
end
multiCol = ceil(rownum/figrows);
if multiCol > 1
  % always make sure that the last row end up on the
  % last row even if we have multiple columns
  if rownum == numrows
    rownum = figrows;
  else
    rownum = rownum-figrows*(multiCol-1);
  end
  colnum = colnum+numcols*(multiCol-1);
end

% get figure position
figpos = get(fignum,'Position');

% set this buttons width
thisButtonWidth = gParams.buttonWidth*uisize+(uisize-1)*gParams.margin;
if ieNotDefined('uisizev'),
  thisButtonHeight = gParams.buttonHeight;
else
  thisButtonHeight = gParams.buttonHeight*uisizev+gParams.margin*(uisizev-1);
end

% set the position for the button
pos(1) = gParams.margin + (gParams.buttonWidth+gParams.margin)*(colnum-1) + gParams.leftMargin;
pos(2) = figpos(4)-thisButtonHeight-gParams.topMargin - (gParams.buttonHeight+gParams.margin)*(rownum-1);
pos(3) = thisButtonWidth;
pos(4) = thisButtonHeight;

