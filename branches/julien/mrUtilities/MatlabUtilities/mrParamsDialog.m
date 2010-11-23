% mrParamsDialog.m
%
%      usage: mrParamsDialog(paramsInfo,<titleString>,<buttonWidth>,<callback>,<callbackArg>,<okCallback>,<cancelCallback>)
%  alt usage: You can also set variables explicitly, with the following syntax:
%             mrParamsDialog(paramsInfo,<titleString>,varargin)
%       e.g.: mrParamsDialog(paramsInfo,'This is the title','buttonWidth=1.5');
%             valid variable names are (buttonWidth,callback,callbackArg,okCallback,cancelCallback)
%             and also ignoreKeys (which keeps mrParamsDialog from allowing ESC to close it)
%         by: justin gardner, modified by julien besle
%       date: 03/13/07
%    purpose: creates a dialog for selection of parameters
%             see wiki for details
%        $Id$
%
function [params params2] = mrParamsDialog(varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7])
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
  
elseif isnumeric(varargin{1})% if it is 1 or 2 numbers then an entry field has been updated
  if length(varargin) == 1 
    buttonHandler(varargin{1});
  elseif length(varargin) == 4
    buttonHandler(varargin{1},varargin{2},varargin{3},varargin{4});
  end
  
else
  mrWarnDlg('(mrParamsDialog) unknown input parameter type');
end

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
% set up figure in first place
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params params2] = initFigure(vars,otherParams)

% close any existing params window
closeHandler;

global gParams; %gParams will only have non-graphical info about the parameters, including handles to uicontrols and figures
% dParams will be a local structure with graphical information specific to the dialog box
% genParms will be a local structure with general graphical information and will be cleared after the box is drawn

% parse the input parameter string
[gParams.vars gParams.varinfo numrows numcols] = mrParamsParse(vars);

% parse the otherParams. The first otherParams is always the title
if length(otherParams) > 1
  titleStr = otherParams{2};
else
  titleStr = '';
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
% button width is in fact a button scaling parameter    
if ~isempty(buttonWidth)
  genParams.buttonScale = buttonWidth;
end

% default to using a modal dialog if there is no callback
if isempty(modal)
  if ~isempty(callback)
    modal = 0;
  else
    modal = 1;
  end
end

% get the figure
if ~isfield(gParams,'fignum') || (gParams.fignum == -1);
  % open figure, and turn off menu
  gParams.fignum = figure;
else
  figure(gParams.fignum);
end
gParams.figlocstr{1} = sprintf('mrParamsDialog_%s',fixBadChars(titleStr));
set(gParams.fignum,'MenuBar','none');
set(gParams.fignum,'NumberTitle','off');
set(gParams.fignum,'closeRequestFcn',@closeHandler);
if isempty(titleStr)
  set(gParams.fignum,'Name','Set parameters');
else
  set(gParams.fignum,'Name',titleStr);
end

% some basic info about location of controls
genParams.minEntriesWidth = 200; %the minimum width of all the parameter entries
genParams.maxEntriesWidth = 500; %the maximum width of all the parameter entries
genParams.maxSingleFieldWidth = 100;
minVarNameWidth = 70;
genParams.margin = 5;
genParams.fontsize = 12;
genParams.fontname = 'Helvetica';
genParams.leftMargin = 10;
genParams.topMargin = 10;
%Matlab doesn't return the extent of multiline text, so we have to guess how much smaller the text height is 
%relative to the height of a textbox, in order to avoid making text boxes that are too large when text wraps
%(although there probably is a way to get this information)
genParams.lineHeightRatio = .67; %approximate height ot multiline text. 


% Collect information for uicontrol
%initialize varname field info(first column)
genParams.varNameWidth = 0;
genParams.varName = repmat({''},1,length(gParams.vars));
%initialize entry field info (second column)
dParams.entryWidth = -inf(1,length(gParams.vars));
dParams.entryValue = zeros(1,length(gParams.vars));
dParams.entryString = repmat({{''}},1,length(gParams.vars));
dParams.testString = repmat({''},1,length(gParams.vars));
dParams.entryStyle = repmat({''},1,length(gParams.vars));
dParams.entryNumCols = ones(1,length(gParams.vars));
dParams.entryNumRows = ones(1,length(gParams.vars));
dParams.incdec = zeros(length(gParams.vars),2);
dParams.numLines = ones(1,length(gParams.vars));
for i = 1:length(gParams.vars)
  if ~gParams.varinfo{i}.visible 
    dParams.numLines(i)=0; %no line for parameters that are not visible
  end
  %get variable name width
  genParams.varName{i} = [gParams.vars{i}{1} '  ']; %add spaces on the right
  h = uicontrol(gParams.fignum,'Style','text','String',genParams.varName{i},'FontSize',genParams.fontsize,'FontName',genParams.fontname);
  thisExtent = get(h,'extent');
  genParams.varNameWidth = max(minVarNameWidth,max(thisExtent(3)+2,genParams.varNameWidth));
  delete(h);

  %get info about the entries
  switch(gParams.varinfo{i}.type)
    case 'pushbutton' 
      dParams.entryStyle{i} = 'pushbutton';
      if isfield(gParams.varinfo{i},'buttonString')
        dParams.entryString{i} = {['  ' gParams.varinfo{i}.buttonString '  ']};%we need to allow some space for the button features
      end
      dParams.testString(i) = dParams.entryString{i};

    case 'popupmenu' 
      dParams.entryStyle{i} = 'popupmenu';
      dParams.entryValue(i) = 1;
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      %make up a string of Xs of lengh equal to the longest string in the menu list
      dParams.testString{i} =repmat('X',1,size(char(dParams.entryString{i}{1}),2));

    case 'statictext'
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      dParams.testString(i) = dParams.entryString{i};
      dParams.entryStyle{i} = 'text';

    case 'checkbox'
      if isnumeric(gParams.varinfo{i}.value)
        dParams.entryValue(i) = gParams.varinfo{i}.value;
      else
        dParams.entryValue(i) = str2num(gParams.varinfo{i}.value);
      end
      dParams.entryStyle{i} = 'checkbox';

    case 'string'
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
        dParams.testString(i) = dParams.entryString{i};
      else
        dParams.entryStyle{i} = 'edit';
      end

    case 'numeric'
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
        dParams.testString(i) = dParams.entryString{i};
      else
        dParams.entryStyle{i} = 'edit';
      end

    case 'array'
      dParams.entryString{i} = num2cell(gParams.varinfo{i}.value);
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
      else
        dParams.entryStyle{i} = 'edit';
      end
      dParams.entryNumCols(i) = size(dParams.entryString{i},2);
      dParams.entryNumRows(i) = size(dParams.entryString{i},1);

    case 'stringarray'
      dParams.entryString{i} = gParams.varinfo{i}.value;
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
      else
        dParams.entryStyle{i} = 'edit';
      end
      dParams.entryNumCols(i) = size(dParams.entryString{i},2);
      dParams.entryNumRows(i) = size(dParams.entryString{i},1);

    otherwise
       keyboard %unknown type...
%         dParams.entryString{i} = gParams.varinfo{i}.value;
%         dParams.entryStyle{i} = gParams.varinfo{i}.type;
  end


  if isfield(gParams.varinfo{i},'incdec')
    dParams.incdec(i,:)=gParams.varinfo{i}.incdec;
  end

end  

%optimize figure dimensions
[figpos,dParams,genParams] = optimizeFigure(gParams.fignum,gParams.figlocstr{1},dParams,genParams);
figWidth = figpos(3);
figHeight = figpos(4);

%cap widths that are more than the max
dParams.entryWidth(dParams.entryWidth>dParams.allEntriesWidth)=dParams.allEntriesWidth;
%if it's gonna be an array compute the field width 
for i = 1:length(dParams.entryStyle)
  if dParams.entryNumCols(i)>1 || dParams.entryNumRows(i)>1
      dParams.entryWidth(i) = min(dParams.allEntriesWidth/dParams.entryNumCols(i)-genParams.margin,genParams.maxSingleFieldWidth);
  end
end

%set control dimensions to normalized so that the figure resizes
set(gParams.fignum,'defaultUicontrolUnits','normalized'); 
%computing the normalized positions is taken care of by getUIControlPos


%--------------------------------- make entry buttons-----------------------------------
for i = 1:length(gParams.varinfo)
  % make ui for varname
  gParams.ui.varname(i) = makeUIcontrol(i,gParams.fignum,dParams,genParams,'varname');
  % make ui for entry
  [gParams.ui.varentry{i} gParams.ui.incdec{i}{1} gParams.ui.incdec{i}{2}] =...
     makeUIcontrol(i,gParams.fignum,dParams,genParams,'varentry');
  if isfield(gParams.varinfo{i},'incdec')
    enableArrows(mrStr2num(gParams.varinfo{i}.value),i);
  end
  % check enable/visible options
  if isfield(gParams.varinfo{i},'enable') && isequal(gParams.varinfo{i}.enable,0)
      set(gParams.ui.varentry{i},'enable','off');
  end
  if isfield(gParams.varinfo{i},'visible') && isequal(gParams.varinfo{i}.visible,0)
    set(gParams.ui.varentry{i},'visible','off');
    set(gParams.ui.varname(i),'visible','off');
  end
end

% for each value that controls another one, call the buttonHandler to
% set up the correct dependency
for i = 1:length(gParams.varinfo)
  if isfield(gParams.varinfo{i},'controls')
    buttonHandler(i);
  end
end

%--------------------------------- make Help/Ok/Cancel buttons-----------------------------------
gParams.callback = [];
makeOkButton = 1;
makeCancelButton = 1;
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
  % if another argument is specified then put up 
  % an ok button with the callback
  if ~isempty(okCallback)
    gParams.okCallback = okCallback;
  else
    makeOkButton = 0;
  end
  % if a final argument is specified then put up 
  % an ok button with the callback
  if ~isempty(cancelCallback)
    gParams.cancelCallback = cancelCallback;
  else
    makeCancelButton = 0;
  end
else
  gParams.callback = [];
end

% position of buttons
totalColWidth = 1/dParams.multiCols;
thisButtonWidth = min(100/figWidth,totalColWidth/3);
bottomMargin = genParams.topMargin/figHeight;
thisButtonHeight = genParams.buttonHeight/figHeight;
intervalBetweenButtons = (totalColWidth - thisButtonWidth*3)/(4);
leftPosition = (dParams.multiCols - 1)/dParams.multiCols + intervalBetweenButtons;

%Help Button
gParams.fignum(2) = figure('visible','off');
set(gParams.fignum(2),'userdata',0); %this is just to tell helpHandler if the help figure has been drawn or not
gParams.figlocstr{2} = sprintf('mrParamsDialogHelp_%s',fixBadChars(titleStr));
gParams.helpButton = uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',{@helpHandler,gParams.fignum(2),genParams},'String','Show help',...
  'Position',[leftPosition bottomMargin thisButtonWidth thisButtonHeight],...
  'FontSize',genParams.fontsize,'FontName',genParams.fontname);

%Cancel Button
if makeCancelButton
  uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',@cancelHandler,'String','Cancel',...
  'Position',[leftPosition+(intervalBetweenButtons+thisButtonWidth) bottomMargin thisButtonWidth thisButtonHeight],...
  'FontSize',genParams.fontsize,'FontName',genParams.fontname);
end

%Ok Button
if makeOkButton
  uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',@okHandler,'String','OK',...
    'Position',[leftPosition+(intervalBetweenButtons+thisButtonWidth)*2 bottomMargin thisButtonWidth thisButtonHeight],...
    'FontSize',genParams.fontsize,'FontName',genParams.fontname);
end

%if non-modal, quit here
if ~modal,return,end


% set the input control to the first field that is editable
focusSet = 0;
if isfield(gParams,'ui') && isfield(gParams.ui,'singleEntry')
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
  figure(gParams.fignum(1));
end

% set keyboard function
if ieNotDefined('ignoreKeys')
  set(gParams.fignum(1),'KeyPressFcn',@mrParamsKeyPressFcn);
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


%%%%%%%%%%%%%%%%%%%%
% callback for button handler
%%%%%%%%%%%%%%%%%%%%
function buttonHandler(varnum,incdec,entryRow,entryCol)

global gParams;

% if this is a push button then call its callback
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
if ~any(strcmp(gParams.varinfo{varnum}.type,{'string','array','stringarray'}))
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
      set(gParams.ui.incdec{varnum}{1},'Enable','off');
    else
      set(gParams.ui.incdec{varnum}{1},'Enable','on');
    end
    % turn on or off inc arrow
    if (val+gParams.varinfo{varnum}.incdec(2)) > gParams.varinfo{varnum}.minmax(2)
      set(gParams.ui.incdec{varnum}{2},'Enable','off');
    else
      set(gParams.ui.incdec{varnum}{2},'Enable','on');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for help
%%%%%%%%%%%%%%%%%%%%
function helpHandler(handle,event,fignum,genParams)

if strcmp(get(fignum,'visible'),'on')
  set(fignum,'visible','off');
  set(handle,'string','Show Help');
else
  set(fignum,'visible','on');
  set(handle,'string','Hide Help');
end

if get(fignum,'userdata')

else  
  global gParams

  % turn off menu/title etc.
  set(fignum,'MenuBar','none');
  set(fignum,'NumberTitle','off');
  set(fignum,'Name','Parameter help');
  set(fignum,'closeRequestFcn',@helpcloseHandler);

  set(fignum,'defaultUicontrolUnits','pixels'); 
  
  dParams.entryWidth = zeros(1,length(gParams.varinfo));
  dParams.entryValue = zeros(1,length(gParams.varinfo));
  dParams.entryString = repmat({{''}},1,length(gParams.varinfo));
  dParams.entryStyle =  repmat({'text'},1,length(gParams.varinfo));
  dParams.incdec = zeros(length(gParams.vars),2);
  dParams.testString = repmat({''},1,length(gParams.varinfo));
  dParams.entryNumCols = ones(1,length(gParams.varinfo));
  dParams.entryNumRows = ones(1,length(gParams.varinfo));
  dParams.numLines = ones(1,length(gParams.varinfo));
  for i = 1:length(gParams.varinfo)
    if ~gParams.varinfo{i}.visible %no need to display the help is parameter is not visible
      dParams.numLines(i)=0;
    end
    dParams.entryString{i} = {[' ' gParams.varinfo{i}.description]}; %add 1 space on the left
    dParams.testString(i) = dParams.entryString{i};
  end

  %gParams.fignum(2) = fignum;
  
  %compute figure dimensions based on number of rows and colums
  [figpos,dParams, genParams] = optimizeFigure(fignum,gParams.figlocstr{2},dParams,genParams);
  
  %set the all the entry widths to the max 
  dParams.entryWidth(:)=dParams.allEntriesWidth;
  
  %set control dimensions to normalized so that the figure resizes
  set(fignum,'defaultUicontrolUnits','normalized'); 
  % put up the info
  for i = 1:length(gParams.varinfo)
    if gParams.varinfo{i}.visible %no need to display the help is parameter is not visible
      makeUIcontrol(i,gParams.fignum(2),dParams,genParams,'varname');
      set(makeUIcontrol(i,gParams.fignum(2),dParams,genParams,'varentry'),'HorizontalAlignment','Left');
    end
  end

  % make close button
  uicontrol(gParams.fignum(2),'Style','pushbutton','Callback',@helpcloseHandler,'String','Close',...
    'Position',getUIControlPos(fignum,dParams,genParams,dParams.numrows,1,2,dParams.allEntriesWidth,1),...
    'FontSize',genParams.fontsize,'FontName',genParams.fontname);
  
  set(fignum,'userdata',1)
end

%%%%%%%%%%%%%%%%%%%%
% callback for helpclose
%%%%%%%%%%%%%%%%%%%%
function helpcloseHandler(varargin)

global gParams;

set(gParams.fignum(2),'visible','off')
set(gParams.helpButton,'string','Show Help');

%%%%%%%%%%%%%%%%%%%%
% callback for close
%%%%%%%%%%%%%%%%%%%%
function closeHandler(varargin)

global gParams;
if isempty(gParams),return,end

if isfield(gParams,'fignum') 
  if isfield(gParams,'figlocstr')
  % save figure locations .mrDefaults
    for iFig = 1:length(gParams.fignum)
      mrSetFigLoc(fixBadChars(gParams.figlocstr{iFig}),get(gParams.fignum(iFig),'Position'));
    end
  end
  % close figure
  delete(gParams.fignum);
end
saveMrDefaults;

clear global gParams;
drawnow


%%%%%%%%%%%%%%%%%%%%
% callback for ok
%%%%%%%%%%%%%%%%%%%%
function okHandler(varargin)

global gParams;
gParams.ok = 1;
if isfield(gParams,'okCallback')
  feval(gParams.okCallback);
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%
% callback for cancel
%%%%%%%%%%%%%%%%%%%%
function cancelHandler(varargin)

global gParams;
gParams.ok = 0;
if isfield(gParams,'cancelCallback')
  feval(gParams.cancelCallback);
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeUIcontrol makes an uicontrol of any type %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hEntry,hMinus,hPlus] = makeUIcontrol(varnum,fignum,dParams,genParams,columnType)

hMinus = [];
hPlus = [];

switch(columnType)
  case 'varname'
    colnum=1;
    entryWidth = genParams.varNameWidth;
    fieldnums = 1;
    rownums = sum(dParams.numLines(1:varnum-1).*dParams.entryNumRows(1:varnum-1))+1;
    numLines = 	dParams.numLines(varnum)*dParams.entryNumRows(varnum);
    entryString = genParams.varName(varnum);
    style = 'text'; 
    hAlignment = 'right';
    incdec = [0 0];
    
  case 'varentry'
    colnum =2;
    entryWidth = dParams.entryWidth(varnum);
    fieldnums =1:dParams.entryNumCols(varnum);
    rownums = sum(dParams.numLines(1:varnum-1).*dParams.entryNumRows(1:varnum-1))+(1:dParams.entryNumRows(varnum));
    numLines = 	dParams.numLines(varnum);
    entryString = dParams.entryString{varnum};
    style = dParams.entryStyle{varnum};
    hAlignment = 'center';
    incdec = dParams.incdec(varnum,:);
    
end
if ~numLines %if numLines==0, that means the control is invisble
  numLines =1;  %set it to 1 to avoid an error
end

for i=1:length(rownums)
  for j=fieldnums
    hEntry(i,j) = uicontrol(fignum,...
    'Style',style,...
    'Callback',sprintf('mrParamsDialog(%f)',varnum),...  %callback has no effect if textbox
    'String',entryString{i,j},...
    'Value',dParams.entryValue(varnum),...
    'Position',getUIControlPos(fignum,dParams,genParams,rownums(i),numLines,colnum,entryWidth,j),...
    'HorizontalAlignment',hAlignment,...
    'FontSize',genParams.fontsize,'FontName',genParams.fontname);
    if any(incdec)
      % make callback string
      deccallback = sprintf('mrParamsDialog(%f,%f,%f,%f)',varnum,incdec(1),i,j);
      inccallback = sprintf('mrParamsDialog(%f,%f,%f,%f)',varnum,incdec(2),i,j);
      hMinus(i,j) = uicontrol(fignum,'Style','pushbutton','Callback',deccallback,'String','-',...
        'Position',getUIControlPos(fignum,dParams,genParams,rownums(i)+.2,numLines*.75,colnum,20,j+.05),... %these dimensions are totally empirical
        'FontSize',genParams.fontsize,'FontName',genParams.fontname);
      hPlus(i,j) = uicontrol(fignum,'Style','pushbutton','Callback',inccallback,'String','+',...
        'Position',getUIControlPos(fignum,dParams,genParams,rownums(i)+.02,numLines*.6,colnum,20,j+.05),... %these dimensions are totally empirical
        'FontSize',genParams.fontsize,'FontName',genParams.fontname);
    end
  end
end

if colnum==2 && strcmp(dParams.entryStyle{varnum},'checkbox') %on windows, the backgroud of checkboxes is colored
  set(hEntry,'BackgroundColor',get(gcf,'color')); %even without string, which is ugly
end
if strcmp(dParams.entryStyle{varnum},'text') && isempty(entryString{1}) %if it's an empty textbox
  set(hEntry,'visible','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol %
%   dealing with multicolumns and margins            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(fignum,dParams,genParams,rownum,numLines,colnum,entryWidth,fieldNum)

numcols = 2;
figrows = dParams.figrows;
numrows = dParams.numrows;
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

% set the horizontal position and width for the button
if colnum - (multiCol-1)*numcols == 1
  pos(1) = genParams.leftMargin + ... %position
            (multiCol - 1) * (genParams.varNameWidth + genParams.margin) + ...
            (multiCol - 1) * (dParams.allEntriesWidth+genParams.margin) + ...
            genParams.margin; 
else
  pos(1) = genParams.leftMargin + ...
    multiCol*(genParams.margin + genParams.varNameWidth) + ...
    (colnum-multiCol-1) * (dParams.allEntriesWidth+genParams.margin) + ...
    (fieldNum-1)*(entryWidth+genParams.margin)+...
    genParams.margin;
end
pos(3) = entryWidth; 

% set the vertical position and height for the button
pos(4) = genParams.buttonHeight*numLines+genParams.margin*(numLines-1);
pos(2) = figpos(4)-pos(4)-genParams.topMargin - (genParams.buttonHeight+genParams.margin)*(rownum-1);

%normalize position
pos([1 3]) = pos([1 3])/figpos(3);
pos([2 4]) = pos([2 4])/figpos(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizeFigure optimizes the number of rows and columns as well as the dimensions of the figure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [figpos,dParams,genParams] = optimizeFigure(fignum,figLocStr,dParams,genParams)

%compute figure dimensions based on number of rows and colums
figpos = mrGetFigLoc(fixBadChars(figLocStr));
if isempty(figpos)
  figpos = get(fignum,'Position');
end


for i = 1:length(dParams.testString)  
  if ~isempty(dParams.testString{i}) && dParams.numLines(i)~=0
    %compute number of lines using string width if it's gonna be displayed using a text box, a popupmenu or a pushbutton
    h = uicontrol(fignum,'Style',dParams.entryStyle{i},'String',dParams.testString{i},'FontSize',genParams.fontsize,'FontName',genParams.fontname);
    thisExtent = get(h,'extent');
    dParams.entryWidth(i) = thisExtent(3)+20; %we need to allow some space for the button features
    delete(h);
  end
end
if ieNotDefined('thisExtent')
  h = uicontrol(fignum,'Style','Text','String','X','FontSize',genParams.fontsize,'FontName',genParams.fontname);
  thisExtent = get(h,'extent');
  delete(h);
end
genParams.buttonHeight = thisExtent(4);

%For edit boxes and buttons on MACs, this height will be too small because of their large borders
if strcmp(computer,'MACI') || strcmp(computer,'MACI64') 
  genParams.buttonHeight = genParams.buttonHeight*1.3;
end
% global mrDEFAULTS;
% mver = ver('matlab');mver = str2num(mver.Version);
% if strcmp(computer,'MACI') || strcmp(computer,'MACI64') || (mver > 7.4)
%   genParams.buttonHeight = 26;
% else
%   genParams.buttonHeight = 22;
% end  


maxEntryNumCols = max(dParams.entryNumCols);
%the total field width is whatever field has the largest width, within the min and max parameters
dParams.allEntriesWidth = max(max(dParams.entryWidth),min(maxEntryNumCols*genParams.maxSingleFieldWidth,genParams.maxEntriesWidth));
dParams.allEntriesWidth = max(genParams.minEntriesWidth,min(genParams.maxEntriesWidth,dParams.allEntriesWidth));
%add space for margins
dParams.allEntriesWidth = dParams.allEntriesWidth + 2*genParams.margin;
%replace non-set widths by the max width
dParams.entryWidth(dParams.entryWidth<0)= dParams.allEntriesWidth;


% get  number of lines for fields that might wrap
for i = 1:length(dParams.entryStyle)  
  if ~isempty(dParams.testString{i}) && dParams.numLines(i)~=0
    dParams.numLines(i) = ceil(ceil(dParams.entryWidth(i)/dParams.allEntriesWidth)*genParams.lineHeightRatio);
  end
end

%numLines is an array of line number per parameter
numrows = sum(dParams.numLines.*dParams.entryNumRows)+1; %we add one for the help/ok/cancel buttons

%start with one multicolumn and set other parameters accordingly
dParams.multiCols=1;
dParams.figrows = numrows;
figHeight = 2*genParams.topMargin+dParams.figrows*genParams.buttonHeight+(dParams.figrows-1)*genParams.margin;
figWidth = 2*genParams.leftMargin+genParams.margin+dParams.multiCols*(genParams.varNameWidth+dParams.allEntriesWidth+genParams.margin);

%optimize figure dimensions 
screenSize = get(0,'MonitorPositions');
thresholdRatio = 1.7; 
%while one of the dimensions is larger than the screen or the height/width is over the threshold, resize
while figHeight/figWidth>thresholdRatio || figHeight>screenSize(1,4) || figWidth>screenSize(1,3)
  %if height/width > thresholdRatio or if height > screenheight, we add a column
  if figHeight/figWidth>thresholdRatio || figHeight>screenSize(1,4)
    dParams.multiCols = dParams.multiCols+1;
    %compute new number of rows per columns, but make sure we're not cutting a field
    dParams.figrows = ceil(numrows/dParams.multiCols);
    while ~ismember(dParams.figrows,cumsum(dParams.numLines.*dParams.entryNumRows))
      dParams.figrows = dParams.figrows+1;
    end
  elseif figWidth > screenSize(1,3) %else if width>screen width, we reduce the button width
    dParams.allEntriesWidth = (screenSize(1,3)-dParams.multiCols*(genParams.varNameWidth+genParams.margin) -2*genParams.leftMargin) / dParams.multiCols;
  end
  
  %compute the new dimensions
  figHeight = 2*genParams.topMargin+dParams.figrows*genParams.buttonHeight+(dParams.figrows-1)*genParams.margin;
  figWidth = 2*genParams.leftMargin+genParams.margin+dParams.multiCols*(genParams.varNameWidth+dParams.allEntriesWidth+genParams.margin);
end
dParams.numrows = sum(dParams.numLines.*dParams.entryNumRows)+1; %this is the total numbe of rows including one for help/ok/cancel buttons

% set the figure position
figpos(4) = figHeight;
figpos(3) = figWidth;
%make sure the figure is not outside the screen
figpos(1) = min(figpos(1),sum(screenSize([1 3]))-1-figWidth);
figpos(2) = min(figpos(2),sum(screenSize([2 4]))-1-figHeight);

set(fignum,'Position',figpos);




