% mlrAdjustGUI.m
%
%      usage: mlrAdjustGUI(v,command,varargin)
%         by: justin gardner
%       date: 10/31/10
%    purpose: Adjusts the MLR GUI to allow plug-in code
%             to change the behavior of the GUI. Note
%             that this serves a different function then
%             mlrGuiSet which is used to set user interface
%             items during normal operations (like graying
%             out a menu item, or changing the value of a control).
%             This is intended to be used only when MLR is loaded
%             to adjust the GUI to suit the needs of different sites.
%
%             Every UI interface item on the MLR window
%             can be specified by either a tag (a string that
%             identifies it), or by the label for a menu item.
%             e.g., the menu item under ROI called Show
%             is /ROI/Show. or e.g., the cortical depth slider
%             has the tag corticalDepth. These identifiers are
%             what is called the "controlName" below and you
%             can get a full list of all possible controlNames 
%             by doing:
%             mlrAdjustGUI(getMLRView,'list');
%
%             To set a property of a control: 
%             mlrAdjustGUI(v,'set',controlName,propertyName,propertyValue);
%      e.g.:  mlrAdjustGUI(getMLRView,'set','baseGammaSlider','Visible','off');
% 
%             Similarly, to set the callback for a control
%      e.g.:  mlrAdjustGUI(getMLRView,'set','baseGammaSlider','Callback',@testCallback);
%             Where testCallback is a function that takes two arguments and
%             usually begins with a snippet of code that gets the view:
%             function testCallback(hObject,eventdata)
%             v = viewGet(getfield(guidata(hObject),'viewNum'),'view');
%
%             To add a new menu item:
%             mlrAdjustGUI(v,'add','menu',menuName,menuLocation,propertyName1,propertyValue1,...);
%      e.g.:  mlrAdjustGUI(getMLRView,'add','menu','Plugin','Plots','Callback',@testCallback,'Separator','on');
%             The above will add the menu item Plugin after the menu identified
%             as Plots. If you wanted instead to put it at the top
%             of the Plots menu, then set the menuLocation to /Plots/
%             
%             To add an interrogator function as a default one
%             (that shows up in the GUI)
%             mlrAdjustGUI(v,'add','interrogator',interrogatorName)
%      e.g.:  mlrAdjustGUI(getMLRView,'add','interrogator','eventRelatedPlot');
%
%             To add colormap functions which will show up in
%             /Edit/Overlay
%             mlrAdjustGUI(v,'add','colormap',colormapName)
%       e.g.: mlrAdjustGUI(getMLRView,'add','colormap','gray');
%
function retval = mlrAdjustGUI(v,command,varargin)

% test code

% default return empty
retval = [];

% check arguments
if nargin < 2
  help mlrAdjustGUI
  return
end

% get figure
if isempty(v),disp(sprintf('(mlrAdjustGUI) Empty view. Is MLR closed?'));return,end
f = viewGet(v,'fignum');
if isempty(f),disp(sprintf('(mlrAdjustGUI) Passed in view does not have a figure associated with it'));return;end

% get controls and menus
uiControls = getUiControls(f);
menuControls = getMenuControls(f);

% do the commandedaction
switch command
 case {'set'}
  setItemProperty(varargin,uiControls,menuControls)
 case {'list'}
  listControlNames(uiControls,menuControls);
 case {'add'}
  switch varargin{1}
    case {'menu'}
     addMenu({varargin{2:end}},menuControls);
    case {'interrogator','interrogators'}
     addInterrogator(v,varargin{2});
    case {'colormap','colormaps'}
     addColormap(v,varargin{2});
    case {'control'}
     addControl(f,{varargin{2:end}},uiControls);
    otherwise
      mrWarnDlg(['(mlrAdjustGUI) Unknow object type ' varargin{1}]);
  end
 case {'remove'}
  switch varargin{1}
    case {'menu'}
     removeMenu(varargin(2),menuControls);
  end
 case {'get'}
  % return handle, check menu items and ui controls
  retval = getHandle(varargin{1},menuControls,uiControls);
 otherwise
  disp(sprintf('(mlrAdjustGUI) Unknown command: %s',command));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  listControlNames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function listControlNames(uiControls,menuControls)

% list uiControl tags
disp('============== UI Controls ============== ');
for i = 1:length(uiControls)
  disp(sprintf('(%i) tag: %s Style: %s',i,uiControls(i).tag,get(uiControls(i).h,'Style')));
end

% list menu labels and tags
disp('============== Menu Items ============== ');
for i = 1:length(menuControls)
  disp(sprintf('(%i) Menu label: ''%s'' tag: %s',i,menuControls(i).fullLabel,menuControls(i).tag));
end


%%%%%%%%%%%%%%%%%%%
%%%   getHandle  %%
%%%%%%%%%%%%%%%%%%%
function h = getHandle(itemName,controls,controls2)

h = [];

% check if there is a match in the following fields
searchFields = {'tag','label','fullLabel'};

for i = 1:length(searchFields)
  % check for existence of field
  if isfield(controls,searchFields{i})
    % and if there is a match
    [tf itemNum] = find(strcmp(itemName,{controls.(searchFields{i})}));
    if tf
      h = controls(itemNum).h;
      return
    end
  end
end
  
% if we didn't find anything and we were passed two
% sets of controls, then check second set
if isempty(h) && (nargin ==3)
  h = getHandle(itemName,controls2);
end

% if we didn't find anything and the thing ends with '/' then
% try without the '/'
if isempty(h) && (length(itemName)>0) && (itemName(end) == '/')
  if nargin == 3
    h = getHandle(itemName(1:end-1),controls,controls2);
  else
    h = getHandle(itemName(1:end-1),controls);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  addInterrogator   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function addInterrogator(v,interrogatorList)

% make into a cell array and set in view
interrogatorList = cellArray(interrogatorList);
viewSet(v,'defaultInterrogators',interrogatorList);

% display what we are doing
disp(sprintf('(mlrAdjustGUI) Adding default interrogators: %s',cellToCommaDelimited(interrogatorList)));

%%%%%%%%%%%%%%%%%%%%%
%%%  addColormap   %%
%%%%%%%%%%%%%%%%%%%%%
function addColormap(v,colormapList)

% make into a cell array and set in view
colormapList = cellArray(colormapList);
viewSet(v,'colormaps',colormapList);

% display what we are doing
disp(sprintf('(mlrAdjustGUI) Adding colormaps: %s',cellToCommaDelimited(colormapList)));

%%%%%%%%%%%%%%%%%%%%%
%%%   addControl  %%%
%%%%%%%%%%%%%%%%%%%%%
function addControl(f,args,uiControls)

% check length of arguments
if length(args) < 1
  disp(sprintf('(mlrAdjustGUI:addControl) Requires at least arguments: controlName'));
  return
else
  % name the arguments
  controlName = args{1};
  controlProperties = {args{2:end}};
  if isodd(length(controlProperties) )
    disp(sprintf('(mlrAdjustGUI:addControl) Properties must all have a matching property value'));
    return
  end
end

% check to see if it has already been added
if ~isempty(getHandle(controlName,uiControls))
  disp(sprintf('(mlrAdjustGUI:addControl) Already added menu item: %s',controlName));
  return
end

% get the gui data
h = guidata(f);

h.(controlName)=uicontrol(f);
set(h.(controlName),'unit','normalized');
% add all the properties
for i = 1:2:length(controlProperties)
  set(h.(controlName),controlProperties{i},controlProperties{i+1});
end

guidata(f,h);

%%%%%%%%%%%%%%%%%%
%%%   addMenu  %%%
%%%%%%%%%%%%%%%%%%
function addMenu(args,menuControls)

% check length of arguments
if length(args) < 2
  disp(sprintf('(mlrAdjustGUI:addMenu) Requires 2 arguments: menuName, menuLocation'));
  return
else
  % name the arguments
  menuName = args{1};
  menuLocation = args{2};
  menuProperties = {args{3:end}};
  if isodd(length(menuProperties) )
    disp(sprintf('(mlrAdjustGUI:addMenu) Properties must all have a matching property value'));
    return
  end
end

% go look for the item location
h = getHandle(menuLocation,menuControls);
% if not found, then print warning, return
if isempty(h)
  disp(sprintf('(mlrAdjustGUI:addMenu) Could not find menu location: %s',menuLocation));
  return
end

% check to see if it has already been added
if ~isempty(getHandle(menuName,menuControls))
  disp(sprintf('(mlrAdjustGUI:addMenu) Already added menu item: %s',menuName));
  return
end

% check to see if the location has a / on the end of it, which
% means to add it *underneath* the location specified
if menuLocation(end) == '/'
  % add the menu to the parent
  hAdded = uimenu(h,'Label',menuName);
  % and reorder to top
  hChildren = get(h,'Children');
  hChildren = [hChildren(2:end)' hAdded]';
  set(h,'Children',hChildren);
else
  
  % get the parent
  hParent = get(h,'Parent');

  % add the menu to the parent
  hAdded = uimenu(hParent,'Label',menuName);

  % reorder the children so that the item created is below the menuLocation
  hChildren = get(hParent,'Children');

  % get where these menu items live
  hAddedNum = find(hAdded==hChildren);
  hLocationNum = find(h==hChildren);

  % now create a reordered children list
  hNewChildren = [];
  for i = 1:length(hChildren)
    % if this one is the location menu item, then add the new child after it
    if i == hLocationNum
      hNewChildren(end+1) = hChildren(hAddedNum);
    end
    % add all the menu items except for the added one to the new children list
    if i ~= hAddedNum
      hNewChildren(end+1) = hChildren(i);
    end
  end
  
  % now set the children so that everything will be in the right order
  set(hParent,'Children',hNewChildren');
end

% add all the properties
for i = 1:2:length(menuProperties)
  set(hAdded,menuProperties{i},menuProperties{i+1});
end

% display what we have done
disp(sprintf('(mlrAdjustGUI:addMenu) Added menu: %s',menuName));

%%%%%%%%%%%%%%%%%%%%%
%%%   removeMenu  %%%
%%%%%%%%%%%%%%%%%%%%%
function removeMenu(args,menuControls)

% check length of arguments
if length(args) < 1
  disp(sprintf('(mlrAdjustGUI:removeMenu) Requires 1 argument: menuLocation'));
  return
else
  % name the arguments
  menuLocation = args{1};
end

% go look for the item location
h = getHandle(menuLocation,menuControls);
% if not found, then print warning, return
if isempty(h)
  disp(sprintf('(mlrAdjustGUI:removeMenu) Could not find menu location: %s',menuLocation));
  return
end

set(h,'visible','off')

% display what we have done
disp(sprintf('(mlrAdjustGUI:removeMenu) Removed menu: %s',menuLocation));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  setItemProperty   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function setItemProperty(args,uiControls,menuControls)

% check arguments
if length(args) ~= 3
  disp(sprintf('(mlrAdjustGUI:setItemProperty) Requires 3 arguments: controlName, propertyName, propertyValue'));
  return
else
  % name the arguments
  controlName = args{1};
  propertyName = args{2};
  propertyValue = args{3};
end

% go look for the control
h = getHandle(controlName,uiControls,menuControls);

% if not found, then print warning, return
if isempty(h),
  disp(sprintf('(mlrAdjustGUI:setItemProperty) Could not find control: %s',controlName));
  return
end

controlType = get(h,'type');
if ~isempty(propertyName)
    propertyName = lower(propertyName);
end    

% check if the property exists
fieldNames = lower(fieldnames(get(h)));

if ~ismember(propertyName,fieldNames)
  % if not, warn and continue
  if isstr(propertyName)
      disp(sprintf('(mlrAdjustGUI) *** Could not find property %s of %s %s ***',propertyName,controlType,controlName));
  else
    disp(sprintf('(mlrAdjustGUI) *** Could not find correct property of %s %s. Did you pass in a string indicating a property to set ***',controlType,controlName));
  end
  return
end

% if it does, make sure we have a valid property to set it to
%if ~any(strcmp(propertyName,{'Callback','TooltipString'}))
validValues = set(h,propertyName);
if ~isempty(validValues)
  tf = false;
  for j = 1:length(validValues) 
    if isequal(propertyValue,validValues{j}),tf = true;end
  end
  if ~tf
    disp(sprintf('(mlrAdjustGUI) Property value:'))
    disp(propertyValue);
    disp(sprintf('is invalid for %s %s property %s. Valid values are: ',controlType,controlName,propertyName));
    disp(validValues);
    return
  end
end

% if we got here, then the value is ok so set it
disp(sprintf('(mlrAdjustGUI) Setting property %s of %s %s',propertyName,controlType,controlName));
set(h,propertyName,propertyValue);

%%%%%%%%%%%%%%%%%%%%%%%
%%%  getUiControls  %%%
%%%%%%%%%%%%%%%%%%%%%%%
function uiControls = getUiControls(f,verbose)

if nargin == 1,verbose = 0;end
uiControls = [];

% get the gui data
h = guidata(f);
itemNames = fieldnames(h);

% for each guidata, check if it is a handle and
% not a menu item
for i = 1:length(itemNames)
  for j = 1:length(h.(itemNames{i}))
    if ishandle(h.(itemNames{i})(j))
      % check if it is not a menu item
      if isequal(get(h.(itemNames{i})(j),'Type'),'uicontrol')
	% keep track of this one.
	uiControls(end+1).tag = itemNames{i};
	uiControls(end).fieldNum = j;
	uiControls(end).h = h.(itemNames{i})(j);
      end
    end
  end
end

% print out info
if verbose
  for i = 1:length(uiControls)
    if uiControls(i).fieldNum == 1
      disp(sprintf('(mlrAdjustGUI) Found uicontrol: %i:%s',i,uiControls(i).tag));
      disp(get(uiControls(i).h,'Callback'));
    else
      disp(sprintf('(mlrAdjustGUI) Found uicontrol: %i:%s %i',i,uiControls(i).tag,uiControls(i).fieldNum));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  getMenuControls  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
function menuControls = getMenuControls(f,menuControls)

% function can be called recursively, if being called
% at top level, default lists.
if nargin == 1
  menuControls = {};
  topLabel = [];
else
  topLabel = menuControls(end).fullLabel;
end

% get list of all children
c = get(f,'Children');
% see if there are anu menus
for i = 1:length(c)
  % if it is a menu item
  if isequal(get(c(i),'Type'),'uimenu')
    % then get its label
    menuControls(end+1).fullLabel = sprintf('%s/%s',topLabel,get(c(i),'Label'));
    menuControls(end).label = get(c(i),'Label');
    menuControls(end).tag = get(c(i),'Tag');
    menuControls(end).h = c(i);
    % and go look for submenus
    menuControls = getMenuControls(c(i),menuControls);
  end
end

