% mlrAnatomyPlugin
%
%        $Id:$ 
%      usage: mlrAnatomyPlugin(action,<v>)
%         by: justin gardner & franco pestilli
%       date: 09/09/2014
%    purpose: Plugin function for LiFE
%
function retval = mlrAnatomyPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help mlrAnatomyPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrAnatomyPlugin) Need a valid view to install plugin'));
  else
    % create a panel - this will be used for adding some UI
    % controls to the right side of the figure for displaying
    % multiple surfaces at once.
    mlrAdjustGUI(v,'add','panel','Multiple base display',.5);
    % add the popup with the names of bases
    mlrAdjustGUI(v,'add','control','multiBaseListbox','panel','Multiple base display','style','popupmenu','position', [0.01    0.92    0.98   0.07 ],'Callback',@multiBaseListboxSelect);
    % add checkbox for multi base viewing
    mlrAdjustGUI(v,'add','control','multiBaseCheckbox','panel','Multiple base display','style','checkbox','position', [0.01    0.84    0.98   0.07 ],'String','MultiDisplay','Callback',@multiBaseCheckbox);
    % add slider for alpha
    mlrAdjustGUI(v,'add','control','multiBaseAlphaText','panel','Multiple base display','style','text','position', [0.01    0.76    0.2   0.07 ],'String','Alpha');
    mlrAdjustGUI(v,'add','control','multiBaseAlphaSlider','panel','Multiple base display','style','slider','position', [0.22    0.76    0.57   0.07 ],'String','Alpha','SliderStep',[0.1 0.25],'Callback',@multiBaseAlpha);
    mlrAdjustGUI(v,'add','control','multiBaseAlphaEdit','panel','Multiple base display','style','edit','position', [0.8    0.76    0.19   0.07 ],'Callback',@multiBaseAlpha);
    % add overlay popup (for setting the base pseudo color - or overlay)
    colors = color2RGB;
    colors = {'none' colors{:}};
    mlrAdjustGUI(v,'add','control','multiBaseOverlayText','panel','Multiple base display','style','text','position', [0.01    0.68    0.2   0.07 ],'String','Overlay');
    mlrAdjustGUI(v,'add','control','multiBaseOverlay','panel','Multiple base display','style','popup','position', [0.22    0.68    0.72   0.07 ],'String',colors,'Callback',@multiBaseOverlay);
    % add color alpha
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaText','panel','Multiple base display','style','text','position', [0.01    0.6    0.2   0.07 ],'String','Overlay Alpha');
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaSlider','panel','Multiple base display','style','slider','position', [0.22    0.6    0.57   0.07 ],'String','Alpha','SliderStep',[0.1 0.25],'Callback',@multiBaseOverlayAlpha);
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaEdit','panel','Multiple base display','style','edit','position', [0.8    0.6    0.19   0.07 ],'Callback',@multiBaseOverlayAlpha);
    

    % ROI controls
    %mlrAdjustGUI(v,'add','control','roiBaseListBox','panel','Multiple base display','style','listbox','position', [0.02    0.1    0.96   0.38 ],'Callback',@roiListboxSelect,'Max',2);

    % add the callback that will tell the above listbox when new
    % bases have been added
    v = viewSet(v,'callback','baseChange',@mlrAnatomyBaseChange);
    % also register a change when someone switches the curBase
    v = viewSet(v,'callback','curBaseChange',@mlrAnatomyBaseChange);

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(mlrAnatomyPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%
%    baseChange    %
%%%%%%%%%%%%%%%%%%%%
function v = mlrAnatomyBaseChange(v)

% get control
baseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');

% get current base names
baseNames = viewGet(v,'baseNames');

% no bases
if isempty(baseNames)
  % set its values
  set(baseListbox,'String','');

  % and set which ones are selected
  set(baseListbox,'Value',[]);
  return
end
  
% get the curBase and what type it is
curBase = viewGet(v,'curBase');
curBaseType = viewGet(v,'baseType');
curBaseName = viewGet(v,'baseName');

% get the current selected base
baseListboxNames = get(baseListbox,'String');
baseListboxValue = get(baseListbox,'Value');
if ~isempty(baseListboxValue) && (baseListboxValue>=1) && (baseListboxValue<=length(baseListboxNames))
  selectedBaseName = baseListboxNames{baseListboxValue};
else
  selectedBaseName = '';
end

% get the types for everyone
baseType = [];
for iBase = 1:viewGet(v,'numBase')
  baseType(iBase) = viewGet(v,'baseType',iBase);
end
    
% get surfaces
baseNames = {baseNames{(baseType==2)|(baseType==3)}};
baseType = baseType((baseType==2)|(baseType==3));

% set values
set(baseListbox,'String',baseNames);

% and set which ones are selected
selectedBaseNum = find(strcmp(selectedBaseName,baseNames));
if ~isempty(selectedBaseNum)
  set(baseListbox,'Value',selectedBaseNum);
else
  set(baseListbox,'Value',1);
end

% call multiBaseListboxSelect to set the base info
multiBaseListboxSelect(baseListbox,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiBaseListboxSelect     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function multiBaseListboxSelect(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get control
multiBaseListbox = hObject;

% get select value
baseNames = get(multiBaseListbox,'String');
selectedVal = get(multiBaseListbox,'Value');

% validate selection val
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames))
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % get base properties
  base = viewGet(v,'base',baseNum);

  % set the various properties
  mlrAdjustGUI(v,'set','multiBaseCheckbox','Value',base.multiDisplay);
  mlrAdjustGUI(v,'set','multiBaseAlphaSlider','Value',base.alpha);
  mlrAdjustGUI(v,'set','multiBaseAlphaEdit','String',base.alpha);
  mlrAdjustGUI(v,'set','multiBaseOverlayAlphaSlider','Value',base.overlayAlpha);
  mlrAdjustGUI(v,'set','multiBaseOverlayAlphaEdit','String',base.overlayAlpha);
  if ~isempty(base.overlay) && isstr(base.overlay)
    % for overlays that are strings, this is usually a color name
    % so first get the color names in the multiBaseOverlay control
    % (These originate from color2RGB
    multiBaseOverlay = mlrAdjustGUI(v,'get','multiBaseOverlay');
    baseColorNames = get(multiBaseOverlay,'String');
    % find a match
    if isempty(base.overlay)
      baseMatch = 1;
    else
      baseMatch = find(strcmp(base.overlay,baseColorNames));
    end
      
    if ~isempty(baseMatch)
      % set the value
      set(multiBaseOverlay,'Value',baseMatch);
    else
      % add the name if it does not already exist
      baseColorNames{end+1} = base.overlay;
      set(multiBaseOverlay,'String',baseColorNames);
      set(multiBaseOverlay,'Value',length(baseColorNames));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiBaseCheckbox    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function multiBaseCheckbox(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get the multiBaselistbox and figure out what base is selected
multiBaseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseNames = get(multiBaseListbox,'String');
selectedVal = get(multiBaseListbox,'Value');

% validate selection val
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames))
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % set the base multiBase
  v = viewSet(v,'baseMultiDisplay',get(hObject,'Value'),baseNum);

  % redisplay
  refreshMLRDisplay(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    multiBaseAlpha    %
%%%%%%%%%%%%%%%%%%%%%%%%
function multiBaseAlpha(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get the multiBaselistbox and figure out what base is selected
multiBaseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseNames = get(multiBaseListbox,'String');
selectedVal = get(multiBaseListbox,'Value');

% validate selection val
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames))
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % now see whether we are being called from slider or edit box
  if isequal(get(hObject,'Style'),'edit')
    % get baseAlpha from edit
    baseAlpha = str2num(get(hObject,'String'));
  else
    % get baseAlpha from slider
    baseAlpha = get(hObject,'Value');
    % round to nearest 1/100
    baseAlpha = round(baseAlpha*100)/100;
  end
  if ~isempty(baseAlpha) && (baseAlpha>=0) && (baseAlpha<=1)
    % set base alpha
    v = viewSet(v,'baseAlpha',baseAlpha,baseNum);
    % update the controls (depends on who changed)
    if isequal(get(hObject,'Style'),'edit')
      % set slider
      mlrAdjustGUI(v,'set','multiBaseAlphaSlider','Value',baseAlpha);
    else
      % or edit
      mlrAdjustGUI(v,'set','multiBaseAlphaEdit','String',baseAlpha);
    end
      
    % redisplay
    refreshMLRDisplay(v);
  else
    % bad value, reset
    baseAlpha = viewGet(v,'baseAlpha',baseNum);
    if isequal(get(hObject,'Style'),'edit')
      set(hObject,'String',baseAlpha);
    else
      set(hObject,'Value',baseAlpha);
    end
      
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiBaseOverlay   %
%%%%%%%%%%%%%%%%%%%%%%%%%
function multiBaseOverlay(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get the multiBaselistbox and figure out what base is selected
multiBaseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseNames = get(multiBaseListbox,'String');
selectedVal = get(multiBaseListbox,'Value');

% validate selection val
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames))
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % set the base color
  baseColors = get(hObject,'String');
  baseColor = baseColors{get(hObject,'Value')};
  if strcmp(baseColor,'none') baseColor = [];end
  
  % and set the base with that color
  v = viewSet(v,'baseOverlay',baseColor,baseNum);

  % redisplay
  refreshMLRDisplay(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiBaseOverlayAlpha   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function multiBaseOverlayAlpha(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get the multiBaselistbox and figure out what base is selected
multiBaseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseNames = get(multiBaseListbox,'String');
selectedVal = get(multiBaseListbox,'Value');

% validate selection val
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames))
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % now see whether we are being called from slider or edit box
  if isequal(get(hObject,'Style'),'edit')
    % get baseOverlayAlpha from edit
    baseOverlayAlpha = str2num(get(hObject,'String'));
  else
    % get baseOverlayAlpha from slider
    baseOverlayAlpha = get(hObject,'Value');
    % round to nearest 1/100
    baseOverlayAlpha = round(baseOverlayAlpha*100)/100;
  end
  
  % set value
  v = viewSet(v,'baseOverlayAlpha',baseOverlayAlpha);
  
  % validate value
  if ~isempty(baseOverlayAlpha) && (baseOverlayAlpha>=0) && (baseOverlayAlpha<=1)
    % set baseOverlayAlpha
    v = viewSet(v,'baseOverlayAlpha',baseOverlayAlpha,baseNum);
    % update the controls (depends on who changed)
    if isequal(get(hObject,'Style'),'edit')
      % set slider
      mlrAdjustGUI(v,'set','multiBaseOverlayAlphaSlider','Value',baseOverlayAlpha);
    else
      % or edit
      mlrAdjustGUI(v,'set','multiBaseOverlayAlphaEdit','String',baseOverlayAlpha);
    end
      
    % redisplay
    refreshMLRDisplay(v);
  else
    % bad value, reset
    baseOverlayAlpha = viewGet(v,'baseOverlayAlpha',baseNum);
    if isequal(get(hObject,'Style'),'edit')
      set(hObject,'String',baseOverlayAlpha);
    else
      set(hObject,'Value',baseOverlayAlpha);
    end
      
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    roiListboxSelect     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiListboxSelect(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get control
roiListbox = hObject;

% get current selection
selectedBases = get(roiListbox,'Value');
