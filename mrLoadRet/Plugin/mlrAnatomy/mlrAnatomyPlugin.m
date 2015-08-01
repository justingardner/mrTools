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

    % add controls for fascicles
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleIntersect','panel','Multiple base display','style','pushbutton','position', [0.01    0.52    0.98   0.07 ],'String','Calculate intersect','Callback',@mlrAnatomyFascicleIntersect);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleN','panel','Multiple base display','style','text','position', [0.01    0.44    0.98   0.07 ],'HorizontalAlignment','Center');
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleDisplay','panel','Multiple base display','style','popupmenu','position', [0.01    0.36    0.98   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleDisplay);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleRestrict','panel','Multiple base display','style','popupmenu','position', [0.01    0.28    0.98   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleRestrict);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMinText','panel','Multiple base display','style','text','position', [0.01    0.20    0.1   0.07 ],'HorizontalAlignment','Center','String','Min');
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMinSlider','panel','Multiple base display','style','slider','position', [0.11    0.20    0.25   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleMinmax);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMinEdit','panel','Multiple base display','style','edit','position', [0.37    0.20    0.12   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleMinmax);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMaxText','panel','Multiple base display','style','text','position', [0.5    0.20    0.1   0.07 ],'HorizontalAlignment','Center','String','Max');
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMaxSlider','panel','Multiple base display','style','slider','position', [0.61    0.20    0.25   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleMinmax);
    mlrAdjustGUI(v,'add','control','mlrAnatomyFascicleMaxEdit','panel','Multiple base display','style','edit','position', [0.87    0.20    0.12   0.07 ],'HorizontalAlignment','Center','Callback',@mlrAnatomyFascicleMinmax);

    
    % add plane rotation
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundText','panel','Multiple base display','style','text','position', [0.01    0.48    0.2   0.07 ],'String','Rotate around');
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundX','panel','Multiple base display','style','radio','position', [0.22    0.52    0.15   0.07 ],'String','X','Callback',@mlrAnatomyRotateAround);
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundY','panel','Multiple base display','style','radio','position', [0.42    0.52    0.15   0.07 ],'String','Y','Callback',@mlrAnatomyRotateAround);
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundZ','panel','Multiple base display','style','radio','position', [0.62    0.52    0.15   0.07 ],'String','Z','Callback',@mlrAnatomyRotateAround);
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundSlider','panel','Multiple base display','style','slider','position', [0.22    0.44    0.57   0.07 ],'String','Rotate','SliderStep',[1 15]/360,'Callback',@mlrAnatomyRotateAround,'Max',360,'Min',0);
    mlrAdjustGUI(v,'add','control','mlrAnatomyRotateAroundEdit','panel','Multiple base display','style','edit','position', [0.8    0.44    0.19   0.07 ],'Callback',@mlrAnatomyRotateAround);

    % add plane center position
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterText','panel','Multiple base display','style','text','position', [0.01    0.3    0.2   0.07 ],'String','Center');
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterX','panel','Multiple base display','style','radio','position', [0.22    0.36    0.15   0.07 ],'String','X','Callback',@mlrAnatomyCenter);
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterY','panel','Multiple base display','style','radio','position', [0.42    0.36    0.15   0.07 ],'String','Y','Callback',@mlrAnatomyCenter);
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterZ','panel','Multiple base display','style','radio','position', [0.62    0.36    0.15   0.07 ],'String','Z','Callback',@mlrAnatomyCenter);
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterSlider','panel','Multiple base display','style','slider','position', [0.22    0.28    0.57   0.07 ],'String','Rotate','SliderStep',[1 10]/200,'Callback',@mlrAnatomyCenter,'Max',100,'Min',-100);
    mlrAdjustGUI(v,'add','control','mlrAnatomyCenterEdit','panel','Multiple base display','style','edit','position', [0.8    0.28    0.19   0.07 ],'Callback',@mlrAnatomyCenter);

    % ROI controls
    %mlrAdjustGUI(v,'add','control','roiBaseListBox','panel','Multiple base display','style','listbox','position', [0.02    0.1    0.96   0.38 ],'Callback',@roiListboxSelect,'Max',2);

    % add a menu item to import rois from freesurfer
    mlrAdjustGUI(v,'add','menu','Import Freesurfer Label','/File/ROI/Import','Callback',@mlrAnatomyImportFreesurferLabel);

    % add a menu item to make a planer base anatomy
    mlrAdjustGUI(v,'add','menu','Make plane','/File/Base anatomy/Import surface','Callback',@mlrAnatomyMakePlaneBase,'Separator','on');
    
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
  set(baseListbox,'String',{'empty'});

  % and set which ones are selected
  set(baseListbox,'Value',1);
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
if isempty(baseNames)
  % set its values
  set(baseListbox,'String',{'empty'});

  % and set which ones are selected
  set(baseListbox,'Value',1);
  return
else
  
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
end

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
if ~isempty(selectedVal) && (selectedVal>0) && (selectedVal <= length(baseNames)) && ~isequal(baseNames{selectedVal},'empty')
  % get base num
  baseNum = viewGet(v,'baseNum',baseNames{selectedVal});

  % get base properties
  base = viewGet(v,'base',baseNum);

  % fascicle controls
  if isempty(base.fascicles)
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleIntersect','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleN','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleDisplay','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleRestrict','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinText','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinEdit','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxText','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxEdit','Visible','off');
  else
    % set display range if not set
    if ~isfield(base.fascicles,'displayMin');
      base.fascicles.displayMin = 1;
      base.fascicles.displayMax = base.fascicles.n;
    end
    % check to see if we are doing a subset
    restrictValue = 1;
    if (base.fascicles.displayMin > 1) || (base.fascicles.displayMax > base.fascicles.n)
      % see if we are doing one fascicle
      if base.fascicles.displayMin == base.fascicles.displayMax
	restrictValue = 2+base.fascicles.displayMin;
      else
	restrictValue = 2;
      end
    end
    % set gui
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleIntersect','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleN','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleN','String',sprintf('N=%i',base.fascicles.n));
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleDisplay','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleDisplay','String',{'Fascicle'});
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleDisplay','Value',1);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleRestrict','Visible','on');
    % make a restriction list with all possible fascicles
    restrictList = {'All','Subset'};
    for iFascicle = 1:base.fascicles.n
      restrictList{end+1}= sprintf('fascicle %04i',iFascicle');
    end
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleRestrict','String',restrictList);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleRestrict','Value',restrictValue);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinText','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Min',1);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Max',base.fascicles.n);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Value',base.fascicles.displayMin);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinEdit','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinEdit','String',sprintf('%i',base.fascicles.displayMin));
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxText','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Min',1);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Max',base.fascicles.n);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Value',base.fascicles.displayMax);
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxEdit','String',sprintf('%i',base.fascicles.displayMax));
    mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxEdit','Visible','on');
  end
  
  %  if the base is not a plane, then turn off all the plane controls
  if isempty(base.plane)
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundText','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundX','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundY','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundZ','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundSlider','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','Visible','off');
  
    mlrAdjustGUI(v,'set','mlrAnatomyCenterText','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterX','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterY','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterZ','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterSlider','Visible','off');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','Visible','off');
else
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundText','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundX','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundY','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundZ','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundSlider','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','Visible','on');
    % set the radio to default to x
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundX','Value',1);
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundY','Value',0);
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundZ','Value',0);
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundSlider','Value',base.plane.xRot);
    mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','String',sprintf('%0.1f',base.plane.xRot));
  
    mlrAdjustGUI(v,'set','mlrAnatomyCenterText','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterX','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterY','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterZ','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterSlider','Visible','on');
    mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','Visible','on');
    % set the radio to default to x
    mlrAdjustGUI(v,'set','mlrAnatomyCenterX','Value',1);
    mlrAdjustGUI(v,'set','mlrAnatomyCenterY','Value',0);
    mlrAdjustGUI(v,'set','mlrAnatomyCenterZ','Value',0);
    mlrAdjustGUI(v,'set','mlrAnatomyCenterSlider','Value',base.plane.x);
    mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','String',sprintf('%0.1f',base.plane.x));
end
  
  % set the various properties
  mlrAdjustGUI(v,'set','multiBaseCheckbox','Value',base.multiDisplay);
  % set alpha
  if isempty(base.alpha) alpha = 1;else alpha = base.alpha;end
  mlrAdjustGUI(v,'set','multiBaseAlphaSlider','Value',alpha);
  mlrAdjustGUI(v,'set','multiBaseAlphaEdit','String',alpha);
  % set overlayAlpha
  if isempty(base.overlayAlpha) overlayAlpha = 1;else overlayAlpha = base.overlayAlpha;end
  mlrAdjustGUI(v,'set','multiBaseOverlayAlphaSlider','Value',overlayAlpha);
  mlrAdjustGUI(v,'set','multiBaseOverlayAlphaEdit','String',overlayAlpha);
  % get overlay
  if ~isempty(base.overlay) && isstr(base.overlay)
    overlay = base.overlay;
  else
    overlay = 'none';
  end
  % for overlays that are strings, this is usually a color name
  % so first get the color names in the multiBaseOverlay control
  % (These originate from color2RGB
  multiBaseOverlay = mlrAdjustGUI(v,'get','multiBaseOverlay');
  baseColorNames = get(multiBaseOverlay,'String');
  % find a match
  if isempty(overlay)
    baseMatch = 1;
  else
    baseMatch = find(strcmp(overlay,baseColorNames));
  end
      
  if ~isempty(baseMatch)
    % set the value
    set(multiBaseOverlay,'Value',baseMatch);
  else
    % add the name if it does not already exist
    baseColorNames{end+1} = overlay;
    set(multiBaseOverlay,'String',baseColorNames);
    set(multiBaseOverlay,'Value',length(baseColorNames));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrAnatomyImportFreesurferLabel   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyImportFreesurferLabel(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

disp(sprintf('(mlrAnatomyImportFreesurferLabel) Not yet implemented'));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrAnatomMakePlaneBase   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyMakePlaneBase(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% find possible canonicals
baseNames = {};
baseNumVoxels = [];
for iBase = 1:viewGet(v,'numBase')
  if viewGet(v,'baseType',iBase) == 0
    baseNames{end+1} = viewGet(v,'baseName',iBase);
    baseNumVoxels(end+1) = prod(viewGet(v,'baseDims',iBase));
  end
end
if isempty(baseNames)
  mrWarnDlg('(mlrAnatomyMakePlaneBase) Could not find any volume anatomies to make a plane base from.');
  return
end

% put the base with the largest number of voxels on top of list (this should
% be the canonical)
[maxVoxels maxIndex] = max(baseNumVoxels);
baseNames = putOnTopOfList(baseNames{maxIndex},baseNames);

if length(baseNames) > 1
  % now put up a dialog box for subject to select parameters of plane
  paramsInfo = {};
  paramsInfo{end+1} = {'baseName',baseNames,'Select the name of the base that you want to use to create the plane base from. This will reslice that base to make the desired plane'};

  % choose parameters
  params = mrParamsDialog(paramsInfo,'Choose base to make plane from');
  if isempty(params),return,end
else
  params.baseName = baseNames{1};
end

% get the base that we are going to make this plane from
canonical = viewGet(v,'base',viewGet(v,'baseNum',params.baseName));

% copy info into to create this new base
b = canonical;
b.plane.canonical = canonical.data;
canonicalDims = size(b.plane.canonical);
b.coordMap.dims = canonicalDims;
b.coordMap.innerSurfaceFileName = '';
b.coordMap.innerCoordsFileName = '';
b.coordMap.outerSurfaceFileName = '';
b.coordMap.outerCoordsFileName = '';
b.coordMap.curvFileName = '';
b.coordMap.anatFileName = canonical.name;
b.coordMap.path = '';
b.type = 2;

% set default position and rotation
b.plane.x = 0;b.plane.y = 0;b.plane.z = 0;
b.plane.zRot = 0;b.plane.yRot = 0;b.plane.xRot = 0;
b.plane.width = canonicalDims(2);
b.plane.height = canonicalDims(1);

% calculate the coordinates and get resliced data
b = calcPlaneCoords(b);

% set name
b.name = sprintf('%sPlane',canonical.name);

% and install the base
v = viewSet(v,'newBase',b);
refreshMLRDisplay(v);

%%%%%%%%%%%%%%%%%%%%%%%%%
%    calcPlaneCoords    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function b = calcPlaneCoords(b)

% get dimensions of canonical
canonicalDims = size(b.plane.canonical);

% make rotation matrix
c = cos(pi*b.plane.zRot/180);
s = sin(pi*b.plane.zRot/180);
rotxy = [c -s 0 0;s  c 0 0;0  0 1 0;0  0 0 1];
c = cos(pi*b.plane.xRot/180);
s = sin(pi*b.plane.xRot/180);
rotyz = [1  0  0 0;0  c -s 0;0  s  c 0;0  0  0 1];
c = cos(pi*b.plane.yRot/180);
s = sin(pi*b.plane.yRot/180);
rotxz = [c  0 -s 0;0  1  0 0;s  0  c 0;0  0  0  1];
r = rotxy*rotyz*rotxz;

% now we make the surface positions
x = 1:b.plane.width;
y = 1:b.plane.height;
[x y] = meshgrid(x,y);
z = ones(size(x))*round(canonicalDims(3)/2);
x = x';y = y';z = z';

% now get coordinates (xy swaping because of image coordinates)
coords = [y(:) x(:) z(:)];
% make homogenous and rotate according to rotation matrix
coords(:,1) = coords(:,1) - canonicalDims(1)/2;
coords(:,2) = coords(:,2) - canonicalDims(2)/2;
coords(:,3) = coords(:,3) - canonicalDims(3)/2;
coords(:,4) = 1;
coords = r * coords';
coords = coords(1:3,:)';
coords(:,1) = coords(:,1) + canonicalDims(1)/2;
coords(:,2) = coords(:,2) + canonicalDims(2)/2;
coords(:,3) = coords(:,3) + canonicalDims(3)/2;

% also move center to user specified center
coords(:,1) = coords(:,1) + b.plane.x;
coords(:,2) = coords(:,2) + b.plane.y;
coords(:,3) = coords(:,3) + b.plane.z;

% and put them into the coordMap strucutre
b.coordMap.coords(1,:,1,1) = coords(:,1);
b.coordMap.coords(1,:,1,2) = coords(:,2);
b.coordMap.coords(1,:,1,3) = coords(:,3);
b.coordMap.outerCoords = b.coordMap.coords;
b.coordMap.innerCoords = b.coordMap.coords;
b.coordMap.innerVtcs = coords;
b.coordMap.outerVtcs = coords;

% now connect all these vertices with appropriate triangles
width = canonicalDims(1);
height = canonicalDims(2);

% from all the triangles that are the upper/left corners
triCoords1a = repmat(1:height-1,width-1,1)+repmat((0:width-2)*height,height-1,1)';
triCoords2a = repmat(2:height,width-1,1)+repmat((0:width-2)*height,height-1,1)';
triCoords3a = repmat(1:height-1,width-1,1)+repmat((1:width-1)*height,height-1,1)';

% from all the triangles that are the bottom/right corners
triCoords1b = repmat(2:height,width-1,1)+repmat((0:width-2)*height,height-1,1)';
triCoords2b = repmat(1:height-1,width-1,1)+repmat((1:width-1)*height,height-1,1)';
triCoords3b = repmat(2:height,width-1,1)+repmat((1:width-1)*height,height-1,1)';

% now put those into the coordmap
b.coordMap.tris = [];
b.coordMap.tris(:,1) = [triCoords1a(:);triCoords1b(:)];
b.coordMap.tris(:,2) = [triCoords2a(:);triCoords2b(:)];
b.coordMap.tris(:,3) = [triCoords3a(:);triCoords3b(:)];

% now get the data for the points
b.data = interp3(b.plane.canonical,coords(:,2),coords(:,1),coords(:,3),'linear',0);
b.data = b.data(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyRotateAround    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyRotateAround(hObject,eventdata)

recompute = false;

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get radio button control values
hX = mlrAdjustGUI(v,'get','mlrAnatomyRotateAroundX');
hY = mlrAdjustGUI(v,'get','mlrAnatomyRotateAroundY');
hZ = mlrAdjustGUI(v,'get','mlrAnatomyRotateAroundZ');

% get selected base  
baseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseList = get(baseListbox,'String');
baseValue = get(baseListbox,'Value');
if baseValue < 1,return,end
baseNum = viewGet(v,'baseNum',baseList{baseValue});
b = viewGet(v,'base',baseNum);

% see if this was a change in the radio button control
if any(hObject == [hX hY hZ])
  set(hX,'Value',0);
  set(hY,'Value',0);
  set(hZ,'Value',0);
  set(hObject,'Value',1);
  
  if hObject==hX
    val = b.plane.xRot;
  elseif hObject==hY
    val = b.plane.yRot;
  else
    val = b.plane.zRot;
  end
  
  % set the edit and the slider
  mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','String',sprintf('%i',val));
% slider control
elseif hObject==mlrAdjustGUI(v,'get','mlrAnatomyRotateAroundSlider')
  % get value
  val = round(get(hObject,'Value'));
  % set edit
  mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','String',sprintf('%0.1f',val));
  % set rotation in base
  if get(hX,'Value')
    b.plane.xRot = val;
  elseif get(hY,'Value')
    b.plane.yRot = val;
  else
    b.plane.zRot = val;
  end
  recompute = true;
% must be edit control then
else
  % get value
  val = round(mrStr2num(get(hObject,'String'))*10)/10;
  val = min(max(0,val),360);
  % set edit and slider
  mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyRotateAroundEdit','String',sprintf('%0.1f',val));
  % set rotation in base
  if get(hX,'Value')
    b.plane.xRot = val;
  elseif get(hY,'Value')
    b.plane.yRot = val;
  else
    b.plane.zRot = val;
  end
  recompute = true;
end

% recompute and redisplay
if recompute  
  % recompute
  b = calcPlaneCoords(b);
  % reset the base and cache
  v = viewSet(v,'base',b,baseNum);
  v = viewSet(v,'baseCache','clear',b.name);
  v = viewSet(v,'overlayCache','clear',b.name);
  v = viewSet(v,'roiCache','clear',b.name);
  % and redisplay
  refreshMLRDisplay(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyCenter    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyCenter(hObject,eventdata)

recompute = false;

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get radio button control values
hX = mlrAdjustGUI(v,'get','mlrAnatomyCenterX');
hY = mlrAdjustGUI(v,'get','mlrAnatomyCenterY');
hZ = mlrAdjustGUI(v,'get','mlrAnatomyCenterZ');

% get selected base  
baseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseList = get(baseListbox,'String');
baseValue = get(baseListbox,'Value');
if baseValue < 1,return,end
baseNum = viewGet(v,'baseNum',baseList{baseValue});
b = viewGet(v,'base',baseNum);

% see if this was a change in the radio button control
if any(hObject == [hX hY hZ])
  set(hX,'Value',0);
  set(hY,'Value',0);
  set(hZ,'Value',0);
  set(hObject,'Value',1);
  
  if hObject==hX
    val = b.plane.x;
  elseif hObject==hY
    val = b.plane.y;
  else
    val = b.plane.z;
  end
  
  % set the edit and the slider
  mlrAdjustGUI(v,'set','mlrAnatomyCenterSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','String',sprintf('%0.1f',val));
% slider control
elseif hObject==mlrAdjustGUI(v,'get','mlrAnatomyCenterSlider')
  % get value
  val = round(get(hObject,'Value'));
  % set edit
  mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','String',sprintf('%0.1f',val));
  % set rotation in base
  if get(hX,'Value')
    b.plane.x = val;
  elseif get(hY,'Value')
    b.plane.y = val;
  else
    b.plane.z = val;
  end
  recompute = true;
% if this is the edit control
else
  % get value
  val = round(mrStr2num(get(hObject,'String'))*10)/10;
  val = min(max(-100,val),100);
  % set edit and slider
  mlrAdjustGUI(v,'set','mlrAnatomyCenterSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyCenterEdit','String',sprintf('%0.1f',val));
  % set rotation in base
  if get(hX,'Value')
    b.plane.x = val;
  elseif get(hY,'Value')
    b.plane.y = val;
  else
    b.plane.z = val;
  end
  recompute = true;
end

% recompute and redisplay
if recompute  
  % recompute
  b = calcPlaneCoords(b);
  % reset the base and cache
  v = viewSet(v,'base',b,baseNum);
  v = viewSet(v,'baseCache','clear',b.name);
  v = viewSet(v,'overlayCache','clear',b.name);
  v = viewSet(v,'roiCache','clear',b.name);
  % and redisplay
  refreshMLRDisplay(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyFascicleIntersect  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyFascicleIntersect(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get all bases in the list box
baseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');
baseSurfaces = get(baseListbox,'String');

% get which one is selected and put it on top of list
selectedVal = get(baseListbox,'Value');
if ~isempty(selectedVal) && (selectedVal>=1) && (selectedVal<=length(baseSurfaces))
  baseSurfaces = putOnTopOfList(baseSurfaces{selectedVal},baseSurfaces);
end

% just get ones that have fascicles
fascicles = {};
for iBase = 1:length(baseSurfaces)
  b = viewGet(v,'base',viewGet(v,'baseNum',baseSurfaces{iBase}));
  if ~isempty(b.fascicles)
    fascicles{end+1} = baseSurfaces{iBase};
  end
end
if isempty(fascicles)
  mrWarnDlg('(mlrAnatomyPluginFascicleIntersect) Could not find any fasicile anatomies');
  return
end

% make combine dialog
paramsInfo = {};
paramsInfo{end+1} = {'fascicle',fascicles,'Fascicle base to compute intersection on. That is, this is the one whose fascicles will change dependent on what action we take below'};
paramsInfo{end+1} = {'intersectSurface',baseSurfaces,'Fascicle base to check intersection with.'};
paramsInfo{end+1} = {'intersectionType',{'intersect','does not intersect'},'Choose what type of intersection to check'};
paramsInfo{end+1} = {'doIntersection',0,'type=pushbutton','buttonString=Do intersection','Press to do the intersection','callback',@mlrAnatomyPluginDoFascicleIntersect,'callbackArg',v,'passParams=1'};

% put up dialog
params = mrParamsDialog(paramsInfo,'Fascicle intersection tool');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyPluginDoFasicleIntersect    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = mlrAnatomyPluginDoFascicleIntersect(v,params)

% get the base we are going to modify
baseNum = viewGet(v,'baseNum',params.fascicle);
b = viewGet(v,'base',baseNum);
f = b.fascicles;

% get intersection surface
intersectBaseNum = viewGet(v,'baseNum',params.intersectSurface);
intersect = viewGet(v,'baseSurface',intersectBaseNum);

% get xform from intersection base to fascicles
base2base = viewGet(v,'base2base',baseNum,intersectBaseNum);
swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];

d = -inf(1,f.n);

% get vertices of intersection surface
v2 = intersect.vtcs;
v2(:,4) = 1;
v2 = swapXY*base2base*swapXY*v2';
v2 = v2(1:3,:)';
v2n = size(v2,1);

%surf1 = f.patches{1};
%surf2.vertices = v2;
%surf2.faces = intersect.tris;
%mlrSmartfig('mlrAnatomyPluginHuh','reuse');clf;
%patch('vertices',surf1.vertices,'faces',surf1.faces);
%hold on
%patch('vertices',surf2.vertices,'faces',surf2.faces);

mlrSmartfig('mlrAnatomyPlugin','reuse');clf;
% draw the intersection surface
patch('Vertices',v2,'Faces',intersect.tris);
hold on
h = [];
for iFascicle = 1:f.n
  % get vertices of this fascicle
  v1 = f.patches{iFascicle}.vertices;
  v1n = size(v1,1);
  % now make a v1nxv2n matrix putting each x,y and z of v2 into rows
  v2x = repmat(v2(:,1),1,v1n);
  v2y = repmat(v2(:,2),1,v1n);
  v2z = repmat(v2(:,3),1,v1n);
  % same for v1, except put x,y and z into columns
  v1x = repmat(v1(:,2)',v2n,1);
  v1y = repmat(v1(:,1)',v2n,1);
  v1z = repmat(v1(:,3)',v2n,1);
  % take difference, then square and sum to get distance
  dist = sqrt((v2x-v1x).^2 + (v2y-v1y).^2 +(v2z-v1z).^2);
  d(iFascicle) = min(dist(:));
  disp(sprintf('%i: min = %0.2f',iFascicle,d(iFascicle)));
  if d(iFascicle)<1,c = 'g';else c='r';end
  if ~isempty(h),set(h,'edgealpha',0);set(h,'facealpha',0.1);end
 h = patch('Vertices',v1,'Faces',f.patches{iFascicle}.faces,'facecolor',c,'edgecolor',c,'facealpha',1,'edgealpha',1);
  if iFascicle==1
    keyboard
  else
    drawnow;
  end
end

% now put all fascicles vertices and triangles into one coordMap
nRunningTotalVertices = 0;
nRunningTotalTris = 0;

% store the intersection
if ~isfield(b.fascicles,'intersect')
  b.fascicles.intersect = [];
  iIntersect = 1;
else
  % find if there is an intersect with matching surface name
  % so we can go replace that one.
  iIntersect = find(strcmp(params.intersectSurface,{b.fascicles.intersect(:).intersectWith}));
  if isempty(iIntersect),iIntersect = length(b.fascicles.intersect)+1;end
end
% and add it
b.fascicles.intersect(iIntersect).intersectWith = params.intersectSurface;
b.fascicles.intersect(iIntersect).d = d;
v = viewSet(v,'base',b,baseNum);

disppercent(-inf,sprintf('(mlrAnatomyPlugin) Converting %i fascicles',f.n));
for iFascicle = 1:f.n
  if d(iFascicle) > 1
    % number of vertices and triangles
    nVertices = size(f.patches{iFascicle}.vertices,1);
    nTris = size(f.patches{iFascicle}.faces,1);
    % the data which is the grayscale value to color the fascicles with (rand for now)
    b.data = [b.data rand(1,nVertices)];
    % convert vertices to a coord map which has one x,y,z element for each possible
    % location on the surface (which actually is just a 1xnVerticesx1 image)
    % add these vertices to existing vertices
    b.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,1) = f.patches{iFascicle}.vertices(:,1);
    b.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,2) = f.patches{iFascicle}.vertices(:,2);
    b.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,3) = f.patches{iFascicle}.vertices(:,3);
    % these are the display vertices which are the same as the coords
    b.coordMap.innerVtcs(nRunningTotalVertices+1:nRunningTotalVertices+nVertices,:) = f.patches{iFascicle}.vertices;
    % triangle faces
    b.coordMap.tris(nRunningTotalTris+1:nRunningTotalTris+nTris,:) = (f.patches{iFascicle}.faces + nRunningTotalVertices);
    % update runing totals
    nRunningTotalVertices = nRunningTotalVertices + nVertices;
    nRunningTotalTris= nRunningTotalTris + nTris;
  end
  disppercent(iFascicle/f.n);
end
disppercent(inf);

% make right length
b.coordMap.tris = b.coordMap.tris(1:nRunningTotalTris,:);
b.coordMap.innerCoords = b.coordMap.innerCoords(:,1:nRunningTotalVertices,:,:,:);
b.coordMap.innerVtcs = b.coordMap.innerVtcs(1:nRunningTotalVertices,:);
b.data = b.data(:,1:nRunningTotalVertices);

% copy the inner to outer since they are all the same for fascicles
b.coordMap.outerCoords = b.coordMap.innerCoords;
b.coordMap.outerVtcs = b.coordMap.innerVtcs;

% make sure it is still a base
[tf b] = isbase(b);

% reset the base
%v = viewSet(v,'base',b,baseNum);
v = viewSet(v,'base',b,baseNum);
v = viewSet(v,'baseCache','clear',b.name);
v = viewSet(v,'overlayCache','clear',b.name);
v = viewSet(v,'roiCache','clear',b.name);
refreshMLRDisplay(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyFascicleDisplay    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyFascicleDisplay(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyFascicleRestrict    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyFascicleRestrict(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get control
baseListbox = mlrAdjustGUI(v,'get','multiBaseListbox');

% get the current selected base
baseListboxNames = get(baseListbox,'String');
baseListboxValue = get(baseListbox,'Value');
if ~isempty(baseListboxValue) && (baseListboxValue>=1) && (baseListboxValue<=length(baseListboxNames))
  selectedBaseName = baseListboxNames{baseListboxValue};
else
  selectedBaseName = '';
end
baseNum = viewGet(v,'baseNum',selectedBaseName);
if isempty(baseNum),'return',end

% get the selected value
val = get(hObject,'Value');

% then particular fascicle has been selected
if val > 2
  val = val-2;
  mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyFascicleMinEdit','String',sprintf('%i',val));
  mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxSlider','Value',val);
  mlrAdjustGUI(v,'set','mlrAnatomyFascicleMaxEdit','String',sprintf('%i',val));

  % get the selected base
  b = viewGet(v,'base',baseNum);
  f = b.fascicles;

  iFascicle = val;
  % number of vertices and triangles
  nVertices = size(f.patches{iFascicle}.vertices,1);
  nTris = size(f.patches{iFascicle}.faces,1);
  b.coordMap.innerCoords = [];
  b.coordMap.innerVtcs = [];
  % the data which is the grayscale value to color the fascicles with (rand for now)
  b.data = rand(1,nVertices);
  % convert vertices to a coord map which has one x,y,z element for each possible
  % location on the surface (which actually is just a 1xnVerticesx1 image)
  % add these vertices to existing vertices
  b.coordMap.innerCoords(1,1:nVertices,1,1) = f.patches{iFascicle}.vertices(:,1);
  b.coordMap.innerCoords(1,1:nVertices,1,2) = f.patches{iFascicle}.vertices(:,2);
  b.coordMap.innerCoords(1,1:nVertices,1,3) = f.patches{iFascicle}.vertices(:,3);
  b.coordMap.outerCoords = b.coordMap.innerCoords;
  % these are the display vertices which are the same as the coords
  b.coordMap.innerVtcs(1:nVertices,:) = f.patches{iFascicle}.vertices;
  b.coordMap.outerVtcs = b.coordMap.innerVtcs;
  % triangle faces
  b.coordMap.tris = f.patches{iFascicle}.faces;

  % reset the base and cache
  v = viewSet(v,'base',b,baseNum);
  v = viewSet(v,'baseCache','clear',b.name);
  v = viewSet(v,'overlayCache','clear',b.name);
  v = viewSet(v,'roiCache','clear',b.name);
  % and redisplay
  refreshMLRDisplay(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatomyFascicleMinmax    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatomyFascicleMinmax(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

keyboard
