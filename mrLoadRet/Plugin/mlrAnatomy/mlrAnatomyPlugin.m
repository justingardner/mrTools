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

    % add a menu item to import rois from freesurfer
    mlrAdjustGUI(v,'add','menu','Import Freesurfer Label','/File/ROI/Import','Callback',@mlrAnatomyImportFreesurferLabel);

    % add a menu item to make a planer base anatomy
    mlrAdjustGUI(v,'add','menu','Make plane base','/File/Base anatomy/Import surface','Callback',@mlrAnatomyMakePlaneBase,'Separator','on');
    
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

disppercent(-inf,'(mlrAnatomyMakePlaneBase) Making plane');

% get the base that we are going to make this plane from
canonical = viewGet(v,'base',viewGet(v,'baseNum',params.baseName));
canonicalDims = size(canonical.data);

% copy info into to create this new base
b = canonical;
b.coordMap.dims = canonicalDims;
b.coordMap.innerSurfaceFileName = '';
b.coordMap.innerCoordsFileName = '';
b.coordMap.outerSurfaceFileName = '';
b.coordMap.outerCoordsFileName = '';
b.coordMap.curvFileName = '';
b.coordMap.anatFileName = canonical.name;
b.coordMap.path = '';

% get center
xCenter = 0;yCenter = 00;zCenter = 0;

% make rotation matrix
params.xyRot = 0;
params.yzRot = 00;
params.xzRot = 0;
c = cos(d2r(params.xyRot));
s = sin(d2r(params.xyRot));
rotxy = [c -s 0 0;s  c 0 0;0  0 1 0;0  0 0 1];
c = cos(d2r(params.yzRot));
s = sin(d2r(params.yzRot));
rotyz = [1  0  0 0;0  c -s 0;0  s  c 0;0  0  0 1];
c = cos(d2r(params.xzRot));
s = sin(d2r(params.xzRot));
rotxz = [c  0 -s 0;0  1  0 0;s  0  c 0;0  0  0  1];
r = rotxy*rotyz*rotxz;

% now we make the surface positions
x = 1:canonicalDims(2);
y = 1:canonicalDims(1);
[x y] = meshgrid(x,y);
z = ones(size(x))*round(canonicalDims(3)/2);
x = x';y = y';z = z';

% now get coordinates (xy swaping because of image coordinates)
% also move center to user specified center
coords = [y(:)+xCenter x(:)+yCenter z(:)+zCenter];
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

% and put them into the coordMap strucutre
b.coordMap.coords(1,:,1,1) = coords(:,1);
b.coordMap.coords(1,:,1,2) = coords(:,2);
b.coordMap.coords(1,:,1,3) = coords(:,3);
b.coordMap.outerCoords = b.coordMap.coords;
b.coordMap.innerCoords = b.coordMap.coords;
b.coordMap.innerVtcs = coords;
b.coordMap.outerVtcs = coords;
b.type = 2;

% now make all the triangles
tris = [];triNum = 0;
for iRow = 1:canonicalDims(1)-1
  for iCol = 1:canonicalDims(2)-1
    triNum = triNum+1;
    tris(triNum,1) = (iRow-1)*canonicalDims(2) + iCol;
    tris(triNum,2) = (iRow-1)*canonicalDims(2) + iCol+1;
    tris(triNum,3) = iRow*canonicalDims(2) + iCol;
    triNum = triNum+1;
    tris(triNum,1) = (iRow-1)*canonicalDims(2) + iCol+1;
    tris(triNum,2) = iRow*canonicalDims(2) + iCol;
    tris(triNum,3) = iRow*canonicalDims(2) + iCol+1;
  end
  disppercent(iRow/canonicalDims(1));
end
b.coordMap.tris = tris;

% now get the data for the points
b.data = interp3(canonical.data,coords(:,2),coords(:,1),coords(:,3),'linear',0);
b.data = b.data(:)';
b.name = 'plane';
v = viewSet(v,'newBase',b);
refreshMLRDisplay(v);
disppercent(inf);

