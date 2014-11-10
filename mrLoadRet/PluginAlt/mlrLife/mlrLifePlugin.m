% mlrLifePlugin
%
%        $Id:$ 
%      usage: mlrLifePlugin(action,<v>)
%         by: justin gardner & franco pestilli
%       date: 09/09/2014
%    purpose: Plugin function for LiFE
%
function retval = mlrLifePlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrLifePlugin) Need a valid view to install plugin'));
  else
    % add the menu item for LiFE - this one should open Dt6 file
    mlrAdjustGUI(v,'add','menu','LiFE','/File/ROI','Callback',@mlrLifeFileOpen);
    % this menu item allows importation of fascicles as a surface
    mlrAdjustGUI(v,'add','menu','Import fascicles','/File/Base anatomy/Import surface','Callback',@mlrLifeImportFascicles);
    % now create a panel - this will be used for adding some UI
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
    mlrAdjustGUI(v,'add','control','multiBaseOverlayText','panel','Multiple base display','style','text','position', [0.01    0.68    0.2   0.07 ],'String','Overlay');
    mlrAdjustGUI(v,'add','control','multiBaseOverlay','panel','Multiple base display','style','popup','position', [0.22    0.68    0.72   0.07 ],'String',color2RGB,'Callback',@multiBaseOverlay);
    % add color alpha
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaText','panel','Multiple base display','style','text','position', [0.01    0.6    0.2   0.07 ],'String','Overlay Alpha');
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaSlider','panel','Multiple base display','style','slider','position', [0.22    0.6    0.57   0.07 ],'String','Alpha','SliderStep',[0.1 0.25],'Callback',@multiBaseOverlayAlpha);
    mlrAdjustGUI(v,'add','control','multiBaseOverlayAlphaEdit','panel','Multiple base display','style','edit','position', [0.8    0.6    0.19   0.07 ],'Callback',@multiBaseOverlayAlpha);
    

    % ROI controls
    mlrAdjustGUI(v,'add','control','roiBaseListBox','panel','Multiple base display','style','listbox','position', [0.02    0.1    0.96   0.38 ],'Callback',@roiListboxSelect,'Max',2);

    % add the callback that will tell the above listbox when new
    % bases have been added
    v = viewSet(v,'callback','baseChange',@mlrLifeBaseChange);
    % also register a change when someone switches the curBase
    v = viewSet(v,'callback','curBaseChange',@mlrLifeBaseChange);
    %mlrAdjustGUI(v,'add','panel','Multiple ROI display',.25);


    % This is a command that could be used to install some default interrogators
    %mlrAdjustGUI(v,'add','interrogator',{'makeTract'});

    % This is a command that could be used to install some default colormaps
    % that will show up when you do /Edit/Overlay
    %mlrAdjustGUI(v,'add','colormap','gray');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(mlrLifePlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeFileOpen    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeFileOpen(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% bring up file dialog to select dt6 file
startPathStr = pwd;
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% now bring up selection dialog
paramsInfo = {{'groupName','DTI'},...
	      {'analysisName','DTI'}};
params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

% convert t1 into a MLR anatomy
% Set required structure fields (additional fields are set to default
% values when viewSet calls isbase).
base.name = stripext(stripext(getLastDir(dt.files.t1)));
base.data = double(t1.img);
% create a header
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;

% set nifti header and permutation matrix in base structure
base.hdr = mlrImageGetNiftiHeader(h);
base.permutationMatrix = getPermutationMatrix(base.hdr);

% make it a full base
[tf base] = isbase(base);

% add it to the view
v = viewSet(v,'newBase',base);

% make group
v = viewSet(v,'newGroup',params.groupName);
v = viewSet(v,'curGroup',params.groupName);

% make nifti header for DTI
h = [];hdr = [];
h.dim = size(dt.dt6);
[tf h] = mlrImageIsHeader(h);
h.pixdim = dt.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = dt.xformToAcpc;
hdr = mlrImageGetNiftiHeader(h);

% save tseries
saveNewTSeries(v,dt.dt6,[],hdr);
scanNum = viewGet(v,'nScans');

% get the overlay
refreshMLRDisplay(v);

% convert RGB into a colormap
% make table
numValues = 8;
values = 0:(256/(numValues-1)):256;
rgbTable = [];
for rValue = 1:length(values)
  for gValue = 1:length(values)
    for bValue = 1:length(values)
      rgbTable(end+1,:) = [values(rValue) values(gValue) values(bValue)];
    end
  end
end
rgbTable = rgbTable(1:2:end,:);

% convert overlay to this table
disppercent(-inf,'Converting vectorRGB to index');
for x = 1:size(o.vectorRGB,1)
  disppercent(x/size(o.vectorRGB,1));
  for y = 1:size(o.vectorRGB,2)
    for z = 1:size(o.vectorRGB,3)
      colorVal = double(squeeze(o.vectorRGB(x,y,z,:)));
      [dummy colorIndex] = min(sum((rgbTable - repmat(colorVal',size(rgbTable,1),1)).^2,2));
      colormapRGB(x,y,z) = colorIndex;
    end
  end
end
disppercent(inf);

% display the overlay
rgbTable = rgbTable/256.0;
% FIX, FIX, FIX Create a proper analysis structure with all the overlays.
mrDispOverlay(colormapRGB,scanNum,params.groupName,v,'overlayName=principleDiffusionDirection','saveName=diffusionAnalysis','cmap',rgbTable);
a = viewGet(v,'curAnalysis');
mrDispOverlay(double(o.b0),scanNum,params.groupName,v,'overlayName=b0','saveName=diffusionAnalysis');

dateString = datestr(now);
b0.name = 'b)';
b0.groupName = params.groupName;
b0.function = '';
b0.reconcileFunction = 'defaultReconcileParams';
b0.data = double(o.b0);
b0.date = dateString;
b0.params = cell(1,viewGet(v,'nScans'));
b0.range = [0 1];
b0.clip = [0 1];
% colormap is made with a little bit less on the dark end
b0.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';
r2.mergeFunction = 'defaultMergeParams';


% create analysis
dti.name = params.analysisName;
dti.type = 'dti';
dti.groupName = params.groupName;
dti.function = '';
dti.reconcileFunction = 'defaultReconcileParams';
dti.mergeFunction = 'defaultMergeParams';
dti.guiFunction = '';
dti.params = params;
dti.overlays = b0;
dti.curOverlay = 1;
dti.date = dateString;
v = viewSet(v,'newAnalysis',dti);

% Save it
saveAnalysis(v,dti.name);


refreshMLRDisplay(v);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeImportFascicles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeImportFascicles(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Load the fascicles from disk
[fg, fgFilename] = fgRead;
if isempty(fg), return,end

% bring up file dialog to select dt6 file
% this is so that we can get the T1 file
% FIX, FIX, FIX - should find some way to automatically
% load the T1 file without having the user have to go find it.
startPathStr = fullfile(fileparts(fgFilename),'..');
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% Get the header from the T1 so that we have the correct sform etc
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;

% convert the X, Y, Z coordinates which are in AC/PC back to 3D image coordinates
huh = h.sform;
h.sform = eye(4);
%huh(1,1) = 1*t1.mmPerVoxel(1);
%huh(2,2) = 1*t1.mmPerVoxel(2);
%huh(3,3) = 1*t1.mmPerVoxel(3);
%huh(1,4) = huh(1,4)*t1.mmPerVoxel(1);
%huh(2,4) = huh(2,4)*t1.mmPerVoxel(2);
%huh(3,4) = huh(3,4)*t1.mmPerVoxel(3);
%scaleXform = diag([t1.mmPerVoxel 1]);
%fg = dtiXformFiberCoords(fg,scaleXform);
%huh = huh*inv(scaleXform);
%fg = dtiXformFiberCoords(fg,inv(huh));

% Build frames from the fascicles
[X, Y, Z] = mbaBuildFascicleFrame(fg.fibers);
%[Y, X, Z] = mbaBuildFascicleFrame(fg.fibers);

% number of fascicles
nFascicles = length(X);

% Build a patch from the frame
nTotalVertices = 0;nTotalTris = 0;
for i = 1:nFascicles
  fasciclePatches{i} = surf2patch(X{i},Y{i},Z{i},'triangles');
  % compute how many vertices and tris we have all together
  nTotalVertices = size(fasciclePatches{i}.vertices,1) + nTotalVertices;
  nTotalTris = size(fasciclePatches{i}.faces,1) + nTotalTris;
end

% create an MLR base structure
fascicleBase.name = fixBadChars(fg.name);

% set this to be a surface
fascicleBase.type = 2; 

% set nifti header and permutation matrix in base structure
fascicleBase.hdr = mlrImageGetNiftiHeader(h);
fascicleBase.permutationMatrix = getPermutationMatrix(fascicleBase.hdr);
fascicleBase.vol2mag = fascicleBase.hdr.sform44;

% set path and names of files
fascicleBase.coordMap.path = fileparts(fgFilename);
fascicleBase.coordMap.innerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.innerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.curvFileName = '';
fascicleBase.coordMap.anatFileName = dt.files.t1;

% initalize fields
fascicleBase.data = [];
% set dimensions for coordMap
fascicleBase.coordMap.dims = h.dim;
fascicleBase.coordMap.innerCoords = nan(1,nTotalVertices,1,3);
fascicleBase.coordMap.innerVtcs = nan(nTotalVertices,3);
fascicleBase.coordMap.tris = nan(nTotalTris,3);
nRunningTotalVertices = 0;
nRunningTotalTris = 0;

% now put all fascicles vertices and triangles into one coordMap
disppercent(-inf,sprintf('(mlrLifePlugin) Converting %i fascicles',nFascicles));
for iFascicle = 1:nFascicles
  % number of vertices and triangles
  nVertices = size(fasciclePatches{iFascicle}.vertices,1);
  nTris = size(fasciclePatches{iFascicle}.faces,1);
  % the data which is the grayscale value to color the fascicles with (rand for now)
  fascicleBase.data = [fascicleBase.data rand(1,nVertices)];
  % convert vertices to a coord map which has one x,y,z element for each possible
  % location on the surface (which actually is just a 1xnVerticesx1 image)
  % add these vertices to existing vertices
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,1) = fasciclePatches{iFascicle}.vertices(:,1);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,2) = fasciclePatches{iFascicle}.vertices(:,2);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,3) = fasciclePatches{iFascicle}.vertices(:,3);
  % these are the display vertices which are the same as the coords
  fascicleBase.coordMap.innerVtcs(nRunningTotalVertices+1:nRunningTotalVertices+nVertices,:) = fasciclePatches{iFascicle}.vertices;
  % triangle faces
  fascicleBase.coordMap.tris(nRunningTotalTris+1:nRunningTotalTris+nTris,:) = (fasciclePatches{iFascicle}.faces + nRunningTotalVertices);
  % update runing totals
  nRunningTotalVertices = nRunningTotalVertices + nVertices;
  nRunningTotalTris= nRunningTotalTris + nTris;
  disppercent(iFascicle/nFascicles);
end
disppercent(inf);

% copy the inner to outer since they are all the same for fascicles
fascicleBase.coordMap.outerCoords = fascicleBase.coordMap.innerCoords;
fascicleBase.coordMap.outerVtcs = fascicleBase.coordMap.innerVtcs;

% make it a full base
[tf fascicleBase] = isbase(fascicleBase);

% add it to the view
v = viewSet(v,'newBase',fascicleBase);

refreshMLRDisplay(v);

keyboard

%%%%%%%%%%%%%%%%%%%%
%    baseChange    %
%%%%%%%%%%%%%%%%%%%%
function v = mlrLifeBaseChange(v)

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
    baseMatch = find(strcmp(base.overlay,baseColorNames));
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
