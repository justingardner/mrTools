function mlrGuiSet(view,field,value,varargin)
%
%        $Id$
% mlrGuiSet(view,field,value);
%
% view can be either a view structure or a viewNum. Either way, sets
% handles in the global variable: guidata(MLR.views{viewNum}.figure).
% 
% djh 6/2005

mrGlobals;

if ieNotDefined('view'), mrErrorDlg('No view specified.'); end
if ieNotDefined('field'), mrErrorDlg('No parameter specified'); end
if ieNotDefined('value'), val = []; end

% viewNum can be either a view structure or a view number.
if isnumeric(view)
  viewNum = view;
  % *** Delete the following line after eliminating references to view ***
  view = MLR.views{viewNum};
elseif isview(view)
  viewNum = view.viewNum;
else
  return
end
if (viewNum < 1) | (viewNum > length(MLR.views))
  mrErrorDlg('Invalid viewNum,');
end

% Return if no gui
if isempty(MLR.views{viewNum}.figure)
  return
else
  % Get the gui handles to
  handles = guidata(MLR.views{viewNum}.figure);
end

switch lower(field)

 case {'grouppopup','groupstring'}
  % Set the groupPopup string array
  set(handles.groupPopup,'String',value);

 case {'group'}
  % Choose the group
  set(handles.groupPopup,'Value',value);

 case {'roipopup','roistring'}
  % Set the roiPopup string array
  set(handles.roiPopup,'String',value);

 case {'roi'}
  % Choose the roi
  if strcmp(get(handles.roiPopup,'style'),'popupmenu')
    value = value(1); %if this is not a listbox, we set only one ROI;
  end
  set(handles.roiPopup,'Value',value);

 case {'basepopup'}
  % Set the basePopup string array
  set(handles.basePopup,'String',value);

 case {'labelrois'}
  % mlrGuiSet(view,'labelrois',value);
  if value
    set(handles.labelsROIsMenuItem,'Checked','on');
  else
    set(handles.labelsROIsMenuItem,'Checked','off');
  end
 case {'showrois'}
  % mlrGuiSet(view,'showrois',value);
  
  % figure out which menu item should be checked
  onItem = find(strcmp(value,{'all','all perimeter','selected','selected perimeter','group','group perimeter','hide'}));
  onOrOff = {'off','off','off','off','off','off','off'};
  onOrOff{onItem} = 'on';
  % turn the check marks on/off
  set(handles.showAllMenuItem,'Checked',onOrOff{1});
  set(handles.showAllPerimeterMenuItem,'Checked',onOrOff{2});
  set(handles.showSelectedMenuItem,'Checked',onOrOff{3});
  set(handles.showSelectedPerimeterMenuItem,'Checked',onOrOff{4});
  set(handles.showGroupMenuItem,'Checked',onOrOff{5});
  set(handles.showGroupPerimeterMenuItem,'Checked',onOrOff{6});
  set(handles.hideROIsMenuItem,'Checked',onOrOff{7});
  
  if isfield(handles,'roiDisplayModePopup')
    switch(value)
      case {'group','selected'}
        roiDisplayMode = 1;
      case {'group perimeter','selected perimeter'}
        roiDisplayMode = 2;
      case {'all'}
        roiDisplayMode = 3;
      case {'all perimeter'}
        roiDisplayMode = 4;
      case 'hide'
        roiDisplayMode = 5;
    end
    set(handles.roiDisplayModePopup,'value',roiDisplayMode);
  end
  
 case {'basetype'}
  % mlrGuiSet(view,'baseType',value);
  % value = 0 for regular or 1 for flat
  if value == 0
    set(handles.sagittalRadioButton,'Visible','on');
    set(handles.coronalRadioButton,'Visible','on');
    set(handles.axialRadioButton,'Visible','on');
    set(handles.corticalDepth,'Visible','off');
    set(handles.corticalDepthSlider,'Visible','off');
    set(handles.corticalDepthText,'Visible','off');
    set(handles.flatViewerMenuItem,'Enable','off');
    set(handles.calcDistMenu, 'Enable', 'off');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','off');
    if isfield(handles,'corticalMaxDepthSlider')
      set(handles.corticalMaxDepthSlider,'Visible','off');
      set(handles.corticalMaxDepthText,'Visible','off');
      set(handles.linkMinMaxDepthCheck,'Visible','off');
    end
  elseif value >= 1
    set(handles.sagittalRadioButton,'Visible','off');
    set(handles.coronalRadioButton,'Visible','off');
    set(handles.axialRadioButton,'Visible','off');
    set(handles.corticalDepth,'Visible','on');
    set(handles.corticalDepthSlider,'Visible','on');
    set(handles.corticalDepthText,'Visible','on');
    set(handles.flatViewerMenuItem,'Enable','on');
    set(handles.flatViewerMenuItem,'Label','Flat Viewer');
    set(handles.calcDistMenu, 'Enable', 'on');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','on');
    if isfield(handles,'corticalMaxDepthSlider')
      set(handles.corticalMaxDepthSlider,'Visible','on');
      set(handles.corticalMaxDepthText,'Visible','on');
      set(handles.linkMinMaxDepthCheck,'Visible','on');
    end
  end
  if value == 2
    set(handles.createRoiMenu,'Enable','on');
      set(handles.createLineMenuItem,'Enable','off');
      set(handles.createRectangleMenuItem,'Enable','off');
      set(handles.addLineMenuItem,'Enable','off');
      set(handles.addRectangleMenuItem,'Enable','off');
      set(handles.removeLineMenuItem,'Enable','off');
      set(handles.removeRectangleMenuItem,'Enable','off');
    set(handles.addRoiMenu,'Enable','on');
    set(handles.removeRoiMenu,'Enable','on');
    set(handles.rotateSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'Visible','on');
    set(handles.baseTiltText,'Visible','on');
    set(handles.baseTilt,'Visible','on');
    set(handles.flatViewerMenuItem,'Enable','on');
    set(handles.flatViewerMenuItem,'Label','Surface Viewer');
    set(handles.calcDistMenu, 'Enable', 'off');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','on');
  else
    set(handles.baseTiltSlider,'Visible','off');
    set(handles.baseTiltText,'Visible','off');
    set(handles.baseTilt,'Visible','off');
    set(handles.createRoiMenu,'Enable','on');
    set(handles.addRoiMenu,'Enable','on');
      set(handles.createLineMenuItem,'Enable','on');
      set(handles.createRectangleMenuItem,'Enable','on');
      set(handles.addLineMenuItem,'Enable','on');
      set(handles.addRectangleMenuItem,'Enable','on');
      set(handles.removeLineMenuItem,'Enable','on');
      set(handles.removeRectangleMenuItem,'Enable','on');
    set(handles.removeRoiMenu,'Enable','on');
    set(handles.rotateSlider,'SliderStep',[1 45]./360);
  end		  
 case {'basevolume'}
  % Choose the baseVolume
  set(handles.basePopup,'Value',value);
  
 case {'basedims'}
  % mlrGuiSet(view,'baseDims',[ydim xdim zdim]);
  newDims = value;
  newCoords = min(handles.coords,newDims);
  handles.coords = min(handles.coords,newCoords);

 case {'basegamma'}
  % mlrGuiSet(view,'baseGamma',value);
  set(handles.baseGammaSlider,'Value',value);
  set(handles.baseGammaText,'String',thisNum2str(value));
  
 case {'basetilt'}
  % mlrGuiSet(view,'baseMax',value);
  set(handles.baseTiltSlider,'Value',value);
  set(handles.baseTiltText,'String',thisNum2str(value));

 case {'analysispopup'}
  % mlrGuiSet(view,'analysisPopup',strings);
  set(handles.analysisPopup,'String',value);

 case {'analysis'}
  % mlrGuiSet(view,'analysis',overlayNum);
  set(handles.analysisPopup,'Value',value);
  
 case {'overlaypopup'}
  % mlrGuiSet(view,'overlayPopup',strings);
  % mlrGuiSet(view,'overlayPopup',strings,overlayList); %optional argument overlayList indicates that this is a subset of strings to change
  if ~strcmp(value,'none') 
    %identify overlays that have been masked by putting a bullet before their name
    epsilon = 1e-7; %value differing by less than epsilon are considered equal
    if ieNotDefined('varargin')
      overlayList = 1:length(value);
    else
      overlayList = varargin{1};
      newStrings=  value;
      value = get(handles.overlayPopup,'String');
      value(overlayList) = newStrings;
    end
    for iOverlay = overlayList
      clip = viewGet(view,'overlayclip',iOverlay);
      minOverlayData = viewGet(view,'minoverlaydata',iOverlay);
      maxOverlayData = viewGet(view,'maxoverlaydata',iOverlay);
      if (~isempty(minOverlayData) && (clip(1)-minOverlayData)>epsilon) ||...
            (~isempty(maxOverlayData) && (maxOverlayData-clip(2))>epsilon) || ...
            clip(1)==clip(2) %if min and max clip values are equal, the whole overlay will be masked
         value{iOverlay} = [char(42) ' ' value{iOverlay}];
      else
         value{iOverlay} = ['  ' value{iOverlay}];
      end
    end
  else
    set(handles.overlayPopup,'value',1);
  end
  set(handles.overlayPopup,'String',value);

 case {'overlay'}
  % mlrGuiSet(view,'overlay',overlayNum);
  if strcmp(get(handles.overlayPopup,'style'),'popupmenu')
    value = value(1); %if this is not a listbox, we set only one overlay;
  end
  set(handles.overlayPopup,'Value',value);
    
  if length(value)==1
    set(handles.overlayMinSlider,'enable','on')
    set(handles.overlayMinText,'enable','on')
    set(handles.overlayMaxSlider,'enable','on')
    set(handles.overlayMaxText,'enable','on')
    if isfield(handles,'clippingOverlaysListbox')
      set(handles.clippingOverlaysListbox,'enable','on');
    end
%     set(handles.alphaText,'enable','on');
%     set(handles.alphaSlider,'enable','on');
  else %if more than one overlay selected, disable overlay controls
    if isfield(handles,'clippingOverlaysListbox')
      set(handles.clippingOverlaysListbox,'enable','off');
%       set(handles.clippingOverlaysListbox,'String',[],'enable','off');
    end
    set(handles.overlayMinSlider,'enable','off')
    set(handles.overlayMinText,'enable','off')
    set(handles.overlayMaxSlider,'enable','off')
    set(handles.overlayMaxText,'enable','off')
%     set(handles.alphaText,'enable','off');
%     set(handles.alphaSlider,'enable','off');
  end
  
 case {'overlaymin'}
  % mlrGuiSet(view,'overlayMin',value);
  if rem(value,1)~=0 %if the value is not an integer
    value = floor(double(value)*1e6)/1e6; %round it down
  end
  value = clipToSlider(handles.overlayMinSlider,value);
  set(handles.overlayMinSlider,'Value',value);
  set(handles.overlayMinText,'String',thisNum2str(value)); 

 case {'overlayminrange'}
  % mlrGuiSet(view,'overlayMinRange',[min,max]);
  if ~all(isfinite(value))
    mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Min Slider because overlay range is not finite [%f %f]',value(1),value(2)));
    set(handles.overlayMinSlider,'Min',max(value(1),-realmax),'Max',min(value(2),realmax),'visible','off');
  elseif value(2) < value(1)
    mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Min Slider because overlay range is not increasing [%f > %f]',value(1),value(2)));
    set(handles.overlayMinSlider,'Min',max(value(1),-realmax),'Max',min(value(2),realmax),'visible','off');
  else
    set(handles.overlayMinSlider,'Min',value(1),'Max',value(2),'visible','on');
  end

 case {'overlaymax'}
  % mlrGuiSet(view,'overlayMax',value);
  if rem(value,1)~=0 %if the value is not an integer
    value = ceil(double(value)*1e6)/1e6; %round it up
  end
  value = clipToSlider(handles.overlayMaxSlider,value);
  set(handles.overlayMaxSlider,'Value',value);
  set(handles.overlayMaxText,'String',thisNum2str(value)); 

 case {'overlaymaxrange'}
  % mlrGuiSet(view,'overlayMinRange',[min,max]);
  if ~all(isfinite(value))
    mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Max Slider because overlay range is not finite [%f %f]',value(1),value(2)));
    set(handles.overlayMaxSlider,'Min',max(value(1),-realmax),'Max',min(value(2),realmax),'visible','off');
  elseif value(2) < value(1)
    mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Max Slider because overlay range is not increasing [%f > %f]',value(1),value(2)));
    set(handles.overlayMaxSlider,'Min',max(value(1),-realmax),'Max',min(value(2),realmax),'visible','off');
  else
    set(handles.overlayMaxSlider,'Min',value(1),'Max',value(2),'visible','on');
  end

 case {'clippingoverlays'}
  % mlrGuiSet(view,'clippingOverlays',overlayList);
  if isfield(handles,'clippingOverlaysListbox') 
    overlayStrings=get(handles.overlayPopup,'String');
    selected=get(handles.clippingOverlaysListbox,'value');
    if selected>length(value)
      selected=1;
    end
    set(handles.clippingOverlaysListbox,'String',overlayStrings(value),'value',selected);
  end

 case {'alpha','overlayalpha'}
  % mlrGuiSet(view,'alpha',value);
  value = clipToSlider(handles.alphaSlider,value);
  set(handles.alphaSlider,'Value',value);
  set(handles.alphaText,'String',thisNum2str(value));
  set(handles.alphaSlider,'sliderStep',[0.1 0.5]);

 case {'nscans'}
  % mlrGuiSet(view,'nscans',value);
  nScans = round(value);
  curScan = round(get(handles.scanSlider,'Value'));
  if (nScans > 1)
    set(handles.scanSlider,'Min',1);
    set(handles.scanSlider,'Max',nScans);
    set(handles.scanSlider,'sliderStep',[1/(nScans-1) 1/(nScans-1)]);          
    set(handles.scanSlider,'Visible','on');
    curScan = min(curScan,nScans);
  else
    set(handles.scanSlider,'Min',0.9);
    set(handles.scanSlider,'Max',1.1);
    set(handles.scanSlider,'Visible','off');
    curScan = 1;
  end
  mlrGuiSet(view,'scan',curScan);

 case {'scan'}
  % mlrGuiSet(view,'scan',value);
  value = clipToSlider(handles.scanSlider,value,1);
  set(handles.scanSlider,'Value',value);
  set(handles.scanText,'String',num2str(value));
  % description
  description = viewGet(view,'description',value);
  set(viewGet(view,'figNum'),'Name',sprintf('%s: %s',getLastDir(MLR.homeDir),description));

 case {'scantext'}
  % mlrGuiSet(view,'scanText',value);
  set(handles.scanText,'String',num2str(value));
  % description
  description = viewGet(view,'description',value);
  set(viewGet(view,'figNum'),'Name',sprintf('%s: %s',getLastDir(MLR.homeDir),description));
  
 case {'nslices'}
  % mlrGuiSet(view,'nslices',value);
  value = round(value);
  if (value > 1)
    % reset range
    set(handles.sliceSlider,'Min',1);
    set(handles.sliceSlider,'Max',value);
    % reset step
    set(handles.sliceSlider,'sliderStep',[1/(value-1),1/(value-1)]);
    set(handles.sliceSlider,'Visible','on');
  else
    set(handles.sliceSlider,'Visible','off');
    set(handles.sliceText,'String',0);
  end

 case {'slice'}
  % mlrGuiSet(view,'slice',value);
  value = clipToSlider(handles.sliceSlider,value,1);
  if ~isempty(value)
    handles.coords(handles.sliceOrientation) = value;
    set(handles.sliceSlider,'Value',value);
    set(handles.sliceText,'String',num2str(value));
  end
 case {'slicetext'}
  % mlrGuiSet(view,'sliceText',value);
  handles.coords(handles.sliceOrientation) = value;
  set(handles.sliceText,'String',num2str(value));
  view = viewSet(view,'curSlice',value);

 case {'corticaldepth'} %this sets a single cortical depth, whether there are one or two sliders
  % mlrGuiSet(view,'corticalDepth',value);
   mlrGuiSet(view,'corticalMinDepth',value);
   if isfield(handles,'linkMinMaxDepthCheck') 
     set(handles.linkMinMaxDepthCheck,'value',true)
   end
   if isfield(handles,'corticalMaxDepthSlider')
     mlrGuiSet(view,'corticalMaxDepth',value);
   end
  
 case {'corticalmindepth'}
  % mlrGuiSet(view,'corticalDepth',value);
  value = min(value,1);value = max(value,0);
  corticalDepthBins = viewGet(view,'corticalDepthBins');
  value = round(value*(corticalDepthBins-1))/(corticalDepthBins-1);
  set(handles.corticalDepthSlider,'Value',value);
  set(handles.corticalDepthText,'String',num2str(value));
  
 case {'corticalmaxdepth'}
  % mlrGuiSet(view,'corticalDepth',value);
  value = min(value,1);value = max(value,0);
  corticalDepthBins = viewGet(view,'corticalDepthBins');
  value = round(value*(corticalDepthBins-1))/(corticalDepthBins-1);
  set(handles.corticalMaxDepthSlider,'Value',value);
  set(handles.corticalMaxDepthText,'String',num2str(value));
  
 case {'sliceorientation'}
  % mlrGuiSet(view,'sliceorientation',value);
  sliceOrientation = value;
  handles.sliceOrientation = sliceOrientation;
  % Set the correct radio button and unset the other two
  switch sliceOrientation
    % axial
   case 1
    % axial
    set(handles.sagittalRadioButton,'Value',0);
    set(handles.coronalRadioButton,'Value',0);
    set(handles.axialRadioButton,'Value',1);
   case 2
    % coronal
    set(handles.sagittalRadioButton,'Value',0);
    set(handles.coronalRadioButton,'Value',1);
    set(handles.axialRadioButton,'Value',0);
   case 3
    % sagittal
    set(handles.sagittalRadioButton,'Value',1);
    set(handles.coronalRadioButton,'Value',0);
    set(handles.axialRadioButton,'Value',0);
  end
 case {'rotate'}
  % mlrGuiSet(view,'rotate',value);
  value = clipToSlider(handles.rotateSlider,value);
  set(handles.rotateText,'String',num2str(value));
  set(handles.rotateSlider,'Value',value);
  viewSet(view,'rotate',value);
 case {'viewnum'}
  % mlrGuiSet(view,'viewnum',value);
  handles.viewNum = value;

 otherwise
  error(['Invalid field: ',field]);
end
guidata(MLR.views{viewNum}.figure,handles);


function value = clipToSlider(slider,value,integerFlag)
% Clips value so that it doesn't exceed slider limits.
% slider is a slider handle
% value must be a number (otherwise, use current slider value)
% integerFlag forces value to be an integer
if ieNotDefined('integerFlag')
  integerFlag = 0;
end
if ~isnumeric(value)
  value = get(slider,'Value');
end
if integerFlag
  if (value < get(slider,'Min'))
    value = ceil(get(slider,'Min'));
  end
  if (value > get(slider,'Max'))
    value = floor(get(slider,'Max'));
  end
else
  if (value < get(slider,'Min'))
    value = get(slider,'Min');
  end
  if (value > get(slider,'Max'))
    value = get(slider,'Max');
  end
end

%modified num2str to increase the number of decimals for reals
function value = thisNum2str(value)

  if rem(value,1)~=0
    value = num2str(value,'%.6f');
  else
    value = num2str(value);
  end

