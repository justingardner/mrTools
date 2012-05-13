% improvedGUIPlugin.m
%
%        $Id: improvedGUIPlugin.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: improvedGUIPlugin(action,<thisView>)
%         by: julien besle
%       date: 13/02/11
%    purpose: 
%

function retval = improvedGUIPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help improvedGUIPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(improvedGUIPlugin) Need a valid view to install plugin'));
  else
    
    % new uicontrols and reposition old ones
    if viewGet(thisView,'baseType')>0
      corticalDepthVisibility = 'on';
    else
      corticalDepthVisibility = 'off';
    end    
    controlFontSize=10;
    checkFontSize=11;
    popupFontSize = 11;
    %---------------------------- Group and Scan controls-----------------------------------------
    mlrAdjustGUI(thisView,'set','group','position',               [0.01    0.96    0.06   0.025]);
    mlrAdjustGUI(thisView,'set','group','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','groupPopup','position',          [0.08    0.95    0.2    0.035]);
    mlrAdjustGUI(thisView,'set','groupPopup','fontSize',popupFontSize);
    mlrAdjustGUI(thisView,'set','scan','position',                [0.01    0.915   0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','scan','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','scanSlider','position',          [0.07    0.92    0.16   0.02 ]);
    mlrAdjustGUI(thisView,'set','scanText','position',            [0.24    0.92    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','scanText','fontSize',controlFontSize);
    %---------------------------- Base controls  -------------------------------------------------
    mlrAdjustGUI(thisView,'set','baseImage','position',           [0.01    0.88    0.06   0.025]);
    mlrAdjustGUI(thisView,'set','baseImage','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','basePopup','position',           [0.07    0.87    0.21   0.035]);
    mlrAdjustGUI(thisView,'set','basePopup','fontSize',popupFontSize);
    mlrAdjustGUI(thisView,'set','slice','position',               [0.02    0.84    0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','slice','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','sliceSlider','position',         [0.08    0.845   0.15   0.02 ]);
    mlrAdjustGUI(thisView,'set','sliceText','position',           [0.24    0.845   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','sliceText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepth','position',       [0.02    0.81    0.07   0.06 ]);
    mlrAdjustGUI(thisView,'set','corticalDepth','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','position', [0.095   0.845   0.135  0.02 ]);
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','SliderStep',min(1/viewGet(thisView,'corticalDepthBins')*[1 3],1));
    mlrAdjustGUI(thisView,'set','corticalDepthText','position',   [0.24    0.845   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','corticalDepthText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepthText','BackgroundColor', get(0,'defaultUicontrolBackgroundColor'));
    corticalDepthBins=viewGet(thisView,'corticaldepthbins');
    corticalDepth=round((corticalDepthBins-1)/2)/(corticalDepthBins-1);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthSlider',...
        'SliderStep',min(1/(viewGet(thisView,'corticalDepthBins')-1)*[1 3],1),...
        'Callback',@corticalMaxDepthSlider_Callback,'visible',corticalDepthVisibility,...
                          'value',corticalDepth,'style','slider','position', [0.095   0.825    0.135  0.02 ]);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthText',...
        'Callback',@corticalMaxDepthText_Callback,'visible',corticalDepthVisibility,'String',num2str(corticalDepth),...
        'fontSize',controlFontSize,'BackgroundColor', get(0,'defaultUicontrolBackgroundColor'),...
                                        'style','edit','position',[0.24    0.825    0.04   0.025 ]);
    mlrGuiSet(thisView,'corticalmindepth',corticalDepth);
    mlrAdjustGUI(thisView,'add','control','linkMinMaxDepthCheck','style','checkbox','value',1,...
        'fontSize', checkFontSize,'visible',corticalDepthVisibility,...
        'String','Fix Depth Range','position',                    [0.095   0.79   0.23   0.025]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','position', [0.01    0.81    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','axialRadioButton','position',    [0.11    0.81    0.07   0.025]);
    mlrAdjustGUI(thisView,'set','axialRadioButton','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','position',  [0.19    0.81    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','baseGamma','string','Gamma');
    mlrAdjustGUI(thisView,'set','baseGamma','position',           [0.02    0.755   0.07   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseGamma','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','baseGammaSlider','position',     [0.10    0.76    0.13   0.02 ]);
    mlrAdjustGUI(thisView,'set','baseGammaText','position',       [0.24    0.76    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','baseGammaText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','rotate','position',              [0.02    0.725   0.06   0.03 ]);
    mlrAdjustGUI(thisView,'set','rotate','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','rotateSlider','position',        [0.09    0.73    0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','rotateText','position',          [0.24    0.73    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','rotateText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','baseTilt','position',            [0.02    0.695   0.03   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseTilt','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','baseTiltSlider','position',      [0.055   0.7     0.175  0.02 ]);
    mlrAdjustGUI(thisView,'set','baseTiltText','position',        [0.24    0.7     0.04   0.025]);
    mlrAdjustGUI(thisView,'set','baseTiltText','fontSize',controlFontSize);
    
    %---------------------------- ROI controls  --------------------------------------------------
    mlrAdjustGUI(thisView,'set','ROI','position',                 [0.01    0.66    0.08    0.025]);
    mlrAdjustGUI(thisView,'set','ROI','string','ROIs');
    mlrAdjustGUI(thisView,'set','ROI','fontWeight','bold');
    mlrAdjustGUI(thisView,'add','control','displayROILabels','style','checkbox','value',viewGet(thisView,'labelROIs'),...
        'callback',@displayROILabels_Callback,'fontSize', checkFontSize,'String',...
        'Display Labels','position',                              [0.09    0.66    0.18   0.025]);
    switch(viewGet(thisView,'showROIs'))
      case {'group','selected'}
        thisView=viewSet(thisView,'showROIs','selected');
        roiDisplayMode = 1;
      case {'group perimeter','selected perimeter'}
        roiDisplayMode = 2;
        thisView=viewSet(thisView,'showROIs','selected perimeter');
      case {'all'}
        thisView=viewSet(thisView,'showROIs','selected');
        roiDisplayMode = 3;
      case {'all perimeter'}
        roiDisplayMode = 4;
        thisView=viewSet(thisView,'showROIs','selected perimeter');
      case 'hide'
        roiDisplayMode = 5;
    end
    mlrAdjustGUI(thisView,'add','control','roiDislayMode','style','text','fontSize', checkFontSize,...
        'String','Display mode','position',                       [0.02    0.63    0.12   0.025]);
    mlrAdjustGUI(thisView,'add','control','roiDisplayModePopup','style','popupmenu','value',roiDisplayMode,...
        'callback',@roiDislayModePopup_Callback,'String',{'Selected ROIs (Voxels)' 'Selected ROIs (Perimeter)' 'All ROIs (Voxels)' 'All ROIs (Perimeter)','Hide'},...
        'fontSize', popupFontSize,'position',                     [0.14    0.625   0.14   0.035]);
    mlrAdjustGUI(thisView,'set','roiPopup','position',            [0.02    0.485   0.26   0.145]);
    mlrAdjustGUI(thisView,'set','roiPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','roiPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','roiPopup','callback',@roiList_Callback);
    mlrAdjustGUI(thisView,'set','roiPopup','value',viewGet(thisView,'currentRoi'));
    mlrAdjustGUI(thisView,'set','roiPopup','fontSize',controlFontSize);
    
    %---------------------------- Analysis and overlays controls  --------------------------------
    mlrAdjustGUI(thisView,'set','analysis','position',            [0.01    0.45    0.08   0.025]);
    mlrAdjustGUI(thisView,'set','analysis','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','analysisPopup','position',       [0.09    0.44    0.19   0.035]);
    mlrAdjustGUI(thisView,'set','analysisPopup','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','overlay','position',             [0.01    0.42    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','overlay','string','Overlays');
    mlrAdjustGUI(thisView,'set','overlay','fontWeight','bold');
    mlrAdjustGUI(thisView,'add','control','clipAcrossOverlays','style','checkbox','value',0,...
        'callback',@clipAcrossOverlays_Callback,'fontSize', checkFontSize,...
        'String','Clip across overlays','position',               [0.11    0.415   0.18   0.025]);
    mlrAdjustGUI(thisView,'add','control','colorBlending','style','text','fontSize', checkFontSize,...
        'String','Blending mode','position',                      [0.02    0.385   0.12   0.025]);
    mlrAdjustGUI(thisView,'add','control','colorBlendingPopup','style','popupmenu','value',2,...
        'callback',@colorBlendingPopup_Callback,'String',{'Alpha blend' 'Additive'},...
        'fontSize', popupFontSize,'position',                     [0.14    0.38    0.14   0.035]);
    mlrAdjustGUI(thisView,'set','overlayPopup','position',        [0.02    0.21    0.26   0.17 ]);
    mlrAdjustGUI(thisView,'set','overlayPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','overlayPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','overlayPopup','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','overlayMin','position',          [0.02    0.175   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','overlayMin','string','Min');
    mlrAdjustGUI(thisView,'set','overlayMin','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','overlayMinText','position',      [0.055   0.18    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','overlayMinText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','mapMax','position',              [0.15    0.175   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','mapMax','string','Max');
    mlrAdjustGUI(thisView,'set','mapMax','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','overlayMaxText','position',      [0.19    0.18    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','overlayMaxText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','overlayMinSlider','position',    [0.02    0.15   0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','overlayMaxSlider','position',    [0.02    0.13    0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','alpha','position',               [0.02    0.095   0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','alpha','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','alphaSlider','position',         [0.08    0.1     0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','alphaText','position',           [0.23    0.1     0.05   0.025]);
    mlrAdjustGUI(thisView,'set','alphaText','fontSize',controlFontSize);
   
    mlrAdjustGUI(thisView,'set','colorbar','position',            [0.35    0.1     0.58   0.07 ]);
    mlrAdjustGUI(thisView,'set','colorbar','fontSize',10);
    mlrAdjustGUI(thisView,'add','axes','colorbarRightBorder',...
      'YaxisLocation','right','XTick',[],'box','off','position',  [0.929   0.1     0.001  0.07 ]);

    % Add some useful menus
    mlrAdjustGUI(thisView,'add','menu','Unlink Stimfile','/Edit/Scan/Link Stimfile','callback',@unlinkStimfileMenuItem_Callback,'tag','unlinkStimfileMenuItem');
    mlrAdjustGUI(thisView,'set','/Edit/Scan/Link Stimfile','separator','on');
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','callback',@exportImage_Callback);
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','label','Export Images');
    mlrAdjustGUI(thisView,'add','menu','Apply MotionComp Transforms','/Analysis/Motion Compensation/Slice Time Correction (only)','callback',@applyMotionCompTransformsCallBack,'tag','applyMotionCompTransformMenuItem');
    mlrAdjustGUI(thisView,'add','menu','Copy sform','/Edit/Base Anatomy/Transforms/','callback',@copyBaseSformCallBack,'tag','copyBaseSformMenuItem');
    mlrAdjustGUI(thisView,'add','menu','Paste sform','/Edit/Base Anatomy/Transforms/','callback',@pasteBaseSformCallBack,'tag','pasteBaseSformMenuItem');
    
    %remove show ROI menus and re-arrange ROI menus
    mlrAdjustGUI(thisView,'set','findCurrentROIMenuItem','location','/ROI/Convert');
    mlrAdjustGUI(thisView,'set','findCurrentROIMenuItem','separator','off');
    mlrAdjustGUI(thisView,'remove','menu','showRoiMenu');
    mlrAdjustGUI(thisView,'set','undoRoiMenuItem','location','/ROI/Restrict');
    mlrAdjustGUI(thisView,'set','convertRoiMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','deleteRoiMenu','separator','on');
    mlrAdjustGUI(thisView,'set','restrictRoiMenuItem','label','Selected ROI(s)');
    mlrAdjustGUI(thisView,'set','deleteRoiMenuItem','label','Selected ROI(s)');
    mlrAdjustGUI(thisView,'set','saveROIMenuItem','label','Save selected ROI(s)');
    mlrAdjustGUI(thisView,'set','copyRoiMenuItem','label','Copy selected ROI(s)');
    mlrAdjustGUI(thisView,'set','pasteRoiMenuItem','label','Paste ROI(s)');
    mlrAdjustGUI(thisView,'set','editRoiMenuItem','label','Edit selected ROI(s)');
    
    %add 3D render viewer
    mlrAdjustGUI(thisView,'add','menu','3D Viewer','flatViewerMenuItem','callback',@renderIn3D,'tag','viewIn3DMenuItem');
    
    % return view
    retval = thisView;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Improvements include: superimposition of several overlays, average over a range of cortical depths, multiple ROI selection, 3D viewer ...';
 otherwise
   disp(sprintf('(improvedGUIPlugin) Unknown command %s',action));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function corticalMaxDepthSlider_Callback(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = get(hObject,'Value');
if get(handles.linkMinMaxDepthCheck,'value')
  minDepth = viewGet(viewNum,'corticalMinDepth');
  maxDepth = str2num(get(handles.corticalMaxDepthText,'string'));
  newMinValue = minDepth+newMaxValue-maxDepth;
  if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
    mlrGuiSet(viewNum,'corticalMinDepth',newMinValue);
    mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
    drawnow;
    refreshMLRDisplay(viewNum);
  else
    set(hObject,'Value',maxDepth);
  end
else
  mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
  refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
function corticalMaxDepthText_Callback(hObject,dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = str2num(get(hObject,'String'));
if isempty(newMaxValue) %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.corticalMaxDepthSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  if get(handles.linkMinMaxDepthCheck,'value')
    minDepth = viewGet(viewNum,'corticalMinDepth');
    maxDepth = get(handles.corticalMaxDepthSlider,'value');
    newMinValue = minDepth+newMaxValue-maxDepth;
    if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
      mlrGuiSet(viewNum,'corticalMinDepth',newMinValue);
      mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
      refreshMLRDisplay(viewNum);
    else
      set(hObject,'string',num2str(maxDepth));
    end
  else
    mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
    refreshMLRDisplay(viewNum);
  end
end

% --------------------------------------------------------------------
function clipAcrossOverlays_Callback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function colorBlendingPopup_Callback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
refreshMLRDisplay(viewNum);


% -------------------------------------------------------------------- Unlink stim menu
function unlinkStimfileMenuItem_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum, 'stimfilename', '');


% -------------------------------------------------------------------- ROI controls
function roiList_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;

roiNames = viewGet(viewNum,'roiNames');% I think we shouldn't have to get the names (but need to modify viewSet call)
% set the roi group
% viewSet(viewNum,'roiGroup',roiNames(get(hObject,'Value')));
viewSet(viewNum,'curROI',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

function roiDislayModePopup_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
switch(get(hObject,'value'))
  case 1
    viewSet(viewNum,'showROIs','selected');
  case 2
    viewSet(viewNum,'showROIs','selected perimeter');
  case 3
    viewSet(viewNum,'showROIs','group');
  case 4
    viewSet(viewNum,'showROIs','group perimeter');
  case 5
    viewSet(viewNum,'showROIs','hide');
end
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function displayROILabels_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum,'labelROIs',get(hObject,'Value'));
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function exportImage_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');

viewType = viewGet(thisView,'basetype');

if viewType<2
  nScans = viewGet(thisView,'nScans');
  paramsInfo = {};
  paramsInfo{end+1} = {'tifFileName', fullfile(viewGet(thisView,'homedir'),'Images/image.tif'),'Where to save the image montage file (TIF format)'};
  paramsInfo{end+1} = {'horizontalRange',[0 1],'minmax=[0 1]','Horizontal coordinates of the image cropped from the slice [left right] (normalized between 0 and 1)'};
  paramsInfo{end+1} = {'verticalRange',[0 1],'minmax=[0 1]','Vertical coordinates of the image cropped from the slice [top bottom] (normalized between 0 and 1)'};
  paramsInfo{end+1} = {'scanList', sprintf('[1:1:%d]',nScans), 'List of scans to export. Must be expressed in a string taht can be evaluated into a row vector'};

  if viewType
    paramsInfo{end+1} = {'depthMode', {'as displayed','one image per depth'}, '''As displayed'' = only one image is exported per scan, with the current depth settings, ''One image per depth''= exports one image per depth/scan within the displayed depth range'};
  else
    basedims = viewGet(thisView,'basedims');
    nSlices = basedims(viewGet(thisView,'basesliceindex'));

    paramsInfo{end+1} = {'sliceList', sprintf('[1:1:%d]',nSlices), 'List of slices to export. Must be expressed in a string that can be evaluated into a row vector'};
  end
  paramsInfo{end+1} = {'nRows', 1, 'minmax=[1 inf]','incdec=[-1 1]','Number of rows in the montage if multiple scans/slices/depths are exported'}; 
  paramsInfo{end+1} = {'exportImage', 0, 'type=pushbutton','','callback',@exportCallback,'callbackArg',thisView,'buttonString=Export Image','passParams=1','Previews and saves the montage without closing the menu'};

  % display dialog
  mrParamsDialog(paramsInfo,'Export Slices Parameters','modal=0');
else
    % if surface, would need to change the ouput of refreshMLRDisplay to get a 2D image rather than the coordinates of the surface...
    % but not sure whether some other parts of mrLoadRet use this output..
    mrWarnDlg('(improvedGUIPlugin) This function is not implemented for surfaces.');
end


% --------------------------------------------------------------------
function exportCallback(thisView,params)

%update the view
thisView = viewGet(thisView.viewNum,'view');

if viewGet(thisView,'basetype')
  switch(params.depthMode)
    case 'as displayed'
      sliceList = [];
    case 'one image per depth'
      depthBin = 1/(viewGet(thisView,'corticaldepthBins')-1);
      depth(1) = viewGet(thisView,'corticalmindepth');
      depth(2) = viewGet(thisView,'corticalmaxdepth');
      depth = sort(depth);
      sliceList = depth(1):depthBin:depth(2);
%       [~,sliceList] = ismember(round(1e5*(depth(1):depthBin:depth(2)))/1e5,round(1e5*(0:depthBin:1))/1e5);
  end
else
  sliceList = eval(params.sliceList);
end

scanList = eval(params.scanList);
mrSliceExport(thisView, [params.horizontalRange params.verticalRange], sliceList, params.tifFileName, params.nRows, scanList)


% --------------------------------------------------------------------
function applyMotionCompTransformsCallBack(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');

applyMotionCompTransform(thisView);


% --------------------------------------------------------------------
function copyBaseSformCallBack(hObject, dump)

mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');
MLR.clipboard = viewGet(thisView,'baseSform');


% --------------------------------------------------------------------
function pasteBaseSformCallBack(hObject, dump)

mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');
if all(size(MLR.clipboard)==4)
    viewSet(thisView,'baseSform',MLR.clipboard);
else
    mrErrorDlg('(paste base sform) Cannot paste. Clipboard does not contain a valid transformation matrix. Use Edit -> Base Anatomy -> Transforms -> Copy sform.')
end
refreshMLRDisplay(viewNum);


