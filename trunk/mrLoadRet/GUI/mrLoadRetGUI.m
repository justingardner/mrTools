function varargout = mrLoadRetGUI(varargin)
% fig = mrLoadRetGUI('viewNum',viewNum)
% $Id$
%
% Creates a new mrLoadRet GUI.
% This function was created along with mrLoadRetGui.fig using GUIDE.
%
% djh, 6/2004
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 09-Aug-2008 14:10:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mrLoadRetGUI_OpeningFcn, ...
    'gui_OutputFcn',  @mrLoadRetGUI_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mrLoadRetGUI is made visible.
function mrLoadRetGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Parse varargin
% mrLoadRetGUI must be called as follows:
%    mrLoadRetGUI('viewNum',viewNum)
% The view structure must already exist in MRL.views{viewNum}
for index = 1:2:length(varargin)
    field = varargin{index};
    val = varargin{index+1};
    switch field
        case 'viewNum'
            viewNum = val;
            handles.viewNum = viewNum;
        otherwise
            warn('mrLoadRetGUI: invalid initialization argument')
    end
end

% Initialize coords and slice orientation
handles.coords = [1,1,1];
handles.sliceOrientation = 3;
guidata(hObject,handles);

% Initialize group popup
set(handles.groupPopup,'String',viewGet([],'groupNames'));

% Initialize rotate slider
set(handles.rotateSlider,'sliderStep',[1/360 10/360]);

% Enable/disable various widgets depending on viewType
% removed the switch here, since we always have one viewType

% Initialize the slice orientation radio buttons
set(handles.sagittalRadioButton,'Value',1);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',0);
set(handles.corticalDepth,'Visible','off');
set(handles.corticalDepthSlider,'Visible','off');
set(handles.corticalDepthSlider,'Value',0.5);
set(handles.corticalDepthText,'Visible','off');

% Choose default command line output for mrLoadRetGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
% Output function
function varargout = mrLoadRetGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
% Button Down functions

function axis_ButtonDownFcn(hObject, eventdata, handles)

function figure_ButtonDownFcn(hObject, eventdata, handles)

function figure_WindowButtonDownFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Resize
function figure_ResizeFcn(hObject, eventdata, handles)
% Change the axis size to fill the figure, leaving room at the bottom for
% the widgets.
if exist('handles','var') & ~isempty(handles)
    figureSize = get(handles.figure,'Position');
    axisSize = get(handles.axis,'Position');
    axisSize(3) = figureSize(3);
    axisSize(4) = figureSize(4) - axisSize(2);
    set(handles.axis,'Position',axisSize);
    refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
% Create function
function figure_CreateFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Delete  and Close functions

function axis_DeleteFcn(hObject, eventdata, handles)

function figure_DeleteFcn(hObject, eventdata, handles)
if ~ieNotDefined('handles') & isfield(handles,'viewNum')
    % Delete the view
    viewNum = handles.viewNum;
    deleteView(viewNum);
end

function figure_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radio buttons

% --- Sagittal
function sagittalRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
view = viewSet(viewNum,'sliceOrientation','sagittal');
refreshMLRDisplay(viewNum);

% --- Coronal
function coronalRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
view = viewSet(viewNum,'sliceOrientation','coronal');
refreshMLRDisplay(viewNum);

% --- Axial
function axialRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
view = viewSet(viewNum,'sliceOrientation','axial');
refreshMLRDisplay(viewNum);

% --- Left
function leftRadioButton_Callback(hObject, eventdata, handles)
% *** Not tested ***
viewNum = handles.viewNum;
mlrGuiSet(viewNum,'slice',1);
refreshMLRDisplay(viewNum);

% --- Right
function rightRadioButton_Callback(hObject, eventdata, handles)
% *** Not tested ***
viewNum = handles.viewNum;
mlrGuiSet(viewNum,'slice',2);
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Popups

% --- Group
function groupPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function groupPopup_Callback(hObject, eventdata, handles)
mrGlobals
viewNum = handles.viewNum;
% Update the group number
viewSet(viewNum,'curGroup',get(hObject,'Value'));
% Delete the overlays
%MLR.views{viewNum}.analyses = [];
%MLR.views{viewNum}.curAnalysis = [];
% Update nScans
mlrGuiSet(viewNum,'nScans',viewGet(viewNum,'nScans'));
% Refresh the display
refreshMLRDisplay(viewNum);

% --- ROI
function roiPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function roiPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curROI',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Analysis
function analysisPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function analysisPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curAnalysis',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Overlay
function overlayPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curOverlay',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Base image
function basePopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function basePopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curBase',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sliders and associated text fields

% --- Scan
function scanSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function scanText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function scanSlider_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'scanText',value);
view = viewSet(view,'curScan',value);
refreshMLRDisplay(viewNum);

function scanText_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(str2num(get(hObject,'String')));
mlrGuiSet(viewNum,'scan',value);
view = viewSet(view,'curScan',value);
viewSet(view,'curScan',get(handles.scanSlider,'Value'));
refreshMLRDisplay(viewNum);

% --- Slice
function sliceSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function sliceText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function sliceSlider_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'sliceText',value);
% set the current slice (get the value from the slider
% since that will be properly clipped by the mlrGuiSet call above)
viewSet(view,'curSlice',get(handles.sliceSlider,'Value'));
refreshMLRDisplay(viewNum);

function sliceText_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(str2num(get(hObject,'String')));
mlrGuiSet(viewNum,'slice',value);
viewSet(view,'curSlice',value);
refreshMLRDisplay(viewNum);

% --- Executes on slider movement.
function corticalDepthSlider_Callback(hObject, eventdata, handles)
% hObject    handle to corticalDepthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
value = get(hObject,'Value');
mlrGuiSet(viewNum,'corticalDepth',value);
refreshMLRDisplay(viewNum);


% --- Executes during object creation, after setting all properties.
function corticalDepthSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corticalDepthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function corticalDepthText_Callback(hObject, eventdata, handles)
% hObject    handle to corticalDepthText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corticalDepthText as text
%        str2double(get(hObject,'String')) returns contents of corticalDepthText as a double

viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
mlrGuiSet(viewNum,'corticalDepth',value);
refreshMLRDisplay(viewNum);

% --- Executes during object creation, after setting all properties.
function corticalDepthText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corticalDepthText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- baseGamma
function baseGammaSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseGammaText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseGammaSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(handles.viewNum,'baseGamma',value);
refreshMLRDisplay(viewNum);

function baseGammaText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'baseGamma',value);
refreshMLRDisplay(viewNum);

% --- baseTilt
function baseTiltSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseTiltText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseTiltSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'tilt',value);
v = viewGet([],'view',viewNum);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

function baseTiltText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if value < 0,value = 0;end
if value > 360,value = 360;end
viewSet(viewNum,'tilt',value);
v = viewGet([],'view',viewNum);
fig = viewGet(v,'figNum');
gui = guidata(fig);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

% --- overlayMax
function overlayMaxSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMaxText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMaxSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'overlayMax',value);
refreshMLRDisplay(viewNum);

function overlayMaxText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'overlayMax',value);
refreshMLRDisplay(viewNum);

% --- overlayMin
function overlayMinSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMinText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMinSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'overlayMin',value);
refreshMLRDisplay(viewNum);

function overlayMinText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'overlayMin',value);
refreshMLRDisplay(viewNum);

% --- alpha
function alphaSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function alphaText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function alphaSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'alpha',value);
refreshMLRDisplay(viewNum);

function alphaText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'alpha',value);
refreshMLRDisplay(viewNum);

% --- rotate
function rotateSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotateText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotateSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
mlrGuiSet(viewNum,'rotate',value);
v = viewGet([],'view',viewNum);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

function rotateText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
mlrGuiSet(viewNum,'rotate',value);
v = viewGet([],'view',viewNum);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function fileBaseMenu_Callback(hObject, eventdata, handles)

function loadAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadAnat(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function loadFromVolumeMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
% get the volume directory from prefs
volumeDirectory = mrGetPref('volumeDirectory');
saveVolumeDirectory = 0;
if isempty(volumeDirectory)
    saveVolumeDirectory = 1
end
% load the anatomy
[view volumeDirectory] = loadAnat(view,'',volumeDirectory);
refreshMLRDisplay(viewNum);
% if volumeDirectory prefMenu was empty before, than save it now
if saveVolumeDirectory
    mrSetPref('volumeDirectory',volumeDirectory);
end

% --------------------------------------------------------------------
function useCurrentScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get the tseries path/filename
tSeriesPathStr = viewGet(v,'tSeriesPathStr',viewGet(v,'curScan'));
% load that as an anatome
if isfile(tSeriesPathStr)
    v = loadAnat(MLR.views{viewNum},getLastDir(tSeriesPathStr),fileparts(tSeriesPathStr));
    refreshMLRDisplay(viewNum);
else
    mrWarnDlg(sprintf('(mrLoadRetGUI) Could not find tSeries %s',tSeriesPathStr));
end

% --------------------------------------------------------------------
function importFlatMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to importFlatMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% v = loadFlat(v);
base = importFlatOFF;
if ~isempty(base)
  viewSet(v, 'newbase', base);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function importSurfaceMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to importSurfaceMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
base = importSurfaceOFF;
if ~isempty(base)
  viewSet(v, 'newbase', base);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function saveAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentBase');
saveAnat(MLR.views{viewNum},n,1);

% --------------------------------------------------------------------
function SaveAsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentBase');
saveAnat(MLR.views{viewNum},n,1,1);

% --------------------------------------------------------------------
function fileAnalysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadAnalysis(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function saveAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentAnalysis');
if ~isempty(n)
  saveAnalysis(MLR.views{viewNum},n);
end

% --------------------------------------------------------------------
function saveAllAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberAnalyses = viewGet(viewNum,'numberofAnalyses');
for n = 1:numberAnalyses
    saveAnalysis(MLR.views{viewNum},n);
end

% --------------------------------------------------------------------
function fileOverlayMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadOverlay(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function saveOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentOverlay');
m = viewGet(viewNum,'currentAnalysis');
saveOverlay(MLR.views{viewNum},n,m);

% --------------------------------------------------------------------
function saveAllOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberOverlays = viewGet(viewNum,'numberofOverlays');
m = viewGet(viewNum,'currentAnalysis');
for n = 1:numberOverlays
    saveOverlay(MLR.views{viewNum},n,m);
end

% --------------------------------------------------------------------
function fileRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadROI(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function loadFromVolumeDirectoryROIMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to loadFromVolumeDirectoryROIMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the volume directory from prefs
volumeDirectory = mrGetPref('volumeDirectory');
% load the rois
v = loadROI(v,[],[],volumeDirectory);
% and refresh
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function importROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = importROI(view);

% --------------------------------------------------------------------
function saveROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentroi');
saveROI(MLR.views{viewNum},n,1);

% --------------------------------------------------------------------
function saveAllROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberROIs = viewGet(viewNum,'numberofrois');
for n = 1:numberROIs
    saveROI(MLR.views{viewNum},n,1);
end

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function importGroupMenuItem_Callback(hObject, eventdata, handles)
importGroupScans;

% --------------------------------------------------------------------
function importTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importTSeries not yet implemented');

% --------------------------------------------------------------------
function importOverlayMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importOverlay not yet implemented');

% --------------------------------------------------------------------
function exportMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function exportAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportAnatomy not yet implemented');

% --------------------------------------------------------------------
function exportTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportTSeries not yet implemented');

% --------------------------------------------------------------------
function exportOverlayMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify name of exported Nifti overlay file','*.hdr');
% pathstr = [] if aborted
if ~isempty(pathstr)
    mrGlobals;
    viewNum = handles.viewNum;
    mrExport2SR(viewNum, pathstr);
end

% --------------------------------------------------------------------
function exportROIMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportROI not yet implemented');

% --------------------------------------------------------------------
function exportImageMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify name of exported image file','*.tif');
% pathstr = [] if aborted
if ~isempty(pathstr)
    img = refreshMLRDisplay(handles.viewNum);
    imwrite(img, pathstr, 'tif');
end


% --------------------------------------------------------------------
function readmeMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
createReadme(MLR.session,MLR.groups);

% --------------------------------------------------------------------
function saveSessionMenuItem_Callback(hObject, eventdata, handles)
saveSession(1);

% --------------------------------------------------------------------
function printMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to printMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrPrint(v);

% --------------------------------------------------------------------
function quitMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrQuit(1,v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function editMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function editSessionMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
session = editSessionGUI('session',MLR.session);
% session is empty if aborted
if ~isempty(session)
    MLR.session = session;
    saveSession(0);
end

% --------------------------------------------------------------------
function editGroupMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
userInput = inputdlg('Enter name for new group: ','New group');
if ~isempty(userInput)
    groupName = userInput{1};
    view = viewSet(view,'newGroup',groupName);
end

% --------------------------------------------------------------------
function deleteGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
view = viewSet(view,'deleteGroup',groupNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
curGroup = viewGet(view,'currentGroup');
group = viewGet(view,'group',curGroup);
if ~isempty(group.scanParams)
    group = editGroupGUI('group',MLR.groups(curGroup));
end
% group is empty if aborted
if ~isempty(group)
    MLR.groups(curGroup) = group;
    saveSession(0);
end

% --------------------------------------------------------------------
function infoGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'curGroup');
groupName = viewGet(view,'groupName',groupNum);
disp(sprintf('\n===== Group Info (%s) =====',groupName));
groupInfo(groupNum,0,1);
disp(sprintf('======================'));

% --------------------------------------------------------------------
function editScanMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function addScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

% Choose nifti file to add
tseriesDir = viewGet(view,'tseriesDir');
pathStr = getPathStrDialog(tseriesDir,'Add scan: choose nifti file',['*',mrGetPref('niftiFileExtension')]);
if isempty(pathStr)
    % Aborted
    return
end
if ~exist(pathStr,'file')
    mrErrorDlg('Invalid nifti file');
end

% Copy the nifti file to the tseries directory
[dir,file,ext,versn] = fileparts(pathStr);
if strcmp(dir,tseriesDir)
    % If tseries file is already in the tseries directory, then use it
    hdr = cbiReadNiftiHeader(pathStr);
    fileName = [file,ext];
else
    % Copy file to the tseries directory
    fileName = ['tseries-',datestr(now,'mmddyy-HHMMSS'),mrGetPref('niftiFileExtension')];
    [data,hdr] = cbiReadNifti(pathStr);
    newPathStr = fullfile(tseriesDir,fileName);
    [bytes,hdr] = cbiWriteNifti(newPathStr,data,hdr);
end

% Add it
scanParams.fileName = fileName;
view = viewSet(view,'newScan',scanParams);

% Open GUI to set description, junk frames, nframes
curGroup = viewGet(view,'currentGroup');
group = editGroupGUI('group',MLR.groups(curGroup));
% group is empty if aborted
if ~isempty(group)
    MLR.groups(curGroup) = group;
    saveSession(0);
end

% --------------------------------------------------------------------
function copyScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'currentScan');
scan.scanParams = viewGet(view,'scanParams',scanNum);
scan.scanParams.fileName = [];
% just save the filename instead of loading the whole tSeries
% since some scans are very large
%scan.tseries = loadTSeries(view,scanNum,'all');
scan.tseries = viewGet(view,'tseriesPathStr',scanNum);
MLR.clipboard = scan;

% --------------------------------------------------------------------
function pasteScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if isfield(MLR.clipboard,'tseries') & isfield(MLR.clipboard,'scanParams') & ...
        isscan(MLR.clipboard.scanParams)
    view = saveNewTSeries(view,MLR.clipboard.tseries,MLR.clipboard.scanParams,MLR.clipboard.scanParams.niftiHdr);
else
    mrErrorDlg('(paste scan) Cannot paste. Clipboard does not contain a valid scan. Use Edit -> Scan -> Copy Scan.')
end

% --------------------------------------------------------------------
function deleteScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if viewGet(view,'nScans') == 0
  disp(sprintf('(mrLoadRetGUI) No scans in group %s to delete',viewGet(view,'groupName')));
  return
end
view = deleteScans(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function transformsMenu_Callback(hObject, eventdata, handles)
%mrWarnDlg('transforms not yet implemented');

% --------------------------------------------------------------------
function sformScanMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to sformScanMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scanSform
scanSform = viewGet(v,'scansform');
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'sform',scanSform,'The sform is usually set by mrAlign to specify the transformation from the scan coordinates to the volume anatomy coordinates. Only change this here if you know what you are doing! Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'scanSform');

% ask the user if they are really sure before actually changing it
if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
  if strcmp(answer,'Yes')
    v = viewSet(v,'scanSform',params.sform);
    v = viewSet(v,'sformCode',params.sformCode);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end

% --------------------------------------------------------------------
function scan2BaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scan2base
scan2base = inv(viewGet(v,'base2scan'));
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'scan2base',scan2base,'This tells you the transformation from the scan coordinates to the base coordinates. If you have set the sfroms properly with mrAlign this should give an easily interpretable value. For instance if you have the same slices, but voxels are twice as big in the scan, then the diagonal elements should have 0.5 in them. This can be fixed here, but should only be done if you really know what you are doing. Otherwise this should be fixed by rerunning mrAlign. Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'editable=0','This gets set to 1 if mrAlign changes the sform (3 for talairach). If it is 0 it means the sform has never been set. If you change the base2scan, it will automatically be set to 1.'}};

params = mrParamsDialog(paramsInfo,'scan2base transformation');

if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
  if strcmp(answer,'Yes')
    v = viewSet(v,'scan2base',params.scan2base);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end

% --------------------------------------------------------------------
function dicomInfoMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
s = viewGet(v,'curScan');
g = viewGet(v,'curGroup');
dicomInfo(s,g);

% --------------------------------------------------------------------
function infoScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'curScan');
groupNum = viewGet(view,'curGroup');
groupName = viewGet(view,'groupName',groupNum);
scanInfo(scanNum,groupNum,1);


% --------------------------------------------------------------------
function linkStimfileMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
% Choose stimfile to add
etcDir = viewGet(view,'etcDir');
pathStr = getPathStrDialog(etcDir,'Link Stimfile: choose matlab file');
if isempty(pathStr)
    % Aborted
    return
end
if ~exist(pathStr,'file')
    mrErrorDlg('Invalid stimfile file');
end

% Copy the nifti file to the tseries directory
[dir,file,ext,versn] = fileparts(pathStr);
if strcmp(dir,etcDir)
    % If stimfile is already in the Etc directory
    fileName = [file,ext];
    viewSet(view, 'stimfilename', fileName);
else
    mrErrorDlg('Stimfile must be in the Etc directory');
end

% --------------------------------------------------------------------
function editAnalysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
userInput = inputdlg('Enter name for new analysis: ','New analysis');
if ~isempty(userInput)
    analysis.name = userInput{1};
    analysis.type = 'dummy';
    analysis.groupName = viewGet(view,'groupName',viewGet(view,'currentGroup'));
    analysis.function = 'dummyAnalysis';
    analysis.reconcileFunction = 'dummyAnalysisReconcileParams';
    analysis.reconcileFunction = 'dummyAnalysisMergeParams';
    analysis.guiFunction = 'dummyAnalysisGUI';
    analysis.params = [];
    analysis.overlays =[];
    analysis.curOverlay = [];
    analysis.date = datestr(now);
    view = viewSet(view,'newanalysis',analysis);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function copyAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'analysis');

% --------------------------------------------------------------------
function pasteAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check analysis] = isanalysis(MLR.clipboard);
if check
    view = viewSet(view,'newAnalysis',analysis);
else
    mrErrorDlg('(paste analysis) Cannot paste. Clipboard does not contain a valid analysis. Use Edit -> Analysis -> Copy Analysis.')
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
%view = editAnalysisGUI(view);

paramsInfo{1} = {'name',viewGet(view,'analysisName'),'Name of analysis'};
paramsInfo{end+1} = {'date',viewGet(view,'analysisDate'),'editable=0','Date analysis was done (not editable)'};
paramsInfo{end+1} = {'function',viewGet(view,'analysisFunction'),'editable=0','Function with which analysis was done (not editable)'};
paramsInfo{end+1} = {'groupName',viewGet(view,'analysisGroupName'),'editable=0','Group on which analysis was done (not editable)'};
paramsInfo{end+1} = {'GUIFunction',viewGet(view,'analysisGUIFunction'),'GUI function for analysis'};
paramsInfo{end+1} = {'reconcileFunction',viewGet(view,'analysisReconcileFunction'),'Reconcile function for analysis'};
paramsInfo{end+1} = {'mergeFunction',viewGet(view,'analysisMergeFunction'),'Merge function for analysis'};
paramsInfo{end+1} = {'type',viewGet(view,'analysisType'),'Type of analysis'};

params = mrParamsDialog(paramsInfo,'Edit analysis');
if ~isempty(params)
  view = viewSet(view,'analysisName',params.name);
  view = viewSet(view,'analysisGUIFunction',params.GUIFunction);
  view = viewSet(view,'analysisReconcileFunction',params.reconcileFunction);
  view = viewSet(view,'analysisMergeFunction',params.mergeFunction);
  view = viewSet(view,'analysisType',params.type);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function EditAnalysisInfoMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to EditAnalysisInfoMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% no current anatomy, just return
if isempty(viewGet(v,'curAnalysis')),return;end

disppercent(-inf,'Gathering analysis info');
% get the current analysis
a = viewGet(v,'Analysis',viewGet(v,'curAnalysis'));

% get fields
fields = fieldnames(a);
fields = setdiff(fields,{'d','overlays','params','curOverlay'});

% make into a display
paramsInfo = {};
for fnum = 1:length(fields)
  if ~isstruct(fields)
    paramsInfo{end+1} = {fields{fnum},a.(fields{fnum}),'editable=0'};
  end
end

% check d
if isfield(a,'d')
  dExists = [];
  for dnum = 1:length(a.d)
    dExists(dnum) = ~isempty(a.d{dnum});
  end
  paramsInfo{end+1} = {sprintf('dScans'),num2str(find(dExists)),'editable=0',sprintf('Scans that d structure exists for')};
end

% check overlays
for onum = 1:length(a.overlays)
  for snum = 1:length(a.overlays(onum).data)
    overlayExists(snum) = ~isempty(a.overlays(onum).data{snum});
  end
  % make params for this
  paramsInfo{end+1} = {sprintf('overlay%i',onum),a.overlays(onum).name,'editable=0',sprintf('Name of overlay %i',onum)};
  paramsInfo{end+1} = {sprintf('overlay%iScans',onum),num2str(find(overlayExists)),'editable=0',sprintf('Scans that overlay %s exists for',a.overlays(onum).name)};
end

paramsInfo{end+1} = {'params',[],'View analysis parameters','type=pushbutton','buttonString=View analysis parameters','callback',@viewAnalysisParams,'callbackArg',v};

disppercent(inf);

% display parameters
mrParamsDialog(paramsInfo,'Analysis Info');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to view params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = viewAnalysisParams(v)

% bogus return value
retval = [];

% get analysis info
curAnalysis = viewGet(v,'curAnalysis');
params = viewGet(v,'analysisParams',curAnalysis);
guiFunction = viewGet(v,'analysisGuiFunction',curAnalysis);
groupName = viewGet(v,'analysisGroupName',curAnalysis);

% check for function
while exist(sprintf('%s.m',stripext(guiFunction))) ~= 2
  paramsInfo = {{'GUIFunction',guiFunction,sprintf('The GUI function for this analysis, %s, was not found. If you want to specify another function name you can enter that here and try again.',guiFunction)}};
  paramsGUIFunction = mrParamsDialog(paramsInfo,sprintf('GUI function: %s not found',guiFunction));
  if isempty(paramsGUIFunction) || strcmp(paramsGUIFunction.GUIFunction,guiFunction)
    return
  else
    guiFunction = paramsGUIFunction.GUIFunction;
  end
end

% params = guiFunction('groupName',groupName,'params',params);
evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params);'];

eval(evalstring);


% --------------------------------------------------------------------
function editOverlayMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'overlay');

% --------------------------------------------------------------------
function pasteOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check overlay] = isoverlay(MLR.clipboard);
if ~check
    mrErrorDlg('(paste overlay) Cannot paste. Clipboard does not contain a valid overlay. Use Edit -> Overlay -> Copy Overlay.')
end
if ~isanalysis(viewGet(view,'analysis'))
    mrErrorDlg('(paste overlay) Overlays must be pasted into an analysis. Use Edit -> Analysis -> New Analysis.')
end
view = viewSet(view,'newOverlay',overlay);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
editOverlayGUImrParams(viewNum);
% view = MLR.views{viewNum};
% view = editOverlayGUI(view);
% view = viewSet(view,'overlayCache','init');
% refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function overlayInfoMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to overlayInfoMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

overlayInfo(v);

% --------------------------------------------------------------------
function editRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'roi');

% --------------------------------------------------------------------
function pasteRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
% Check to see that it is a valid ROI structure and then add it.
[check roi] = isroi(MLR.clipboard);
if check
    view = viewSet(view,'newROI',roi);
else
    mrErrorDlg('(paste ROI) Cannot paste. Clipboard does not contain a valid ROI. Use Edit -> ROI -> Copy ROI.')
end
% Select it and reset view.prevCoords
ROInum = viewGet(view,'numberofROIs');
if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get the roi
roiNum = viewGet(v,'currentROI');
if isempty(roiNum),return,end
roiName = viewGet(v,'roiName',roiNum);
roiColor = viewGet(v,'roiColor',roiNum);
roiNotes = viewGet(v,'roiNotes',roiNum);

% get the list of colors, putting our color on top
colors = putOnTopOfList(roiColor,color2RGB);

% make parameter string
roiParams{1} = {'name',roiName,'Name of roi, avoid using punctuation and space'};
roiParams{2} = {'color',colors,'The color that the roi will display in'};
roiParams{3} = {'notes',roiNotes,'Brief notes about the ROI'};

params = mrParamsDialog(roiParams,'Edit ROI',1.5);

% if not empty, then change the parameters
if ~isempty(params)
  v = viewSet(v,'roiColor',params.color,roiNum);
  v = viewSet(v,'roiName',params.name,roiNum);
  v = viewSet(v,'roiNotes',params.notes,roiNum);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function editManyROIsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to editManyROIsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
nROIs = viewGet(v,'numROIs');

paramsInfo = {};
for roinum = 1:nROIs
  % get name and colors for each roi
  roiNames{roinum} = viewGet(v,'roiName',roinum);
  roiNotes = viewGet(v,'roiNotes',roinum);
  colors = putOnTopOfList(viewGet(v,'roiColor',roinum),color2RGB);
  paramsInfo{end+1} = {sprintf('%sName',fixBadChars(roiNames{roinum})),roiNames{roinum},'Name of roi, avoid using punctuation and space'};
  paramsInfo{end+1} = {sprintf('%sColor',fixBadChars(roiNames{roinum})),colors,sprintf('The color that roi %s will display in',roiNames{roinum})};
  paramsInfo{end+1} = {sprintf('%sNotes',fixBadChars(roiNames{roinum})),roiNotes,sprintf('Note for roi %s',roiNames{roinum})};
end
if isempty(paramsInfo),return,end
params = mrParamsDialog(paramsInfo,'Edit Many ROIs');

% if not empty, then change the parameters
if ~isempty(params)
  for roinum = 1:nROIs
    roiName = fixBadChars(roiNames{roinum});
    v = viewSet(v,'roiColor',params.(sprintf('%sColor',roiName)),roinum);
    v = viewSet(v,'roiName',params.(sprintf('%sName',roiName)),roinum);
    v = viewSet(v,'roiNotes',params.(sprintf('%sNotes',roiName)),roinum);
  end
  refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
function editAllROIsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to editAllROIsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
nROIs = viewGet(v,'numROIs');

if nROIs == 0,return,end
paramsInfo = {};
colors = putOnTopOfList(viewGet(v,'roiColor'),color2RGB);
% get color to set all ROIs to
paramsInfo{end+1} = {'changeColor',0,'type=checkbox','Click to change color of all ROIs'};
paramsInfo{end+1} = {'color',colors,'type=popupmenu','contingent=changeColor','Color for all ROIs'};
paramsInfo{end+1} = {'changeNotes',0,'type=checkbox','Click to change notes of all ROIs'};
paramsInfo{end+1} = {'notes',viewGet(v,'roiNotes'),'contingent=changeNotes','Notes for all ROIs'};
params = mrParamsDialog(paramsInfo,'Edit All ROIs',1.5);

% if not empty, then change the parameters
if ~isempty(params)
  for roinum = 1:nROIs
    if ~isempty(params.color)
      v = viewSet(v,'roiColor',params.color,roinum);
    end
    if ~isempty(params.notes)
      v = viewSet(v,'roiNotes',params.notes,roinum);
    end
  end
  refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
function infoROIMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

roiInfo(v);

% --------------------------------------------------------------------
function editBaseMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'baseAnatomy');

% --------------------------------------------------------------------
function pasteBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check base] = isbase(MLR.clipboard);
if check
    view = viewSet(view,'newBase',base);
else
    mrErrorDlg('(paste base anatomy) Cannot paste. Clipboard does not contain a valid scan. Use Edit -> Base Anatomy -> Copy Base Anatomy.')
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
v = editBaseGUI(v);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function baseTransformsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to baseTransformsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function baseSformMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to baseSformMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the baseSform
baseSform = viewGet(v,'baseSform');
sformCode = viewGet(v,'basesformCode');
% params dialog
paramsInfo = {{'sform',baseSform,'The sform is usually set by mrAlign to specify the transformation from the base coordinates to the volume anatomy coordinates. Only change this here if you really know what you are doing!'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'baseSform');

% ask the user if they are really sure before actually changing it
if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign.');
  if strcmp(answer,'Yes')
    v = viewSet(v,'baseSform',params.sform);
    v = viewSet(v,'baseSformCode',params.sformCode);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end


% --------------------------------------------------------------------
function base2scanMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to base2scanMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the base2scan xform
base2scan = viewGet(v,'base2scan');
baseSformCode = viewGet(v,'baseSformCode');
% params dialog
paramsInfo = {{'base2scan',base2scan,'This tells you the transformation from the base coordinates to the scan coordinates. If you have set the sfroms properly with mrAlign this should give an easily interpretable value. For instance if you have the same slices, but voxels are twice as big in the scan, then the diagonal elements should have 0.5 in them. This can be fixed here, but should only be done if you really know what you are doing. Otherwise this should be fixed by rerunning mrAlign.'},...
    {'baseSformCode',baseSformCode,'editable=0','This gets set to 1 if mrAlign changes the sform (3 for talairach). If it is 0 it means the sform has never been set. If you change the base2scan, it will automatically be set to 1.'}};

params = mrParamsDialog(paramsInfo,'base2scan transformation');

if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign.');
  if strcmp(answer,'Yes')
    v = viewSet(v,'base2scan',params.base2scan);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end


% --------------------------------------------------------------------
function infoBaseAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
baseInfo(view);

% --------------------------------------------------------------------
function prefMenu_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get settings from start
interpMethod = mrGetPref('interpMethod');
selectedROIColor = mrGetPref('selectedROIColor');

% remember old cache sizes
roiCacheSize = mrGetPref('roiCacheSize');
baseCacheSize = mrGetPref('baseCacheSize');
overlayCacheSize = mrGetPref('overlayCacheSize');

% prefs dialog
prefParams = mrEditPrefs;

% if changes, check for change in cache size
if ~isempty(prefParams)
  if (roiCacheSize ~= prefParams.roiCacheSize)
    v = viewSet(v,'roiCache','init');
  end
  if (baseCacheSize ~= prefParams.baseCacheSize)
    v = viewSet(v,'baseCache','init');
  end
  if (overlayCacheSize ~= prefParams.overlayCacheSize)
    v = viewSet(v,'overlayCache','init');
  end

  % check to see if interpolation method has changed
  if ~strcmp(interpMethod,prefParams.interpMethod)
    % dump overlay cache and redraw
    v = viewSet(v,'overlayCache','init');
    refreshMLRDisplay(v.viewNum);
  end

  % check to see if interpolation method has changed
  if ~strcmp(selectedROIColor,prefParams.selectedROIColor)
    refreshMLRDisplay(v.viewNum);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function windowMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function graphMenuItem_Callback(hObject, eventdata, handles)
newGraphWin;

% --------------------------------------------------------------------
function newWindowMenuItem_Callback(hObject, eventdata, handles)
view = mrOpenWindow('',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function motionCompMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function motionCompMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionComp(view);

% --------------------------------------------------------------------
function motionCompWithinMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionCompWithinScan(view);

% --------------------------------------------------------------------
function motionCompBetweenMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionCompBetweenScans(view);

% --------------------------------------------------------------------
function sliceTimeCorrectMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = sliceTimeCorrect(view);

% --------------------------------------------------------------------
function averageTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = averageTSeries(view);

% --------------------------------------------------------------------
function concatenateTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = concatTSeries(view);

% --------------------------------------------------------------------
function tsStatsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = timeSeriesStats(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function corAnalMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = corAnal(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function eventRelatedMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = eventRelated(view);

% --------------------------------------------------------------------
function glmMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = eventRelatedGlm(view);

% --------------------------------------------------------------------
function recomputeAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if viewGet(view,'numAnalyses') > 0
  n = viewGet(view,'currentAnalysis');
  groupName = viewGet(view,'analysisGroupName',n);
  analysisFunction = viewGet(view,'analysisFunction',n);
  guiFunction = viewGet(view,'analysisGuiFunction',n);
  params = viewGet(view,'analysisParams',n);
  % params = guiFunction('groupName',groupName,'params',params);
  evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params);'];
  eval(evalstring);
  % params is empty if GUI cancelled
  if ~isempty(params)
    % view = analysisFunction(view,params);
    evalstring = ['view = ',analysisFunction,'(view,params);'];
    eval(evalstring);
    refreshMLRDisplay(viewNum);
  end
else
  mrWarnDlg(sprintf('(mrLoadRetGUI) No analyses loaded'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function deleteBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
baseNum = viewGet(view,'currentBase');
if ~isempty(baseNum)
  view = viewSet(view,'deleteBase',baseNum);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllBasesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numBases = viewGet(view,'numberofbasevolumes');
for baseNum = numBases:-1:1;
    view = viewSet(view,'deleteBase',baseNum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
analysisNum = viewGet(view,'currentAnalysis');
view = viewSet(view,'deleteAnalysis',analysisNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteAllAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numAnalyses = viewGet(view,'numberofAnalyses');
for analysisNum = numAnalyses:-1:1;
    view = viewSet(view,'deleteAnalysis',analysisNum);
end
% also, go through and delete all "laodedAnalyses", i.e.
% analysis that are loaded in ohter groups
for g = 1:viewGet(view,'numGroups');
  view = viewSet(view,'loadedAnalyses',[],g);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
view = viewSet(view,'deleteOverlay',overlayNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteAllOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numOverlays = viewGet(view,'numberofoverlays');
for overlayNum = numOverlays:-1:1;
    view = viewSet(view,'deleteOverlay',overlayNum);
end
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function createRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'rectangle',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createPolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'polygon',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'line',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function deleteRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roinum = viewGet(view,'currentROI');
view = viewSet(view,'deleteROI',roinum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeManyROIMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to removeManyROIMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
numrois = viewGet(v,'numberofrois');
% put up a dialog with rois to delete
roinames = viewGet(v,'roiNames');
paramsDialog = {};
for roinum = 1:length(roinames)
  paramsDialog{end+1} = {fixBadChars(roinames{roinum}),0,'type=checkbox',sprintf('Remove ROI %i: %s',roinum,roinames{roinum})};
end
params = mrParamsDialog(paramsDialog,'Select ROIs to remove');
if ~isempty(params)
  % now go through and delete anything the user selected
  for roinum = 1:length(roinames)
    if params.(fixBadChars(roinames{roinum}))
      v = viewSet(v,'deleteROI',viewGet(v,'roinum',roinames{roinum}));
    end
  end
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numrois = viewGet(view,'numberofrois');
for roinum = numrois:-1:1;
    view = viewSet(view,'deleteROI',roinum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function addRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'rectangle',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addPolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'polygon',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'line',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function removeRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'rectangle',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removePolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'polygon',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'line',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function combineROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roiNames = viewGet(view,'roiNames');
paramInfo = {...
  {'combineROI',putOnTopOfList(viewGet(view,'roiName'),viewGet(view,'roiNames')),'editable=0','The ROI that will be combined with the otherROI'},...
  {'otherROI',roiNames,'The otherROI is combined with the combineROI and the result is put into combineROI.'},...
  {'action',{'A not B', 'Intersection', 'Union', 'XOR'},'Select action for combining ROIs.'},...
  {'combine',0,'type=pushbutton','callback',@doCombine','passParams=1','callbackArg',viewNum,'buttonString=Do combination','Click this button to do the combination. This is the same as hitting OK but won''t close the dialog so you can continue to do more combinations'}};
params = mrParamsDialog(paramInfo,'Combine ROIs');
if ~isempty(params)
  doCombine(viewNum,params);
end

function retval = doCombine(viewNum,params)

% get the view (it is important to get it from the viewNum
% so that we always have an up-to-date view
v = viewGet([],'view',viewNum);

retval = 1;
disp(sprintf('(doCombine) %s %s %s',params.combineROI,params.action,params.otherROI));
v = combineROIs(v,params.combineROI,params.otherROI,params.action);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function restrictRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function restrictRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roinum = viewGet(view,'currentROI');
scan = viewGet(view,'curscan');
view = restrictROI(view,roinum,scan);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function restrictAllROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scan = viewGet(view,'curscan');
nROIs = viewGet(view,'numberofROIs');
for roinum = 1:nROIs
    view = restrictROI(view,roinum,scan);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function undoRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
curCoords = viewGet(view,'roiCoords');
prevCoords = viewGet(view,'prevROIcoords');
view = viewSet(view,'prevROIcoords',curCoords);
view = viewSet(view,'ROIcoords',prevCoords);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function convertCorticalDepthRoiMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to convertCorticalDepthRoiMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = convertROICorticalDepth(v);


% --------------------------------------------------------------------
function convertRoiToBaseAnatomyMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to convertRoiToBaseAnatomyMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = convertROI(v);


% --------------------------------------------------------------------
function showRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function showAllMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','all');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showAllPerimeterMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','all perimeter');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showSelectedMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','selected');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showSelectedPerimeterMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','selected perimeter');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function hideROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','hide');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function labelsROIsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to labelsROIsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% switch labels on/off
labelROIs = viewGet(v,'labelROIs');
viewSet(v,'labelROIs',~labelROIs');

% redisplay
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function findCurrentROIMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to findCurrentROIMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = findROI(v);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function interrogateOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
% start or stop the interrogator
if strcmp(get(hObject,'Checked'),'on')
    mrInterrogator('end',viewNum);
    set(hObject,'Checked','off');
else
    % start the mrInterrogator
    mrInterrogator('init',viewNum);
    set(hObject,'Checked','on');
end

return

% --------------------------------------------------------------------
function plotMeanTseriesMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function plotMeanTseriesCurrentCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesCurrentSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = selectScans(view);
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesCurrentAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = [1:viewGet(view,'nscans')];
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesAllCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = viewGet(view,'currentScan');
plotMeanFourierAmp(view, groupNum, roiList, scanList,  'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = selectScans(view);
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = selectScans(view);
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMotionCorrectionMatrices_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
d.groupNum = viewGet(view,'currentGroup');
d.scanNum = viewGet(view,'curScan');
d.expname = '';
d.ver = 4;
dispmotioncorrect(d);


% --------------------------------------------------------------------
function plotsDisplayEPIImages_Callback(hObject, eventdata, handles)
% hObject    handle to plotsDisplayEPIImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
mlrDisplayEPI(view);


% --------------------------------------------------------------------
function plotSpikeDetection_Callback(hObject, eventdata, handles)
% hObject    handle to plotSpikeDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
mlrSpikeDetector(view,viewGet(view,'curScan'),viewGet(view,'curGroup'));


% --------------------------------------------------------------------
function flatViewerMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to flatViewerMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get the flat files info which
% is stored in the baseCoordMap
params = viewGet(v,'baseCoordMap');
flatPath = viewGet(v,'baseCoordMapPath');
% switch directories to the flatDir, asking
% the user to find it if it does not exist
thispwd = pwd;
if isdir(flatPath)
  cd(flatPath);
else
  mrWarnDlg(sprintf('Directory %s does not exist, please find the anatomy folder',flatPath));
  pathStr = uigetdir(mrGetPref('volumeDirectory'),'Find anatomy folder from which this flat was created');
  if pathStr == 0,return,end
  cd(pathStr);
  viewSet(v,'baseCoordMapPath',pathStr);
end

baseType = viewGet(v,'baseType');
if baseType == 1
  % now bring up the flat viewer
  mrFlatViewer(params.flatFileName,params.outerCoordsFileName,params.innerCoordsFileName,params.curvFileName,params.anatFileName,viewNum);
else
  mrSurfViewer(params.outerSurfaceFileName,params.outerCoordsFileName,params.innerSurfaceFileName,params.innerCoordsFileName,params.curvFileName,params.anatFileName);
end

  cd(thispwd);

% --------------------------------------------------------------------
function calcDistMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function calcDistMenuItemSeg_Callback(hObject, eventdata, handles)
% hObject    handle to flatViewerMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% [view userCancel] = newROI(view);
% if userCancel,return,end
if ~exist('dijkstra') == 3
  mrWarnDlg('(mrLoadRetGUI) You need to install the mrGray tools to use the function calcDist');
else
  dist = calcDist(v, 'segments'); % ,'polygon',1);
  disp(sprintf('Distance been pts A and B: %2.4f', dist(2)));
  refreshMLRDisplay(viewNum);
end


