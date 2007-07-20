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

% Last Modified by GUIDE v2.5 20-Jul-2007 09:52:22

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
viewType = viewGet(viewNum,'viewType');
switch viewType
    case 'Volume'
        % Make left & right radio buttons invisible
        set(handles.leftRadioButton,'Visible','off');
        set(handles.rightRadioButton,'Visible','off');
        % Initialize the slice orientation radio buttons
        set(handles.sagittalRadioButton,'Value',1);
        set(handles.coronalRadioButton,'Value',0);
        set(handles.axialRadioButton,'Value',0);
    case {'Surface'}
        % Make slice orientation radio buttons invisible
        set(handles.sagittalRadioButton,'Visible','off');
        set(handles.coronalRadioButton,'Visible','off');
        set(handles.axialRadioButton,'Visible','off');
        % Make left & right radio buttons invisible
        set(handles.leftRadioButton,'Visible','off');
        set(handles.rightRadioButton,'Visible','off');
        % Make slice slider invisible
        set(handles.slice,'Visible','off');
        set(handles.sliceText,'Visible','off');
        set(handles.sliceSlider,'Visible','off');
    case 'Flat'
        % Make slice orientation radio buttons invisible
        set(handles.sagittalRadioButton,'Visible','off');
        set(handles.coronalRadioButton,'Visible','off');
        set(handles.axialRadioButton,'Visible','off');
        % Make slice slider invisible
        set(handles.slice,'Visible','off');
        set(handles.sliceText,'Visible','off');
        set(handles.sliceSlider,'Visible','off');
        % Initialize left & right radio buttons
        set(handles.leftRadioButton,'Value',1);
        set(handles.rightRadioButton,'Value',0);
end

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
viewNum = handles.viewNum;
value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'scanText',value);
refreshMLRDisplay(viewNum);

function scanText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
mlrGuiSet(viewNum,'scan',value);
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
viewNum = handles.viewNum;
value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'sliceText',value);
refreshMLRDisplay(viewNum);

function sliceText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
mlrGuiSet(viewNum,'slice',value);
refreshMLRDisplay(viewNum);

% --- baseMin
function baseMinSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseMinText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseMinSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(handles.viewNum,'baseMin',value);
refreshMLRDisplay(viewNum);

function baseMinText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'baseMin',value);
refreshMLRDisplay(viewNum);

% --- baseMax
function baseMaxSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseMaxText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseMaxSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'baseMax',value);
refreshMLRDisplay(viewNum);

function baseMaxText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
viewSet(viewNum,'baseMax',value);
refreshMLRDisplay(viewNum);

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
refreshMLRDisplay(viewNum);

function rotateText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
mlrGuiSet(viewNum,'rotate',value);
refreshMLRDisplay(viewNum);

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
function loadFlatMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to loadFlatMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
v = loadFlat(v);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function saveAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentBase');
saveAnat(MLR.views{viewNum},n,1);

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
saveAnalysis(MLR.views{viewNum},n);

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
function importAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importAnatomy not yet implemented');

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
function importROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = importROI(view);

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
function quitMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
if isfield(MLR,'views') && ~isempty(MLR.views)
    viewNum = handles.viewNum;
    thisView = MLR.views{viewNum};
    % remember figure location
    mrSetFigLoc('mrLoadRetGUI',get(thisView.figure,'Position'));
    % remember GUI settings
    viewSettings.rotate = viewGet(thisView,'rotate');
    viewSettings.curScan = viewGet(thisView,'curScan');
    viewSettings.curSlice = viewGet(thisView,'curSlice');
    viewSettings.curGroup = viewGet(thisView,'curGroup');
    viewSettings.sliceOrientation = viewGet(thisView,'sliceOrientation');
    viewSettings.overlayMin = viewGet(thisView,'overlayMin');
    viewSettings.overlayMax = viewGet(thisView,'overlayMax');
    viewSettings.alpha = viewGet(thisView,'alpha');
    % close view figures
    for viewNum = 1:length(MLR.views)
        view = MLR.views{viewNum};
        if isview(view)
            delete(view.figure);
        end
    end
    % close graph figure, remembering figure location
    if ~isempty(MLR.graphFigure)
        mrSetFigLoc('graphFigure',get(MLR.graphFigure,'Position'));
        close(MLR.graphFigure);
        MLR.graphFigure = [];
    end
    % save the view in the current directory
    view = thisView;
    eval(sprintf('save %s view viewSettings -V6;',fullfile(MLR.homeDir,'mrLastView')));
    % save .mrDefaults in the home directory
    saveMrDefaults;
else
    closereq;
end
clear global MLR

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
groupInfo(groupNum);
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
scan.tseries = loadTSeries(view,scanNum,'all');
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
scanList = selectScans(view);
for iScan = 1:length(scanList)
    % first get time series name for each one of these scans
    % since as we delete them, then numbers stop making sense
    tSeriesFile{iScan} = viewGet(view,'tSeriesFile',scanList(iScan));
end
% now go through and delete
for iScan = 1:length(scanList)
    % get the scan number
    scanNum = viewGet(view,'scanNum',tSeriesFile{iScan});
    if ~isempty(scanNum)
        scanNum = scanNum(1);
        view = viewSet(view,'deleteScan',scanNum);
        disp(sprintf('Scan for file %s deleted.',tSeriesFile{iScan}));
    else
        disp(sprintf('(mrLoadRetGUI) Could not delete scan for file %s',tSeriesFile{iScan}));
    end
    refreshMLRDisplay(viewNum);
end
if ~isempty(scanList)
    disp(sprintf('To remove the nifti files for these deleted scans run mrCleanDir'));
end

% --------------------------------------------------------------------
function transformsMenu_Callback(hObject, eventdata, handles)
mrWarnDlg('transforms not yet implemented');

% --------------------------------------------------------------------
function sformScanMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to sformScanMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scanxform
scanXform = viewGet(v,'scanxform');
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'sform',scanXform,'The sform is usually set by mrAlign to specify the transformation from the scan coordinates to the volume anatomy coordinates. Only change this here if you know what you are doing! Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'scanXform');

% ask the user if they are really sure before actually changing it
if ~isempty(params)
    answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign/mrUpdateNifitHdr. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
    if strcmp(answer,'Yes')
        v = viewSet(v,'scanXform',params.sform);
        v = viewSet(v,'sformCode',params.sformCode);
        saveSession;
        refreshMLRDisplay(viewNum);
    end
end

% --------------------------------------------------------------------
function scan2BaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scanxform
scanXform = viewGet(v,'scanXform');
baseXform = viewGet(v,'baseXform');
shiftXform = shiftOriginXform;
scan2base = inv(shiftXform)*inv(scanXform)*baseXform*shiftXform;
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'scan2base',scan2base,'This tells you the transformation from the scan coordinates to the base coordinates. If you have set the sfroms properly with mrAlign this should give an easily interpretable value. For instance if you have the same slices, but voxels are twice as big in the scan, then the diagonal elements should have 0.5 in them. This can be fixed here, but should only be done if you really know what you are doing. Otherwise this should be fixed by rerunning mrAlign and then saving out the proper transform to this scan file and then running mrUpdateNiftiHdr. Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'scan2base transformation');

if ~isempty(params)
    answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign/mrUpdateNifitHdr. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
    if strcmp(answer,'Yes')
        v = viewSet(v,'scanXform',inv(shiftXform*params.scan2base*inv(shiftXform)*inv(baseXform)));
        v = viewSet(v,'sformCode',params.sformCode);
        saveSession;
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
view = editAnalysisGUI(view);
refreshMLRDisplay(viewNum);

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
view = MLR.views{viewNum};
view = editOverlayGUI(view);
refreshMLRDisplay(viewNum);

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
view = MLR.views{viewNum};
% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end
roiName = viewGet(view,'roiName',roiNum);
roiColor = viewGet(view,'roiColor',roiNum);
colors = {'red','green','blue','yellow','cyan','magenta','white','black'};
% remove our color from list
colors = setdiff(colors,roiColor);
% and add it at the top
colors{end+1} = roiColor;
colors = fliplr(colors);
% make parameter string
roiParams = {{'name',roiName,'Name of roi, avoid using punctuation and space'},{'color',colors,'The color that the roi will display in'}};
params = mrParamsDialog(roiParams);
% if not empty, then change the parameters
if ~isempty(params)
    view = viewSet(view,'roiColor',params.color,roiNum);
    view = viewSet(view,'roiName',params.name,roiNum);
    refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function infoROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roiNum = viewGet(view,'currentROI');
roiName = viewGet(view,'roiName',roiNum);
roiDate = viewGet(view,'roidate',roiNum);
roiColor = viewGet(view,'roicolor',roiNum);
roiVoxelSize = viewGet(view,'roivoxelsize',roiNum);
roiVolume = viewGet(view,'roivolume',roiNum);
roiXform = viewGet(view,'roixform',roiNum);

paramsInfo = {{'name',roiName,'editable=0','The name of the ROI'},...
  {'date',roiDate,'editable=0','The date of creation'},...
  {'color',roiColor,'editable=0','ROI color'},...
  {'voxelsize',roiVoxelSize,'editable=0','Voxel dimensions in mm'},...
  {'volume',roiVolume,'editable=0','Volume of ROI in cubic mm'},...
  {'xform',roiXform,'editable=0','xform matrix specifies the transformation to the base coordinate system'}};
mrParamsDialog(paramsInfo,'ROI information');

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
view = MLR.views{viewNum};
view = editBaseGUI(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function infoBaseAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'curScan');
groupNum = viewGet(view,'curGroup');
baseDims = viewGet(view,'baseDims');
baseQform = viewGet(view,'baseqform');
baseSform = viewGet(view,'baseXform');
baseVolPermutation = viewGet(view,'baseVolPermutation');
baseVoxelSize = viewGet(view,'baseVoxelSize');
baseName = viewGet(view,'baseName');

paramsInfo = {{'baseName',baseName,'editable=0','The name of the base anatomy'},...
    {'voxelSize',baseVoxelSize,'editable=0','Voxel dimensions in mm'},...
    {'baseDims',baseDims,'editable=0','Dimensions of base anatomy'},...
    {'qform',baseQform,'editable=0','Qform matrix specifies the transformation to the scanner coordinate frame'},...
    {'sform',baseSform,'editable=0','Sform matrix is set by mrAlign and usually specifies the transformation to base coordinate system'}};
mrParamsDialog(paramsInfo,'Base anatomy information');

% --------------------------------------------------------------------
function prefMenu_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function windowMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function graphMenuItem_Callback(hObject, eventdata, handles)
newGraphWin;

% --------------------------------------------------------------------
function volumeMenuItem_Callback(hObject, eventdata, handles)
view = mrOpenWindow('Volume');

% --------------------------------------------------------------------
function surfaceMenuItem_Callback(hObject, eventdata, handles)
view = mrOpenWindow('Surface');

% --------------------------------------------------------------------
function flatMenuItem_Callback(hObject, eventdata, handles)
view = mrOpenWindow('Flat');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function deleteBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
baseNum = viewGet(view,'currentBase');
view = viewSet(view,'deleteBase',baseNum);
refreshMLRDisplay(viewNum);

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
view = newROI(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = newROI(view);
view = drawROI(view,'rectangle',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createPolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = newROI(view);
view = drawROI(view,'polygon',1);
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
function roiCoordinatesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

scanNum = viewGet(view,'curScan');

% get the coordinates
coords = getROICoordinates(view,roiNum,scanNum);

disp(sprintf('ROI %s: n=%i',viewGet(view,'roiName',roiNum),size(coords,2)));
% and display them to the buffer
for i = 1:size(coords,2)
    disp(sprintf('%i %i %i',coords(1,i),coords(2,i),coords(3,i)));
end

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
plotMeanFourierAmp(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = selectScans(view);
plotMeanFourierAmp(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpAllCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanFourierAmp(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = selectScans(view);
plotMeanFourierAmp(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList);

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



