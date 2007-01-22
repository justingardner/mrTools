function varargout = mrLoadRetGUI(varargin)
% fig = mrLoadRetGUI('viewNum',viewNum)
%
% Creates a new mrLoadRet GUI.
% This function was created along with mrLoadRetGui.fig using GUIDE.
%
% djh, 6/2004
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 13-Nov-2006 16:56:20

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
MLR.views{viewNum}.analyses = [];
MLR.views{viewNum}.curAnalysis = [];
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
saveAnalysis(MLR.views{viewNum},n,1);

% --------------------------------------------------------------------
function saveAllAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberAnalyses = viewGet(viewNum,'numberofAnalyses');
for n = 1:numberAnalyses
	saveAnalysis(MLR.views{viewNum},n,1);
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
saveOverlay(MLR.views{viewNum},n,m,1);

% --------------------------------------------------------------------
function saveAllOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberOverlays = viewGet(viewNum,'numberofOverlays');
m = viewGet(viewNum,'currentAnalysis');
for n = 1:numberOverlays
	saveOverlay(MLR.views{viewNum},n,m,1);
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
numberROIs = 1:viewGet(viewNum,'numberofrois');
for n = 1:numberROIs
	saveROI(MLR.views{viewNum},n,1);
end

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function importAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importAnatomy not yet implemented');

% --------------------------------------------------------------------
function importTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importTSeries not yet implemented');

% --------------------------------------------------------------------
function importOverlayMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importOverlay not yet implemented');

% --------------------------------------------------------------------
function importROIMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('importROI not yet implemented');

% --------------------------------------------------------------------
function exportMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function exportAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportAnatomy not yet implemented');

% --------------------------------------------------------------------
function exportOverlayMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportOverlay not yet implemented');

% --------------------------------------------------------------------
function exportTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportTSeries not yet implemented');

% --------------------------------------------------------------------
function exportROIMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportROI not yet implemented');

% --------------------------------------------------------------------
function exportImageMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportImage not yet implemented');

% --------------------------------------------------------------------
function readmeMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
createReadme(MLR.session,MLR.groups);

% --------------------------------------------------------------------
function saveSessionMenuItem_Callback(hObject, eventdata, handles)
saveSession;

% --------------------------------------------------------------------
function quitMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
for viewNum = 1:length(MLR.views)
	view = MLR.views{viewNum};
	if isview(view)
		delete(view.figure);
	end
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
function editScanMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function addScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

% Choose nifti file to add
tseriesDir = viewGet(view,'tseriesDir');
pathStr = getPathStrDialog(tseriesDir,'Add scan: choose nifti file',['*',niftiFileExtension]);
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
    fileName = ['tseries-',datestr(now,'mmddyy-HHMMSS'),niftiFileExtension];
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
	mrErrorDlg('Cannot paste: Invalid scan.');
end

% --------------------------------------------------------------------
function deleteScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'currentScan');
view = viewSet(view,'deleteScan',scanNum);
refreshMLRDisplay(viewNum);

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
if isanalysis(MLR.clipboard)
	view = viewSet(view,'newAnalysis',MLR.clipboard);
else
	mrErrorDlg('Cannot paste: Invalid analysis.');
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
if isoverlay(MLR.clipboard)
	view = viewSet(view,'newOverlay',MLR.clipboard);
else
	mrErrorDlg('Cannot paste: Invalid overlay.');
end
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
if isroi(MLR.clipboard)
	view = viewSet(view,'newROI',MLR.clipboard);
else
	mrErrorDlg('Cannot paste: Invalid ROI.');
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
view = editRoiGUI(view);
refreshMLRDisplay(viewNum);

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
if isbase(MLR.clipboard)
	view = viewSet(view,'newBase',MLR.clipboard);
else
	mrErrorDlg('Cannot paste: Invalid base anatomy.');
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = editBaseGUI(view);
refreshMLRDisplay(viewNum);

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
function averageTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = averageTSeries(view);

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
function corAnalMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = corAnal(view);
refreshMLRDisplay(viewNum);

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
function removeRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function removeRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'rectangle',0);
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
mrWarnDlg('restrictAllROIs not yet implemented');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function interrogateOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
interrogator = viewGet(view,'interrogator',overlayNum,analysisNum);
if isempty(interrogator)
    mrErrorDlg('Interrogate Overlay: invalide interrogator function');
end
scan = viewGet(view,'currentScan');
slice = viewGet(view,'currentSlice');
overlayCoords = viewGet(view,'cursliceoverlaycoords');
fig = viewGet(view,'figNum');
gui = guidata(fig);

h = mrMsgBox('Left click to select voxel. Right click to quit.',1);
button = 1;
while (button ~=3)
	% Select main axes of view figure for user input
	set(0,'CurrentFigure',fig);
	set(fig,'CurrentAxes',gui.axis);
	[i,j,button] = ginput(1);
	j = round(j);
	i = round(i);
	x = round(overlayCoords(j,i,1));
	y = round(overlayCoords(j,i,2));
	z = round(overlayCoords(j,i,3));
	% Draw graph
	feval(interrogator,view,overlayNum,scan,x,y,z);
end
mrCloseDlg(h);

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
% hObject    handle to plotMeanFourierAmpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
function plotMeanFourierAmpAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
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
function plotMeanFourierAmpAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = selectScans(view);
plotMeanFourierAmp(view, groupNum, roiList, scanList);


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
function plotMeanTseriesAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);


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
function eventRelatedMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = eventRelated(view);


