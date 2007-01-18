function varargout = editOverlayGUI(varargin)
% view = editOverlayGUI(view);
%
% GUI for editing overlays.
% Created by GUIDE.
%
% djh 7/2003

% Last Modified by GUIDE v2.5 03-Aug-2006 14:31:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editOverlayGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editOverlayGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before editOverlayGUI is made visible.
function editOverlayGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Parse varargin
% Must be called as follows:
%    editOverlayGUI(view)
if length(varargin) < 1
	mrErrorDlg('editOverlayGUI: invalid initialization argument');
end
view = varargin{1};
if ~isview(view)
	mrErrorDlg('editOverlayGUI: invalid initialization argument');
end

% Current overlay
analysisNum = viewGet(view,'currentAnalysis');
overlayNum = viewGet(view,'currentOverlay',analysisNum);
overlayName = viewGet(view,'overlayName',overlayNum,analysisNum);
overlayCmap = viewGet(view,'overlayCmap',overlayNum,analysisNum);
overlayRange = viewGet(view,'overlayRange',overlayNum,analysisNum);

% Initialize handles
handles.view = view;
handles.overlayNum = overlayNum;
set(handles.nameEditText,'String',overlayName);
set(handles.minEditText,'String',num2str(overlayRange(1)));
set(handles.maxEditText,'String',num2str(overlayRange(2)));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editOverlayGUI wait for user response (see UIRESUME)
uiwait(handles.editOverlayDialog);


% --- Outputs from editOverlayGUI
function varargout = editOverlayGUI_OutputFcn(hObject, eventdata, handles) 
% Returns view
if ~isempty(handles)
    varargout{1} = handles.view;
    close(handles.editOverlayDialog);
else
    varargout{1} = [];
end


% --- colormapPopup
function colormapPopup_Callback(hObject, eventdata, handles)

function colormapPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- nameEditText
function nameEditText_Callback(hObject, eventdata, handles)

function nameEditText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- minEditText
function minEditText_Callback(hObject, eventdata, handles)

function minEditText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- maxEditText
function maxEditText_Callback(hObject, eventdata, handles)

function maxEditText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- okPushButton.
function okPushButton_Callback(hObject, eventdata, handles)
view = handles.view;
string = get(handles.colormapPopup,'String');
value = get(handles.colormapPopup,'Value');
cmapName = string{value};
view = viewSet(view,'overlayCmap',cmapName,handles.overlayNum);
nameString = get(handles.nameEditText,'String');
view = viewSet(view,'overlayName',nameString,handles.overlayNum);
minString = get(handles.minEditText,'String');
maxString = get(handles.maxEditText,'String');
view = viewSet(view,'overlayRange',[str2num(minString),str2num(maxString)],handles.overlayNum);
handles.view = view;
guidata(hObject, handles);
uiresume;


% --- cancelPushButton
function cancelPushButton_Callback(hObject, eventdata, handles)
uiresume;

