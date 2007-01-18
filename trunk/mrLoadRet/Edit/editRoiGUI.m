function varargout = editRoiGUI(varargin)
% view = editRoiGUI(view);
%
% GUI for editing overlays.
% Created by GUIDE.
%
% djh 7/2003

% Last Modified by GUIDE v2.5 10-Jul-2006 21:56:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editRoiGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editRoiGUI_OutputFcn, ...
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


% --- Executes just before editRoiGUI is made visible.
function editRoiGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Parse varargin
% Must be called as follows:
%    editRoiGUI(view)
if length(varargin) < 1
	mrErrorDlg('editRoiGUI: invalid initialization argument');
end
view = varargin{1};
if ~isview(view)
	mrErrorDlg('editRoiGUI: invalid initialization argument');
end

% Current ROI
roiNum = viewGet(view,'currentROI');
roiName = viewGet(view,'roiName',roiNum);
roiColor = viewGet(view,'roiColor',roiNum);

% Initialize name text
set(handles.nameEditText,'String',roiName);

% Initialize color popup
colorPopupValue = find(strcmp(roiColor, get(handles.colorPopup,'String')));
set(handles.colorPopup,'Value',colorPopupValue);

% Choose default command line output for editRoiGUI
handles.view = view;
handles.roiNum = roiNum;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editRoiGUI wait for user response (see UIRESUME)
uiwait(handles.editRoiDialog);


% --- Outputs from editRoiGUI
function varargout = editRoiGUI_OutputFcn(hObject, eventdata, handles) 
% Returns view
if ~isempty(handles)
    varargout{1} = handles.view;
    close(handles.editRoiDialog);
else
    varargout{1} = [];
end


% --- colorPopup
function colorPopup_Callback(hObject, eventdata, handles)

function colorPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- nameEditText
function nameEditText_Callback(hObject, eventdata, handles)

function nameEditText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- okPushButton.
function okPushButton_Callback(hObject, eventdata, handles)
view = handles.view;
string = get(handles.colorPopup,'String');
value = get(handles.colorPopup,'Value');
roiColor = string{value};
view = viewSet(view,'roiColor',roiColor,handles.roiNum);
nameString = get(handles.nameEditText,'String');
view = viewSet(view,'roiName',nameString,handles.roiNum);
handles.view = view;
guidata(hObject, handles);
uiresume;

