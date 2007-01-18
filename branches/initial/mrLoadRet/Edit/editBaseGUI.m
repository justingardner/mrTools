function varargout = editBaseGUI(varargin)
% view = editBaseGUI(view);
%
% GUI for editing base volumes.
% Created by GUIDE.
%
% djh 7/2003

% Last Modified by GUIDE v2.5 03-Aug-2006 18:52:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editBaseGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editBaseGUI_OutputFcn, ...
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


% --- Executes just before editBaseGUI is made visible.
function editBaseGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Parse varargin
% Must be called as follows:
%    editBaseGUI(view)
if length(varargin) < 1
	mrErrorDlg('editBaseGUI: invalid initialization argument');
end
view = varargin{1};
if ~isview(view)
	mrErrorDlg('editBaseGUI: invalid initialization argument');
end

% Current base
baseNum = viewGet(view,'curBase');
baseName = viewGet(view,'baseName',baseNum);
baseRange = viewGet(view,'baseRange',baseNum);

% Initialize handles
handles.view = view;
handles.baseNum = baseNum;
set(handles.nameEditText,'String',baseName);
set(handles.minEditText,'String',num2str(baseRange(1)));
set(handles.maxEditText,'String',num2str(baseRange(2)));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editBaseGUI wait for user response (see UIRESUME)
uiwait(handles.editBaseDialog);


% --- Outputs from editBaseGUI
function varargout = editBaseGUI_OutputFcn(hObject, eventdata, handles) 
% Returns view
if ~isempty(handles)
    varargout{1} = handles.view;
    close(handles.editBaseDialog);
else
    varargout{1} = [];
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
nameString = get(handles.nameEditText,'String');
view = viewSet(view,'baseName',nameString,handles.baseNum);
minString = get(handles.minEditText,'String');
maxString = get(handles.maxEditText,'String');
view = viewSet(view,'baseRange',[str2num(minString),str2num(maxString)],handles.baseNum);
handles.view = view;
guidata(hObject, handles);
uiresume;


% --- cancelPushButton
function cancelPushButton_Callback(hObject, eventdata, handles)
uiresume;

