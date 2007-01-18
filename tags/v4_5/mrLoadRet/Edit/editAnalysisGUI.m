function varargout = editAnalysisGUI(varargin)
% view = editAnalysisGUI(view);
%
% GUI for editing analysis (change name and view group, function, date).
% Created by GUIDE.
%
% djh 7/2003

% Last Modified by GUIDE v2.5 03-Aug-2006 20:30:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editAnalysisGUI_OutputFcn, ...
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


% --- Executes just before editAnalysisGUI is made visible.
function editAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Parse varargin
% Must be called as follows:
%    editAnalysisGUI(view)
if length(varargin) < 1
	mrErrorDlg('editAnalysisGUI: invalid initialization argument');
end
view = varargin{1};
if ~isview(view)
	mrErrorDlg('editAnalysisGUI: invalid initialization argument');
end

% Current analysis
analysisNum = viewGet(view,'curAnalysis');
analysisName = viewGet(view,'analysisName',analysisNum);
analysisGroup = viewGet(view,'analysisGroupName',analysisNum);
analysisFunction = viewGet(view,'analysisFunction',analysisNum);
analysisDate = viewGet(view,'analysisDate',analysisNum);

% Initialize handles
handles.view = view;
handles.analysisNum = analysisNum;
set(handles.nameEditText,'String',analysisName);
set(handles.groupText,'String',analysisGroup);
set(handles.functionText,'String',analysisFunction);
set(handles.dateText,'String',analysisDate);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editAnalysisGUI wait for user response (see UIRESUME)
uiwait(handles.editAnalysisDialog);


% --- Outputs from editAnalysisGUI
function varargout = editAnalysisGUI_OutputFcn(hObject, eventdata, handles) 
% Returns view
if ~isempty(handles)
    varargout{1} = handles.view;
    close(handles.editAnalysisDialog);
else
    varargout{1} = [];
end


% --- nameEditText
function nameEditText_Callback(hObject, eventdata, handles)

function nameEditText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- viewParamsPushButton.
function viewParamsPushButton_Callback(hObject, eventdata, handles)
view = handles.view;
n = handles.analysisNum;
params = viewGet(view,'analysisParams',n);
guiFunction = viewGet(view,'analysisGuiFunction',n);
groupName = viewGet(view,'analysisGroupName',n);
% params = guiFunction('groupName',groupName,'params',params);
evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params);'];
eval(evalstring);


% --- okPushButton.
function okPushButton_Callback(hObject, eventdata, handles)
view = handles.view;
nameString = get(handles.nameEditText,'String');
view = viewSet(view,'analysisName',nameString,handles.analysisNum);
handles.view = view;
guidata(hObject, handles);
uiresume;


% --- cancelPushButton
function cancelPushButton_Callback(hObject, eventdata, handles)
uiresume;



