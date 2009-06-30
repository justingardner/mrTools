function varargout = checkAlignmentGUI(varargin)
% CHECKALIGNMENTGUI M-file for checkAlignmentGUI.fig

% Last Modified by GUIDE v2.5 07-Jan-2007 12:07:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkAlignmentGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @checkAlignmentGUI_OutputFcn, ...
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


% --- Executes just before checkAlignmentGUI is made visible.
function checkAlignmentGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Parse varargin
handles.vol1 = varargin{1};
handles.vol2 = varargin{2};
handles.vol3 = varargin{3};
% Set nSlices in slider
nSlices = size(handles.vol1,3);
set(handles.sliceSlider,'Value',1);
set(handles.sliceSlider,'Min',1);
set(handles.sliceSlider,'Max',nSlices);
set(handles.sliceSlider,'sliderStep',[1/(nSlices-1) 10/(nSlices-1)]);
% Initialize display
sliceSlider_Callback(handles.sliceSlider, [], handles);
% Choose default command line output for checkAlignmentGUI
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes checkAlignmentGUI wait for user response (see UIRESUME)
uiwait(handles.checkAlignmentFigure);


% --- Outputs from this function are returned to the command line.
function varargout = checkAlignmentGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
% Finally, close the GUI
close(handles.checkAlignmentFigure);


% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
% Select checkAlignmentGUI figure and display images corresponding to
% current slice.
figure(handles.checkAlignmentFigure)
slice = round(get(hObject,'Value'));
subplot(1,3,1);
imagesc(handles.vol1(:,:,slice));
axis('image'); axis('off'); colormap(gray(256)); 
subplot(1,3,2);
imagesc(handles.vol2(:,:,slice));
axis('image'); axis('off'); colormap(gray(256)); 
subplot(1,3,3);
imagesc(handles.vol3(:,:,slice));
axis('image'); axis('off'); colormap(gray(256)); 


% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% uiresume triggers corAnalGUI_OutputFcn
uiresume;


