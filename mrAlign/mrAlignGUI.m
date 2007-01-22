function varargout = mrAlignGUI(varargin)
% mrAlignGUI M-file for mrAlignGUI.fig
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mrAlignGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mrAlignGUI_OutputFcn, ...
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


% --- Executes just before mrAlignGUI is made visible.
function mrAlignGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global ALIGN

% Initialize ALIGN global variable
ALIGN.volumePath = [];
ALIGN.inplanePath = [];
ALIGN.xform = eye(4);
ALIGN.guiXform = eye(4);
ALIGN.volSize = [64 64 64];
ALIGN.inplaneSize = [64 64 64];

% Initialize GUI
mrAlignGUI('sagittalRadioButton_Callback',hObject, eventdata, handles);
set(handles.transposeButton,'Value',0);
set(handles.flipButton,'Value',0);
set(handles.overlayButton,'Value',1);
set(handles.transparencySlider,'Value',1);
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
refreshAlignDisplay(handles);

% Choose default command line output for mrAlignGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object deletion, before destroying properties.
function axes_DeleteFcn(hObject, eventdata, handles)
clear global ALIGN

% --- Outputs from this function are returned to the command line.
function varargout = mrAlignGUI_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on mouse press over axes background.
function axes_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on button press in sagittalRadioButton.
function sagittalRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',1);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',0);
[m,index] = max(ALIGN.volumePermutation * [1 0 0]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in coronalRadioButton.
function coronalRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',0);
set(handles.coronalRadioButton,'Value',1);
set(handles.axialRadioButton,'Value',0);
[m,index] = max(ALIGN.volumePermutation * [0 1 0]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in axialRadioButton.
function axialRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',0);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',1);
[m,index] = max(ALIGN.volumePermutation * [0 0 1]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in overlayButton.
function overlayButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes on button press in transposeButton.
function transposeButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes on button press in flipButton.
function flipButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
global ALIGN
slice = (round(get(hObject,'Value')));
ALIGN.coords(ALIGN.sliceOrientation) = slice;
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transparencySlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function transparencySlider_Callback(hObject, eventdata, handles)
global ALIGN
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transX_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transY_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transZ_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transZ_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotX_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotY_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotZ_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotZ_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadVolMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user to choose volume.
if (ispref('mrLoadRet','defaultAnatomyPath'))
	initPath = getpref('mrLoadRet','defaultAnatomyPath');
else
    initPath = pwd;
end
pathStr = getPathStrDialog(initPath,'Choose vAnatomy file','*.img');
ALIGN.volumePath = pathStr;

% Load volume and header
h = mrMsgBox('Loading volume. Please wait');
[vData,hdr] = cbiReadNifti(ALIGN.volumePath);
volumeDimension = length(size(vData));
mrCloseDlg(h);

% Handle 4D file
volumeDimension = length(size(vData));
if (volumeDimension == 4)
    buttonName = questdlg('You have selected a 4D file.',...
		'Options for 4D files',...
		'Mean','First frame','Cancel','Mean');
	switch buttonName
		case 'Mean'
			vData = mean(vData,4);
		case 'First frame'
			vData = vData(:,:,:,1);
		case 'Cancel'
			return
	end
end

% Warning if no (qform) alignment information in the header.
% qform is initialized to identity by default in cbiReadNiftiHeader.
if ~(hdr.qform_code)
    mrWarnDlg('No alignment information in the volume header.');
end

% Warning if no (sform) base coordinate frame in the header.
% sform is initialized to identity by default in cbiReadNiftiHeader.
if ~(hdr.sform_code)
    mrWarnDlg('No base coordinate frame in the volume header.');
end

% Extract permutation matrix to keep track of slice orientation
[q,r] = qr(inv(hdr.qform44(1:3,1:3)));
permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

% Update ALIGN structure and GUI
ALIGN.volume = vData;
ALIGN.volumeHdr = hdr;
ALIGN.volumePermutation = permutationMatrix;
ALIGN.volSize = size(ALIGN.volume);
ALIGN.volumeVoxelSize = hdr.pixdim([2,3,4]);
ALIGN.coords = min(ALIGN.coords,ALIGN.volSize);
ALIGN.volumeClip = clipRange(ALIGN.volume);

% If both inplane and volume are loaded, then use the sforms from each for
% the alignment. Otherwise, use identity.
if ~isempty(ALIGN.volumeHdr) & ~isempty(ALIGN.inplaneHdr)
	ALIGN.xform = inv(ALIGN.volumeHdr.sform44) * ALIGN.inplaneHdr.sform44;
else
    ALIGN.xform = eye(4);
end

% Refresh GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
sagittalRadioButton_Callback(hObject, eventdata, handles);
refreshAlignDisplay(handles);

function cRange = clipRange(image)
% Choose clipping based on histogram
histThresh = length(image(:))/1000;
[cnt, val] = hist(image(:),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
cRange = [clipMin,clipMax];

% --------------------------------------------------------------------
function loadInplaneMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user to choose inplanes. 
initPath = pwd;
pathStr = getPathStrDialog(initPath,'Choose inplane anatomy file','*.img');
ALIGN.inplanePath = pathStr;

% Load inplane file and header
h = mrMsgBox('Loading inplanes. Please wait');
[vData,hdr] = cbiReadNifti(ALIGN.inplanePath);
mrCloseDlg(h);

% Handle 4D file
inplaneDimension = length(size(vData));
if (inplaneDimension == 4)
    buttonName = questdlg('You have selected a 4D file.',...
		'Options for 4D files',...
		'Mean','First frame','Cancel','Mean');
	switch buttonName
		case 'Mean'
			vData = mean(vData,4);
		case 'First frame'
			vData = vData(:,:,:,1);
		case 'Cancel'
			return
	end
end

% Warning if no (qform) alignment information in the header.
% qform is initialized to identity by default in cbiReadNiftiHeader.
if ~(hdr.qform_code)
    mrWarnDlg('No alignment information in the inplane header.');
end

% Warning if no (sform) base coordinate frame in the header.
% sform is initialized to identity by default in cbiReadNiftiHeader.
if ~(hdr.sform_code)
    mrWarnDlg('No base coordinate frame in the inplane header.');
end

% Update ALIGN structure
ALIGN.inplanes = vData;
ALIGN.inplanesClip = clipRange(ALIGN.inplanes);
ALIGN.inplaneHdr = hdr;
ALIGN.inplaneSize = size(ALIGN.inplanes);
ALIGN.inplaneVoxelSize = hdr.pixdim([2,3,4]);

% If both inplane and volume are loaded, then use the sforms from each for
% the alignment. Otherwise, use identity.
if ~isempty(ALIGN.volumeHdr) & ~isempty(ALIGN.inplaneHdr)
	ALIGN.xform = inv(ALIGN.volumeHdr.sform44) * ALIGN.inplaneHdr.sform44;
else
	ALIGN.xform = eye(4);
end

% Refresh GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function saveAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

sform = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;
ALIGN.inplaneHdr = cbiSetNiftiSform(ALIGN.inplaneHdr,sform);
hdr = cbiWriteNiftiHeader(ALIGN.inplaneHdr,ALIGN.inplanePath);

% --------------------------------------------------------------------
function setBaseCoordinateFrameMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

sform = ALIGN.volumeHdr.qform44;
ALIGN.volumeHdr = cbiSetNiftiSform(ALIGN.volumeHdr,sform);
hdr = cbiWriteNiftiHeader(ALIGN.volumeHdr,ALIGN.volumePath);

% --------------------------------------------------------------------
function saveAlignToFileMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Extract sform
sform = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;

% Prompt user for filename(s)
pathStr = getPathStrDialog(pwd,'Choose one or more nifti files','*.img','on');
if ~iscell(pathStr)
	pathStr = {pathStr};
end

% Loop through files and add sform to the headers
for p = 1:length(pathStr)
	if exist(pathStr{p},'file')
		hdr = cbiReadNiftiHeader(pathStr{p});
		hdr = cbiSetNiftiSform(hdr,sform);
		hdr = cbiWriteNiftiHeader(hdr,pathStr{p});
	else
		mrWarnDlg(['File ',pathStr{p},' not found.']);
	end
end

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Imports xform from separate file
function importAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user for file and load it
pathstr = getPathStrDialog(pwd,'Choose alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('This alignment file appears to correspond to an older version of mrAlign. You may need to use "Import" from the "File" menu.');
end
% Error if xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
end

% Update ALIGN structure and GUI.
ALIGN.xform = xform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function composeAlignmentMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user for first alignment file and load it.
pathstr = getPathStrDialog(pwd,'Choose first alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('This alignment file appears to correspond to an older version of mrAlign.');
end
% Error if  xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
else
    xform1 = xform;
end

% Prompt user for second alignment file and load it.
pathstr = getPathStrDialog(pwd,'Choose second alignment file','*.mat');
load(pathstr);
if ~exist(pathstr,'file')
    return
end
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('This alignment file appears to correspond to an older version of mrAlign.');
end  
% Error if xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
else
    xform2 = xform;
end

% Compose the alignments. Update ALIGN structure and GUI.
ALIGN.xform = xform2 * xform1;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
% Imports alignment from mrAlign-4.2 or earlier
function importOldAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user and load it
pathstr = getPathStrDialog(pwd,'Choose alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Error if  rot, trans, and scaleFac aren't loaded from the file.
if ~exist('rot','var') | ~exist('trans','var') | ~exist('scaleFac','var')
    myErrorDlg('Invalid alignment file.');
end

% Use loaded rot, trans, and scaleFac to compute xform.
S1 = [diag(scaleFac(1,:)) zeros(3,1); 0 0 0 1];
S2 = [diag(scaleFac(2,:)) zeros(3,1); 0 0 0 1];
Mi = [rot trans(:); 0 0 0 1];
xform = S2*Mi*inv(S1);

% Convert from mrAlign-4.2 convention to current convention (permute the
% dimensions and flip two of them).
% *** This hasn't been exhaustively test. Possible bug is that volSize(2)
% and volSize(3) may need to be swapped below.
xform = [xform(3,:); xform(1,:); xform(2,:); xform(4,:)];
xform(2,:) = -xform(2,:);
xform(2,4) = xform(2,4) + ALIGN.volSize(2);
xform(3,:) = -xform(3,:);
xform(3,4) = xform(3,4) + ALIGN.volSize(3);

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = inv(shiftXform)* xform;

% Update ALIGN structure and GUI
ALIGN.xform = xform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function exportMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Exports xform to separate file
function exportAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

xform = ALIGN.guiXform * ALIGN.xform;
mrAlignVersion = ALIGN.version;
pathstr = putPathStrDialog(pwd,'Specify alignment file','*.mat');
% pathstr = [] if aborted
if ~isempty(pathstr)
    save(pathstr,'xform','mrAlignVersion');
end

% --------------------------------------------------------------------
% Exports alignment for mrLoadRet-3.1 or earlier 
function exportOldAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Extract xform
xform = ALIGN.guiXform * ALIGN.xform;

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = xform * shiftXform;

% Reverse of the logic in importAlignMenuItem_Callback (mrAlignGUI) to
% convert from the current convention back to mrAlign-4.2 convention.
% *** This hasn't been exhaustively test. Possible bug is that volSize(2)
% and volSize(3) may need to be swapped below.
xform(3,4) = xform(3,4) - ALIGN.volSize(3);
xform(3,:) = -xform(3,:);
xform(2,4) = xform(2,4) - ALIGN.volSize(2);
xform(2,:) = -xform(2,:);
xform = [xform(2,:); xform(3,:); xform(1,:); xform(4,:)];

% Get scaleFac from voxel sizes
% *** This hasn't been exhaustively test. Possible bug is that x and y
% voxel sizes need to be swapped.
scaleFac = [1./ALIGN.inplaneVoxelSize'; 1./ALIGN.volumeVoxelSize'];

% compute rot and trans from 4x4
b = (xform(1:3,4))';
A = xform(1:3,1:3);
trans = b ./ scaleFac(2,:);
rot = diag(scaleFac(2,:)) \ A / diag( 1./scaleFac(1,:));

% Label it with version number. 
mrAlignVersion = ['exported from mrAlign ',num2str(ALIGN.version)];

pathstr = putPathStrDialog(pwd,'Specify alignment file','*.mat');
% pathstr = [] if aborted
if ~isempty(pathstr)
    save(pathstr,'rot','trans','scaleFac','mrAlignVersion');
end

% --------------------------------------------------------------------
function writeTifMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify alignment file','*.tif');
% pathstr = [] if aborted
if ~isempty(pathstr)
	img = refreshAlignDisplay(handles);
	imwrite(img,pathstr,'tif');
end

% --------------------------------------------------------------------
function quitMenuItem_Callback(hObject, eventdata, handles)
clear global ALIGN
delete(handles.figure1);


% --------------------------------------------------------------------
function manualAlignmentMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function initializeIdentityMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Set xform to identity, but scaled by voxel sizes
% *** Not yet test/debugged ***
inplaneXform = eye(4);
inplaneXform(1,1) = 1 / ALIGN.inplaneVoxelSize(1);
inplaneXform(2,2) = 1 / ALIGN.inplaneVoxelSize(2);
inplaneXform(3,3) = 1 / ALIGN.inplaneVoxelSize(3);
volumeXform = eye(4);
volumeXform(1,1) = 1 / ALIGN.volumeVoxelSize(1);
volumeXform(2,2) = 1 / ALIGN.volumeVoxelSize(2);
volumeXform(3,3) = 1 / ALIGN.volumeVoxelSize(3);
ALIGN.xform = inv(volumeXform) * inplaneXform;

% Reset GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipXMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [-1 0 0 ALIGN.inplaneSize(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipYMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [1 0 0 0; 0 -1 0 ALIGN.inplaneSize(2); 0 0 1 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipZMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [1 0 0 0; 0 1 0 0; 0 0 -1 ALIGN.inplaneSize(3); 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function computeAlignmentMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function initializeMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Error if there's no alignment information in the header.
% This would happen if these were analyze, not nifti, files.
if isempty(ALIGN.volumeHdr.qform44)
    mrErrorDlg('No alignment information in the volume header.');
end
if isempty(ALIGN.inplaneHdr.qform44)
    mrErrorDlg('No alignment information in the inplane header.');
end

% Compute alignment by composing the xforms from the two nifti headers.
ALIGN.xform = inv(ALIGN.volumeHdr.qform44) * ALIGN.inplaneHdr.qform44;

% Reset GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function cropInplanesMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.crop = selectCropRegion(ALIGN.inplanes);

% --------------------------------------------------------------------
function invertContrastMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

mrErrorDlg('Invert contrast not yet implemented.');

ALIGN.inplanes = - ALIGN.inplanes;
ALIGN.inplanesClip = clipRange(ALIGN.inplanes);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function fineMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
if isempty(ALIGN.volume) | isempty(ALIGN.inplanes)
	mrWarnDlg('Load Volume and Load Inplanes before computing alignment');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('Initialize aligment or load a previously saved alignment before computing.');
	return
end

% Compute alignment
xform = ALIGN.guiXform * ALIGN.xform;
xform = computeAlignment(ALIGN.inplanes, ALIGN.volume, xform, ALIGN.crop, ALIGN.NIter);
ALIGN.xform = xform;

% Reset GUI and refresh display
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------

function coarseMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

global ALIGN
if isempty(ALIGN.volume) | isempty(ALIGN.inplanes)
	mrWarnDlg('Load Volume and Load Inplanes before computing alignment');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('Initialize aligment or load a previously saved alignment before computing.');
	return
end

% reduce images
wbh = mrMsgBox('Reducing volumes. Please wait...');
lpf=[.0625 .25 .375 .25 .0625]';
volFilt = convXYZsep(ALIGN.volume, lpf, lpf, lpf, 'repeat', 'same');
inpFilt = convXYZsep(ALIGN.inplanes, lpf, lpf, lpf, 'repeat', 'same');
volReduce = volFilt(1:2:size(volFilt,1),1:2:size(volFilt,2),1:2:size(volFilt,3));
inpReduce = inpFilt(1:2:size(inpFilt,1),1:2:size(inpFilt,2),1:2:size(inpFilt,3));
mrCloseDlg(wbh);

% reduce xform & crop
xform = ALIGN.guiXform * ALIGN.xform;
xform(1:3,4) = xform(1:3,4)/2;
crop = floor(ALIGN.crop/2);

% compute alignment & expand new xform
xform = computeAlignment(inpReduce, volReduce, xform, crop, ALIGN.NIter);

% expand xform
xform(1:3,4) = xform(1:3,4)*2;
ALIGN.xform = xform;

% Reset GUI and refresh display
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function mutualInformationMenuItem_Callback(hObject, eventdata, handles)

mrErrorDlg('Mutual information registration not yet implemented.');

global ALIGN
% Options
% opts = optimset('fminunc');
opts = optimset('MaxFunEvals',100,'TolFun',1e-2,'DiffMaxChange',1,'DiffMinChange',0.1);
% Initial values (rot in deg then trans in voxels)
% *** Should be initialized to current xform, extracted from:
%     ALIGN.xform * ALIGN.guiXform;
x0 = zeros(1,6);
% Search
x = fminunc(@mutualInformationFun,x0,opts);
ALIGN.xform = extractXform(x);
% Reset GUI and refresh display
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

function val = mutualInformationFun(x)
global ALIGN
display(x)
% Interpolate the volume
xform = extractXform(x);
[NyI NxI NzI] = size(ALIGN.inplanes);
interpolatedVolume = regInplanes(ALIGN.volume, NxI, NyI, NzI, xform);
% Crop
if ~isempty(ALIGN.crop)
    crop = ALIGN.crop;
    inpCrop = ALIGN.inplanes([crop(1,2):crop(2,2)],[crop(1,1):crop(2,1)],:);
    volCrop = interpolatedVolume([crop(1,2):crop(2,2)],[crop(1,1):crop(2,1)],:);
else
    inpCrop = ALIGN.crop;
    volCrop = interpolatedVolume;
end
% Compute mutual information
val = - mutualInformation(inpCrop,volCrop,1);
display(val)

function xform = extractXform(x)
global ALIGN
addXform = eye(4);
x(1:3) = x(1:3)*pi/180;
cosx = cos(x(1));		sinx = sin(x(1));
cosy = cos(x(2));		siny = sin(x(2));
cosz = cos(x(3));		sinz = sin(x(3));
addXform(1,1) = cosz*cosy+sinz*sinx*siny;
addXform(1,2) = sinz*cosy-cosz*sinx*siny;
addXform(1,3) = cosx*siny;
addXform(2,1) = -sinz*cosx;
addXform(2,2) = cosz*cosx;
addXform(2,3) = sinx;
addXform(3,1) = sinz*sinx*cosy-cosz*siny;
addXform(3,2) = -cosz*sinx*cosy-sinz*siny;
addXform(3,3) = cosx*cosy;
addXform(1,4) = x(4);
addXform(1,5) = x(5);
addXform(1,6) = x(6);
xform = addXform * ALIGN.xform;

% --------------------------------------------------------------------
function checkAlignmentMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function checkWithCorrectionMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
if isempty(ALIGN.volume) | isempty(ALIGN.inplanes)
	mrWarnDlg('Load Volume and Load Inplanes before checking alignment.');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('Load, Initialize, or Compute the alignment before checking it.');
	return
end
checkAlignment(1)

% --------------------------------------------------------------------
function checkWithoutCorrectionMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
if isempty(ALIGN.volume) | isempty(ALIGN.inplanes)
	mrWarnDlg('Load Volume and Load Inplanes before checking alignment.');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('Load, Initialize, or Compute the alignment before checking it.');
	return
end
checkAlignment(0)

% --------------------------------------------------------------------
function jointHistogramMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Interpolate the volume
[NyI NxI NzI] = size(ALIGN.inplanes);
h = mrMsgBox('Wait while interpolating the inplanes...'); drawnow
interpolatedVolume = regInplanes(ALIGN.volume, NxI, NyI, NzI, ALIGN.xform);
mrCloseDlg(h);
% Compute histogram and display it
histogram = flipud(hist2(ALIGN.inplanes,interpolatedVolume));
FF = figure;
imagesc(log(histogram));
axis('image'); colormap('gray'); axis('off');
disp('Press a key to continue...')
pause
close(FF)

