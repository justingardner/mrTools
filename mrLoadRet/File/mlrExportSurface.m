% mlrExportSurface.m
%
%        $Id:$ 
%      usage: mlrExportSurface(v)
%         by: justin gardner
%       date: 12/23/18
%    purpose: Export surface to wavefront off format
%
function retval = mlrExportSurface(varargin)

% check arguments
if ~any(nargin == [1 2])
  help mlrExportSurface
  return
end

% check if first argument is not view, that
% means that we have been called from the gui
% so we got an hObject
if ishandle(varargin{1})
  % get view
  v = viewGet(getfield(guidata(varargin{1}),'viewNum'),'view');
elseif isview(varargin{1})
  v = varargin{1};
else
  disp(sprintf('(mlrExportForAnalysis) Call with a view'));
  return
end

% get filename to save under
[filename, pathname] = uiputfile({'*.obj','Wavefront obj (*.obj)'},'Export surface to Wavefront obj file',setext(viewGet(v,'baseName'),'obj'));
if isequal(filename,0),return,end

% get the surface
s = viewGet(v,'baseSurface');

% get the colors
vertexColors = squeeze(refreshMLRDisplay(v));

% save it
mlrExportOFF(fullfile(pathname,filename),s,vertexColors);

  
  