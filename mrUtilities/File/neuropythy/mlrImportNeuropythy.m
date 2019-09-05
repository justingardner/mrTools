% mlrImportNeuropythy.m
%
%      usage: mlrImportNeuropythy(lh_retino,rh_retino)
%         by: justin gardner
%       date: 09/04/19
%    purpose: Import anatomy predicted retinotopy from neuropythy into
%             mlr rois and overlays. See: http://gru.stanford.edu/doku.php/mrtools/atlas
%
function retval = mlrImportNeuropythy(lh_retino,rh_retino,varargin)

% check arguments
if ~any(nargin == [2])
  help mlrImportNeuropythy
  return
end

% parse arguments
getArgs(varargin,{'surfPath',[]});

% check names
[tf lh_retino] = checkNames(lh_retino);if ~tf,return,end
[tf rh_retino] = checkNames(rh_retino);if ~tf,return,end

% load the retino files
lh = load(lh_retino);lh.name = lh_retino;
rh = load(rh_retino);rh.name = rh_retino;

% check that these structures have expected fields
if ~checkExpectedFields(lh,{'varea','eccen','angle','sigma'}),return,end
if ~checkExpectedFields(rh,{'varea','eccen','angle','sigma'}),return,end

% make the overlays
[lOverlays lOverlayNames] = makeOverlays(lh);
[rOverlays rOverlayNames] = makeOverlays(rh);

% get volume directory
if isempty(surfPath)
  titleStr = 'Choose left gray matter (outer) surface for this subject';
  disp(sprintf('(mlrImportNeuropythy) %s',titleStr));
  surfPath = mlrGetPathStrDialog(mrGetPref('volumeDirectory'),titleStr,'*.off');
  if isempty(surfPath),return,end
end

% make sure this is a GM (not WM or Inf)
surfFilename = getLastDir(surfPath);
surfPath = fileparts(surfPath);
surfFilename = replaceStr(surfFilename,'WM.','GM.');
surfFilename = replaceStr(surfFilename,'Inf.','GM.');

% make this into a left
leftSurfFilename = replaceStr(surfFilename,'right','left');
rightSurfFilename = replaceStr(surfFilename,'left','right');

% bring up in viewer
curpwd = pwd;
try
  cd(surfPath);
  leftSurfaces = mrSurfViewer(leftSurfFilename,'','','','','',lOverlays,lOverlayNames);
  if isempty(leftSurfaces),
    disp(sprintf('(mlrImportNeuropythy) Left overlays have been rejected. Quitting'));
    cd(curpwd);
    return;
  end
  rightSurfaces = mrSurfViewer(rightSurfFilename,'','','','','',rOverlays,rOverlayNames);
  if isempty(rightSurfaces),
    disp(sprintf('(mlrImportNeuropythy) Right overlays have been rejected. Quitting'));
    cd(curpwd);
    return;
  end
catch
  cd(curpwd);
  keyboard
end

cd(curpwd);

keyboard

%%%%%%%%%%%%%%%%%%%%
%    replaceStr    %
%%%%%%%%%%%%%%%%%%%%
function s = replaceStr(s,searchstr,replacestr)

foundloc = strfind(s,searchstr);
if ~isempty(foundloc)
  s = sprintf('%s%s%s',s(1:(foundloc(1)-1)),replacestr,s((foundloc(1)+length(searchstr)):end));
end

%%%%%%%%%%%%%%%%%%%%%%
%    makeOverlays    %
%%%%%%%%%%%%%%%%%%%%%%
function [overlays overlayNames] = makeOverlays(x)

% number of vertices
nVertex = length(x.varea);

% make visual area overlay
% first get all area labels
areaLabels = setdiff(unique(x.varea),0);
% make a color table for each of these areas
areaColors = hsv(length(areaLabels));
% addisgn name
overlayNames{1} = 'Visual area';
% default to nan for the color
overlays{1} = nan(nVertex,3);
% set all values in overlay approriately
visualAreaVertices = find(x.varea~=0);
overlays{1}(visualAreaVertices,:) = areaColors(x.varea(visualAreaVertices),:);

% make eccentricity overlay
overlayNames{2} = 'eccentricity';
overlays{2} = makeOverlay(x.eccen,visualAreaVertices,'parula',true);

% polar angle
overlayNames{3} = 'polar angle';
overlays{3} = makeOverlay(x.angle,visualAreaVertices,'hsv',false);

% sigma
overlayNames{4} = 'sigma';
overlays{4} = makeOverlay(x.sigma,visualAreaVertices,'parula',true);

%%%%%%%%%%%%%%%%%%%%%
%    makeOverlay    %
%%%%%%%%%%%%%%%%%%%%%
function overlay = makeOverlay(data,visualAreaVertices,colorMapName,logtransform)

% get number of vertices
nVertex = length(data);

% set up colormap
nColors = 1024;
colorMap = eval(sprintf('%s(%i)',colorMapName,nColors));

% get min and max
minVal = min(data(visualAreaVertices));
maxVal = max(data(visualAreaVertices));
valRange = maxVal - minVal;

% set overlay
overlay = nan(nVertex,3);
% get values
vals = (data(visualAreaVertices)-minVal)/valRange;
% log transform if needed
if logtransform
  vals = log(vals+0.0001);
  vals = (vals-min(vals))/(max(vals)-min(vals));
end
% and set into color range
vals = max(1,round(vals*nColors));
overlay(visualAreaVertices,:) = colorMap(vals,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkExpectedFields    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = checkExpectedFields(x,expectedFields)

tf = 0;

% check
missingFields = setdiff(expectedFields,fieldnames(x));
if ~isempty(missingFields)
  % make string of what is missing
  missingFieldsStr = '';
  for iField = 1:length(missingFields)
    missingFieldsStr = sprintf('%s%s ',missingFieldsStr,missingFields{iField});
  end
  % display
  disp(sprintf('(mlrImportNeuropythy) Structure is missing expected fields: %s',x.name,missingFieldsStr));
  return
end

tf = 1;

%%%%%%%%%%%%%%%%%%%%
%    checkNames    %
%%%%%%%%%%%%%%%%%%%%
function [tf name] = checkNames(name)

tf = 0;

if ~isstr(name) 
  disp(sprintf('(mlrImortNeuropythy) Input should be the filename of the mat file saved from neuropythy'));
  return
else
  name = setext(name,'mat');
  if ~isfile(name)
    disp(sprintf('(mlrImortNeuropythy) Could not find file %s',name));
    return
  end
end

tf = 1;
  

