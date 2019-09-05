% mlrImportNeuropythy.m
%
%      usage: mlrImportNeuropythy(lh_retino,rh_retino)
%         by: justin gardner
%       date: 09/04/19
%    purpose: Import anatomy predicted retinotopy from neuropythy into
%             mlr rois and overlays. See: http://gru.stanford.edu/doku.php/mrtools/atlas
%
function retval = mlrImportNeuropythy(lh_retino,rh_retino)

% check arguments
if ~any(nargin == [2])
  help mlrImportNeuropythy
  return
end

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

mrSurfViewer('jg_left_GM.off','','','','','',lOverlays,lOverlayNames);
keyboard

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
overlays{2} = makeOverlay(x.eccen,visualAreaVertices,'parula');

% polar angle
overlayNames{3} = 'polar angle';
overlays{3} = makeOverlay(x.angle,visualAreaVertices,'hsv');

% sigma
overlayNames{4} = 'sigma';
overlays{4} = makeOverlay(x.sigma,visualAreaVertices,'parula');

%%%%%%%%%%%%%%%%%%%%%
%    makeOverlay    %
%%%%%%%%%%%%%%%%%%%%%
function overlay = makeOverlay(data,visualAreaVertices,colorMapName)

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
overlay(visualAreaVertices,:) = colorMap(max(1,round(nColors*(data(visualAreaVertices)-minVal)/valRange)),:);


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
  

