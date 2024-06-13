% mlrImportNeuropythy.m
%
%      usage: [lh rh] = mlrImportNeuropythy(lh_retino,rh_retino)
%         by: justin gardner
%       date: 09/04/19
%    purpose: Import anatomy predicted retinotopy from neuropythy into
%             mlr rois and overlays. 
%
%             See: http://gru.stanford.edu/doku.php/mrtools/atlas#benson_visual_field_atlas
%
function [lh rh] = mlrImportNeuropythy(lh_retino,rh_retino,varargin)

% check arguments
if nargin < 2
  help mlrImportNeuropythy
  return
end

% parse arguments
getArgs(varargin,{'surfPath',[],'labels',{'V1','V2','V3','hV4','V01','V02','LO1','LO2','TO1','TO2','V3b','V3a'},'lh_labels','','rh_labels','','doTestInMLR=0','saveDir=[]'});

% check names
[tf lh_retino] = checkNames(lh_retino);if ~tf,return,end
[tf rh_retino] = checkNames(rh_retino);if ~tf,return,end

% load the retino files
lh = load(lh_retino);lh.name = lh_retino;
rh = load(rh_retino);rh.name = rh_retino;

% check that these structures have expected fields
if ~checkExpectedFields(lh,{'varea','eccen','angle','sigma'}),return,end
if ~checkExpectedFields(rh,{'varea','eccen','angle','sigma'}),return,end

% check labels
[tf lh.labels rh.labels] = checkLabels(labels,lh_labels,rh_labels);
if ~tf,return,end

% make the overlays
lh = makeOverlays(lh);
rh = makeOverlays(rh);

% get surface names
[lh.surfaces rh.surfaces] = mlrGetSurfaceNames('surfPath',surfPath);
surfPath = lh.surfaces.path;

% bring up in viewer
curpwd = pwd;
try
  cd(surfPath);
  % left surface
  lh.surfaces = mrSurfViewer(lh.surfaces.outerSurface,lh.surfaces.outerCoords,lh.surfaces.innerSurface,lh.surfaces.innerCoords,'','',lh.overlays,lh.overlayNames);
  if isempty(lh.surfaces)
    disp(sprintf('(mlrImportNeuropythy) Left overlays have been rejected. Quitting'));
    cd(curpwd);
    return;
  end
  % right surface
  rh.surfaces = mrSurfViewer(rh.surfaces.outerSurface,rh.surfaces.outerCoords,rh.surfaces.innerSurface,rh.surfaces.innerCoords,'','',rh.overlays,rh.overlayNames);
  if isempty(rh.surfaces)
    disp(sprintf('(mlrImportNeuropythy) Right overlays have been rejected. Quitting'));
    cd(curpwd);
    return;
  end
catch
  cd(curpwd);
  keyboard
end
cd(curpwd);

% set path in surfaces structure
lh.surfaces.path = surfPath;
rh.surfaces.path = surfPath;

% ok. Now make the ROIS
lh.rois = makeROIs(lh);
rh.rois = makeROIs(rh);

% test in MLR
if doTestInMLR
  mlrTestROIsInMLR({lh.rois{:} rh.rois{:}},{lh.surfaces rh.surfaces});
end

% save to a directory
if ~isempty(saveDir)
  saveROI(saveDir,{lh.rois{:} rh.rois{:}});
end

%%%%%%%%%%%%%%%%%%
%    makeROIs    %
%%%%%%%%%%%%%%%%%%
function roi = makeROIs(x)

% colors for rois
areaColors = hsv(length(x.areas));

roi = {};
for iArea = x.areas
  roi{iArea} = mlrMakeROIFromSurfaceVertices(x.varea==iArea,x.surfaces);
  % make a label
  if length(x.labels) >= iArea
    roi{iArea}.name = x.labels{iArea};
  else
    roi{iArea}.name = 'Unknown';
  end
  % set color
  roi{iArea}.color = areaColors(iArea,:);
end

%%%%%%%%%%%%%%%%%%%%%%
%    makeOverlays    %
%%%%%%%%%%%%%%%%%%%%%%
function x = makeOverlays(x)

% number of vertices
x.nVertex = length(x.varea);

% make visual area overlay
% first get all area labels
x.areas = setdiff(unique(x.varea),0);
% make a color table for each of these areas
areaColors = hsv(length(x.areas));
% addisgn name
x.overlayNames{1} = 'Visual area';
% default to nan for the color
x.overlays{1} = nan(x.nVertex,3);
% set all values in overlay approriately
x.visualAreaVertices = find(x.varea~=0);
x.overlays{1}(x.visualAreaVertices,:) = areaColors(x.varea(x.visualAreaVertices),:);

% make eccentricity overlay
x.overlayNames{2} = 'eccentricity';
x.overlays{2} = makeOverlay(x.eccen,x.visualAreaVertices,'parula',true);

% polar angle
x.overlayNames{3} = 'polar angle';
x.overlays{3} = makeOverlay(x.angle,x.visualAreaVertices,'hsv',false);

% sigma
x.overlayNames{4} = 'sigma';
x.overlays{4} = makeOverlay(x.sigma,x.visualAreaVertices,'parula',true);

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
  if ~mlrIsFile(name)
    disp(sprintf('(mlrImortNeuropythy) Could not find file %s',name));
    return
  end
end

tf = 1;
  
%%%%%%%%%%%%%%%%%%%%%
%    checkLabels    %
%%%%%%%%%%%%%%%%%%%%%
function [tf lh_labels rh_labels] = checkLabels(labels,lh_labels,rh_labels)

tf = 0;

% first check whether lh_labels was passed in
if ~isempty(lh_labels)
  % if it is a string, then it should be a filename
  if isstr(lh_labels)
    % replace tilde if it is there
    lh_labels = mlrReplaceTilde(lh_labels);
    % check if it is afile
    if ~mlrIsFile(lh_labels)
      disp(sprintf('(mlrImportNeuropythy) %s is not a text file with the area labels in it',lh_labels));
      return
    end
    % if so, read it
    lh_labels = textread(lh_labels,'%s');
  end
else
  % use default labels
  lh_labels = labels;
end

% now do the same for rh_labels
if ~isempty(rh_labels)
  % if it is a string, then it should be a filename
  if isstr(rh_labels)
    % replace tilde if it is there
    rh_labels = mlrReplaceTilde(rh_labels);
    % check if it is afile
    if ~mlrIsFile(rh_labels)
      disp(sprintf('(mlrImportNeuropythy) %s is not a text file with the area labels in it',lh_labels));
      return
    end
    % if so, read it
    rh_labels = textread(rh_labels,'%s');
  end
else
  rh_labels = labels;
end

% now make sure that the labels have l and r properly set
for i = 1:length(lh_labels)
  lh_labels{i} = sprintf('l%s',lh_labels{i});
end
for i = 1:length(rh_labels)
  rh_labels{i} = sprintf('r%s',rh_labels{i});
end

tf = 1;


%%%%%%%%%%%%%%%%%
%    saveROI    %
%%%%%%%%%%%%%%%%%
function saveROI(saveDir,roi)

% make into cell array
roi = cellArray(roi);

% make the directory if it does not exist
if ~isdir(saveDir)
  mkdir(saveDir);
end

% for each roi, go and save
for iROI = 1:length(roi)
  % get the save name
  savename = fixBadChars(roi{iROI}.name,[],{'.','_'});
  eval(sprintf('%s = roi{iROI};',savename));
  save(fullfile(saveDir,savename),savename);
end
