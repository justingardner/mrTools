% mlrImportFreesurferLabel.m
%
%      usage: roi = mlrImportFreesurferLabel(labelFilename)
%         by: justin gardner
%       date: 09/05/19
%    purpose: Imports a freesurfer label as an roi. 
%             roi = mlrImportFreesurfer('lh_.V1_exvivo.label');
%
%             You can also load a whole directory of ROIs
%             roi = mlrImportFreesurfer('~/data/atlas/freesurfer/jg/label');
%
%             You can also have it show you the ROIS in a temporary MLR
%             roi = mlrImportFreesurfer('~/data/atlas/freesurfer/jg/label','doTestInMLR=1');
%
%             See: http://gru.stanford.edu/doku.php/mrtools/atlas#freesurfer_labels
%
function roi = mlrImportFreesurferLabel(filename,varargin)

% check arguments
if nargin < 1
  help mlrImportLabel
  return
end

roi = [];

% parse arguments
getArgs(varargin,{'leftSurfaceNames=[]','rightSurfaceNames=[]','hemi=[]','doTestInMLR=0','saveDir=[]'});

% check if this is a directory, and if so import all labels found in directory
if isdir(filename)
  % get directory listing
  listing = dir(filename);
  labelNames = {};
  % look for label files in directory
  for iFile = 1:length(listing)
    % see if it is a label
    if strcmp('label',getext(listing(iFile).name));
      labelNames{end+1} = listing(iFile).name;
    end
  end

  % see if we found any
  if isempty(labelNames)
    disp(sprintf('(mlrImportFreesurferLabel) Could not find any .label files in directory: %s',filename));
    return
  end
  % if we found some, then ask get the surface names
  if isempty(leftSurfaceNames) || isempty(rightSurfaceNAmes)
    [leftSurfaceNames rightSurfaceNames] = mlrGetSurfaceNames;
  end
  % now recursively run to get each roi
  colors = hsv(length(labelNames));
  for iROI = 1:length(labelNames)
    % display what we are doing
    disp(sprintf('(mlrImportFreesurferLabel) Importing %s',labelNames{iROI}));
    % import
    roi{iROI} = mlrImportFreesurferLabel(fullfile(filename,labelNames{iROI}),'leftSurfaceNames',leftSurfaceNames,'rightSurfaceNames',rightSurfaceNames);
    % set the color
    roi{iROI}.color = colors(iROI,:);
  end
  % test in MLR
  if doTestInMLR
    mlrTestROIsInMLR(roi,{leftSurfaceNames rightSurfaceNames});
  end
  % if we are asked to save
  if ~isempty(saveDir),saveROI(saveDir,roi),end
  % all done
  return
end

% check the file
if ~isfile(filename)
  disp(sprintf('(mlrImportFreesurferLabel) Could not find file %s',filename));
  return
end


% open file
f = fopen(filename);

% get header
header = fgetl(f);

% get number of vertices
nVertices = str2num(fgetl(f));
% get each vertex row-by-row
for iVertex = 1:nVertices
  vertex(iVertex,:) = str2num(fgetl(f));
end

% close file
fclose(f);


% decide if this is left or right
if isempty(hemi)
  splitFilename = strsplit(getLastDir(filename),'.');
  if (length(splitFilename) < 1) || ~any(strcmp(splitFilename{1},{'lh','rh'}))
    % ask user which hemisphere this is
    paramsInfo{1} = {'hemisphere',{'right','left'}};
    params = mrParamsDialog(paramsInfo,sprintf('Which hemisphere is roi %s for?',filename));
    if isempty(params),return,end
    if isequal(params.hemisphere,'right')
      hemi = 'rh';
    else
      hemi = 'lh';
    end
  else
    hemi = splitFilename{1};
  end
end

% choose which surface names to use
if strcmp(hemi,'rh')
  surfaceNames = rightSurfaceNames;
  if isempty(surfaceNames)
    [~,surfaceNames] = mlrGetSurfaceNames;
  end
else
  surfaceNames = leftSurfaceNames;
  if isempty(surfaceNames)
    [surfaceNames,~] = mlrGetSurfaceNames;
  end
end

% load the surface
surf = loadSurfOFF(fullfile(surfaceNames.path,surfaceNames.outerCoords));
roiVertices = zeros(1,surf.Nvtcs);
roiVertices(vertex(:,1)+1) = 1;

% make the roi
roi = mlrMakeROIFromSurfaceVertices(roiVertices,surfaceNames);
roi.name = getLastDir(filename);

if ~isempty(saveDir),saveROI(saveDir,roi),end

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
  eval(sprintf('%s = roi{iROI}',savename));
  save(fullfile(saveDir,savename),savename);
end
