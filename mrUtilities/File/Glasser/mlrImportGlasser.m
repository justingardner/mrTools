% mlrImportGlasser.m
%
%      usage: roi = mlrImportGlasser(labelFilename)
%         by: akshay jagadeesh
%       date: 06/01/2021
%    purpose: Imports a glasser label as an roi. 
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
function mlrImportGlasser(filename, varargin)

if ieNotDefined('filename')
    filename = '/Users/gru/data/freesurfer/s0423/atlas/glasser.mat';
end


getArgs(varargin,{'doTestInMLR=0','saveDir=[]'});


% check the file
if ~isfile(filename)
  disp(sprintf('(mlrImportGlasser) Could not find file %s',filename));
  return
end

%% open file and fix some stuff
f = load(filename);
f.labels = f.labels(2:end,:);
f.rh = round(f.rh);
f.lh = round(f.lh);

%% Load surface names
[lh_surfaces,rh_surfaces] = mlrGetSurfaceNames;

%% Load surfaces
lh_surf = loadSurfOFF(fullfile(lh_surfaces.path,lh_surfaces.outerCoords));
rh_surf = loadSurfOFF(fullfile(rh_surfaces.path, rh_surfaces.outerCoords));

%% Now assign each vertex to the value specified in the loaded file above.
areaColors = hsv(size(f.labels,1));

% Go through each area and create a vertex map for that area alone.
lh_rois = {};
rh_rois = {};

for i = 1:size(f.labels,1)
    % LEFT 
    lh_roiVertices = zeros(1,lh_surf.Nvtcs);
    lh_roiVertices(f.lh == i) = 1;
    lh_rois{i} = mlrMakeROIFromSurfaceVertices(lh_roiVertices, lh_surfaces);
    lh_rois{i}.name = ['lh_' strtrim(f.labels(i,:))];
    lh_rois{i}.color = areaColors(i,:);
    
    % RIGHT
    rh_roiVertices = zeros(1, rh_surf.Nvtcs);
    rh_roiVertices(f.rh == i) = 1;
    rh_rois{i} = mlrMakeROIFromSurfaceVertices(rh_roiVertices, rh_surfaces);
    rh_rois{i}.name = ['rh_' strtrim(f.labels(i,:))];
    rh_rois{i}.color = areaColors(i,:);
end
  
%% Test
if doTestInMLR
    mlrTestROIsInMLR({lh_rois{:} rh_rois{:}},{lh_surfaces rh_surfaces});
end

%% Save
if ~isempty(saveDir)
    saveROI(saveDir,{lh_rois{:} rh_rois{:}});
end

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

