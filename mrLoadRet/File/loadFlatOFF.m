% loadFlatOFF.m
%
%        $Id$	
%      usage: loadFlatOFF('flatPatch.off')
%         by: eli merriam
%       date: 09/11/07
%    purpose: load a surfRelax flatpatch, and everything that goes with it
%  
%    Three ways to use this code:
function base = loadFlatOFF(flatFileName)

% check arguments
if ~any(nargin == [0 1 2])
  help loadFlatOFF
  return
end
base = [];
if ~exist('mrParamsDialog')
  disp(sprintf('(loadFlatOFF) You must have mrTools in your path to run this'));
  return
end

% init arguments
if nargin == 0
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax off file';'*lat.off','SurfRelax off flat file';'*.*','All files'};
  title = 'Choose flat OFF file';
  flatFileName = getPathStrDialog(startPathStr,title,filterspec,'on');
  % Aborted
  if isempty(flatFileName)
    disp(sprintf('(loadFlatOFF) loading flat patch aborted'));
    return
  end
  % make into a cell array
  flatFileName = cellArray(flatFileName);
end

% only take one flat file for now
if iscell(flatFileName)
  flatFileName =   flatFileName{1};
end

% if we are passed in a string, then this is a file name
if isstr(flatFileName)
  basename = sprintf('%s.off', stripext(flatFileName));
  % check for file
  if isfile(flatFileName)
    % load the file
    surf.flat = loadSurfOFF(flatFileName);
  else
    disp(sprintf('(loadFlatOFF) %s does not exist', flatFileName));
    return;
  end
  % check that it is indeed a patch
  if ~isfield(surf.flat, 'nPatch')
    disp(sprintf('(loadFlatOFF) %s is not a flat patch file', flatFileName))
    return;
  end
end

% guess which hemisphere
if regexp(flatFileName, 'left') 
  params.whichHemi = 'left';
elseif regexp(flatFileName, 'right') 
  params.whichHemi = 'right';
else
  disp('(loadFlatOFF) Cannot guess hemisphere.  The user must choose')
  params.whichHemi = {'left', 'right'};
end

% contents of the current directory
[params.flatDir params.flatFileName] = fileparts(flatFileName);
dirContents = dir(params.flatDir);


% grab the white matter file name from the parent name
if exist(fullfile(params.flatDir, surf.flat.parentSurfaceName), 'file')
  params.wmFileName{1} = surf.flat.parentSurfaceName;
else
  % if it the parent surface doesn't exist, guess the WM name
  disp(sprintf('(loadFloatPatch) %s does not exist.  Guessing the WM surface name', surf.flat.parentSurfaceName))
  params.wmFileName = {};
  for i=2:length(dirContents)
    if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'WM')) & (regexp(dirContents(i).name, '.off'))
      params.wmFileName{end+1} = dirContents(i).name;
    end
  end
end
% ATTN still need to add recursive directory searching


% guess the GM surface
params.gmFileName = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'GM')) & (regexp(dirContents(i).name, '.off'))
    params.gmFileName{end+1} = dirContents(i).name;
  end
end
% ATTN still need to add recursive directory searching


% need to guess the Curv surface
% this may be empty, in which case, we'll calculate one
params.curvFileName = {};
params.calcCurvFlag = 1;
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, params.whichHemi)) & (regexp(dirContents(i).name, 'Curv')) & (regexp(dirContents(i).name, '.vff'))
    params.curvFileName{end+1} = dirContents(i).name;
    params.calcCurvFlag = 0;
  end
end

% guess the anatomy file
% there may be many suitable anatomy files; list them all
params.anatFileName = {};
for i=2:length(dirContents)
  if (regexp(dirContents(i).name, '.hdr'))
    params.anatFileName{end+1} = dirContents(i).name;
  end
end

% parse the parameters
paramsInfo = {};
paramsInfo{end+1} = {'flatDir', params.flatDir, 'editable=0', 'directory path for the flat .off surface file'};
paramsInfo{end+1} = {'whichHemi', params.whichHemi, 'editable=0', 'the hemisphere that the patch comes from -- not really important'};
paramsInfo{end+1} = {'flatFileName', params.flatFileName, 'editable=0', 'name of the flat patch, either input on the command line or chosen by the file selector'};
paramsInfo{end+1} = {'wmFileName', params.wmFileName, 'name of the surface at the white/gray boundary (i.e., the inner surface'};
paramsInfo{end+1} = {'gmFileName', params.gmFileName, 'name of the surface at the gray/pial boundary (i.e., the outer surface)'};
if params.calcCurvFlag == 1;
  paramsInfo{end+1} = {'curvFileName', 'will calculate curvature on the fly', 'editable=0', 'name of the curvatue file'};
  paramsInfo{end+1} = {'calcCurvFlag', 1, 'type=checkbox', 'whether or not to recalculate the curvature on the file -- requires surffilt'};
else
  paramsInfo{end+1} = {'curvFileName', params.curvFileName, 'name of the curvature file'};
  paramsInfo{end+1} = {'calcCurvFlag', 0, 'type=checkbox', 'whether or not to recaculate the curvature on the fly -- requires surffilt'};
end
paramsInfo{end+1} = {'anatFileName', params.anatFileName, 'base anatomy file, must be a valid nifti file in register with the session'};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]', 'the resolution of the flat patch -- a value of 2 doubles the resolution'};
paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'thresholding the surface makes the background two-tone (binary curvature)'};

params = mrParamsDialog(paramsInfo, 'loadFlatPatch', []);

if isempty(params)
  return;
end

% load the rest of the surfaces
disppercent(-inf,'Loading surfaces');
[surf, params] = loadSurfHandler(surf, params);
disppercent(inf);

% calculate the base anatomy structure
disppercent(-inf,'Computing curvature image and coordinates');
base = calcFlatBase(surf, params);
disppercent(inf);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surf, params] = loadSurfHandler(surf, params)
% we have already loaded the flat patch
% now load the rest of the surfaces

% load the white matter surface
surf.wm = loadSurfOFF(fullfile(params.flatDir, params.wmFileName));

% load the gray matter surface
surf.gm = loadSurfOFF(fullfile(params.flatDir, params.gmFileName));

% load the curvature
if params.calcCurvFlag==1
  params.baseName = stripext(params.wmFileName);
  params.curvFileName = sprintf('%s_Curv.vff', params.baseName);
  % check for the SurfRelax program called 'surffilt'
  [status,result] = system('surffilt');
  if status ~= 0
    disp(sprintf('(loadFlatOFF) Could not run the SurfRelax program surffilt. Make sure that you have SurfRelax correctly installed on your system.'))
    return;
  else
    % calculate the curvature
    disppercent(-inf, sprintf('(loadFlatOFF) calculating surface curvature'));
    system(sprintf('surffilt -mcurv -iter 1 %s %s', ...
                   fullfile(params.flatDir, params.wmFileName), ...
                   fullfile(params.flatDir, params.curvFileName)));
    disppercent(inf);
  end
end
% read in the curvature file
[surf.curv, hdr] = loadVFF(fullfile(params.flatDir, params.curvFileName));
surf.curv = surf.curv';           % needs to be transposed;

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = cbiReadNifti(fullfile(params.flatDir, params.anatFileName));

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(surf.anat.hdr.qform44(1:3,1:3)));
surf.anat.permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base] = calcFlatBase(surf, params)
% creates a flat patch base anatomy
% returns a base, which can then be either saved or installed in MLR

flat.whichInx  = surf.flat.patch2parent(:,2);
flat.locsFlat  = surf.flat.vtcs;
flat.curvature = surf.curv(flat.whichInx);
flat.locsInner = surf.wm.vtcs(flat.whichInx,:);
flat.locsOuter = surf.gm.vtcs(flat.whichInx,:);
flat.hdr       = surf.anat.hdr;

% voxScale = [params.flatRes params.flatRes params.flatRes];
% flat.hdr = cbiSetNiftiQform(flat.hdr, flat.hdr.qform44*diag([1./voxScale 1]));
% flat.hdr = cbiSetNiftiSform(flat.hdr, flat.hdr.sform44*diag([1./voxScale 1]));


% this X-Y swaping only changes the orientation of the image
% isn't a crucial step
flat.locsFlat = [flat.locsFlat(:,2) flat.locsFlat(:,1) flat.locsFlat(:,3)];

% % make the lookup table
% % first get coordinates
% xmin = min(flat.locsFlat(:,1));
% xmax = max(flat.locsFlat(:,1));
% ymin = min(flat.locsFlat(:,2));
% ymax = max(flat.locsFlat(:,2));
% x = xmin:(1/params.flatRes):xmax;
% y = ymin:(1/params.flatRes):ymax;

% % now we need to interp the curvature on to a regular
% % grid. we do this by finding the x,y point that is closest
% % to the one on the grid.
% disppercent(-Inf, 'interpolating surface')
% for i = 1:length(x)
%   disppercent(i/length(x));
%   for j = 1:length(y)
%     % find nearest point in curvature
%     dist = (flat.locsFlat(:,1)-x(i)).^2+(flat.locsFlat(:,2)-y(j)).^2;
%     flat.pos(i,j) = first(find(min(dist)==dist));
%     % points that are greater than a distance of 5 away are
%     % probably not in the flat patch so mask them out
%     if (min(dist) < 5)
%       flat.mask(i,j) = 1;
%       flat.baseCoordsInner(i,j,:) = flat.locsInner(flat.pos(i,j),:);
%       flat.baseCoordsOuter(i,j,:) = flat.locsOuter(flat.pos(i,j),:);
%       flat.map(i,j) = flat.curvature(flat.pos(i,j), :);
%     else
%       flat.mask(i,j) = 0;
%       flat.baseCoordsInner(i,j,:) = [0 0 0];
%       flat.baseCoordsOuter(i,j,:) = [0 0 0];
%       flat.map(i,j) = 0;
%     end
%   end
% end
% disppercent(inf)

flat.minLocsFlat = min(flat.locsFlat);
flat.locsFlat(:,1) = flat.locsFlat(:,1) - flat.minLocsFlat(1) + 1;
flat.locsFlat(:,2) = flat.locsFlat(:,2) - flat.minLocsFlat(2) + 1;

imSize = round(max(flat.locsFlat));

x = flat.locsFlat(:,1);
y = flat.locsFlat(:,2);
xi = [1:(1/params.flatRes):imSize(1)];
yi = [1:(1/params.flatRes):imSize(2)]';

flat.map = griddata(x,y,flat.curvature,xi,yi,'linear');

% grid the 3d coords
for i=1:3
  flat.baseCoordsInner(:,:,i) =  griddata(x,y, flat.locsInner(:,i), xi, yi,'linear');
  flat.baseCoordsOuter(:,:,i) =  griddata(x,y, flat.locsOuter(:,i), xi, yi,'linear');
end

% mask out out-of-brain coords
flat.baseCoordsInner(isnan(flat.map)) = 0;
flat.baseCoordsOuter(isnan(flat.map)) = 0;

% get rid of any NaN's
flat.baseCoordsInner(~isfinite(flat.baseCoordsInner)) = 0;
flat.baseCoordsOuter(~isfinite(flat.baseCoordsOuter)) = 0;

% base.map = mrUpSample(base.map);

% now blur image
flat.blurMap(:,:) = blur(flat.map(:,:));
% threshold image
flat.thresholdMap(:,:) = (flat.blurMap(:,:)>nanmedian(flat.blurMap(:)))*0.5+0.25;
% flat.medianVal = nanmedian(flat.blurMap(:));
% flat.blurMap(flat.blurMap<flat.medianVal) = -1;
% flat.blurMap(flat.blurMap>flat.medianVal) = 1;
flat.thresholdMap = blur(flat.thresholdMap);
flat.thresholdMap(isnan(flat.map)) = 0;


% now generate a base structure
clear base;
base.hdr = flat.hdr;
[pathstr, base.name] = fileparts(params.flatFileName);
base.permutationMatrix = surf.anat.permutationMatrix;

base.coordMap.flatDir = params.flatDir;
base.coordMap.whichHemi = params.whichHemi;
base.coordMap.flatFileName = params.flatFileName;
base.coordMap.innerFileName = params.wmFileName;
base.coordMap.outerFileName = params.gmFileName;
base.coordMap.curvFileName = params.curvFileName;
base.coordMap.anatFileName = params.anatFileName;
base.coordMap.flatRes = params.flatRes;
base.coordMap.threshold = params.threshold;

% load all the flat maps into the base. We
% need to make all the flat images have
% the same width and height.
if params.threshold
  base.data(:,:,1) = flat.thresholdMap;
  base.range = [0 1];
  base.clip = [0 1];
else
  flat.map(flat.map>1) = 1;
  flat.map(flat.map<-1) = -1;
  flat.map = 32*(flat.map+1);
  base.data(:,:,1) = flat.map;
end    

base.coordMap.coords(:,:,1,1) = flat.baseCoordsInner(:,:,2);
base.coordMap.coords(:,:,1,2) = flat.baseCoordsInner(:,:,1);
base.coordMap.coords(:,:,1,3) = flat.baseCoordsInner(:,:,3);

base.coordMap.innerCoords(:,:,1,1) = flat.baseCoordsInner(:,:,2);
base.coordMap.innerCoords(:,:,1,2) = flat.baseCoordsInner(:,:,1);
base.coordMap.innerCoords(:,:,1,3) = flat.baseCoordsInner(:,:,3);

base.coordMap.outerCoords(:,:,1,1) = flat.baseCoordsOuter(:,:,2);
base.coordMap.outerCoords(:,:,1,2) = flat.baseCoordsOuter(:,:,1);
base.coordMap.outerCoords(:,:,1,3) = flat.baseCoordsOuter(:,:,3);

base.coordMap.dims = flat.hdr.dim([2 3 4])';

base.range = [min(min(base.data)) max(max(base.data))];
base.clip = base.range;

% save(savename, 'base');
% v = viewSet(v,'newBase',base);



