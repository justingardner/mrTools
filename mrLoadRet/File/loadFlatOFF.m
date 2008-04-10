% loadFlatOFF.m
%
%        $Id$	
%      usage: loadFlatOFF('flatPatch.off')
%         by: eli merriam
%       date: 09/11/07
%    purpose: load a surfRelax flatpatch, and everything that goes with it
%  
%    Three ways to use this code:
function base = loadFlatOFF(flatFile)

% check arguments
if ~any(nargin == [0 1 2])
  help loadFlatOFF
  return;
end
base = [];
if ~exist('mrParamsDialog');
  disp(sprintf('(loadFlatOFF) You must have mrTools in your path to run this'));
  return;
end

% init arguments
if nargin == 0
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax off file';'*lat.off','SurfRelax off flat file';'*.*','All files'};
  title = 'Choose flat OFF file';
  flatFile = getPathStrDialog(startPathStr,title,filterspec,'on');
  flatFile = cellArray(flatFile);
  flatFile = flatFile{1};
end

if isstr(flatFile);
  % check to see if we are passed in a file name
  if ~isfile(flatFile);
    disp(sprintf('(loadFlatOFF) %s does not exist', flatFile));
    return;
  end
  % load the with mrFlatViewer
  params = mrFlatViewer(flatFile);
  if isempty(params),return,end

  % this returns a param structure, which unfortunally has all of the wrong names
  % for now, just translate them...
  params.flatDir = params.flatPath;
  params.flatFileName = params.flatFile;
  params.innerFileName = params.inner;
  params.outerFileName = params.outer;
  params.curvFileName = params.curv;
  params.anatFileName = params.anatomy;
  % mrFlatViewer  will create the gFlatViewer global structure
  % we can then use that for a lot of the further processing
  global gFlatViewer
end

% could be passed in a params structure, rather than a flat file
% this is the case when loadFlatOFF is called from makeFlat
if isstruct(flatFile);
  params = flatFile;
  if isfile(fullfile(params.flatDir, params.flatFileName));
    flatFile = loadSurfOFF(fullfile(params.flatDir, params.flatFileName));
    surf.inner = loadSurfOFF(fullfile(params.flatDir, params.innerFileName));
    surf.outer = loadSurfOFF(fullfile(params.flatDir, params.outerFileName));
    surf.curv = loadVFF(fullfile(params.flatDir, params.curvFileName));
    anat.hdr = cbiReadNiftiHeader(fullfile(params.flatDir, params.anatFileName));
    % now do necessary world2array xformation here
    surf.outer = xformSurfaceWorld2Array(surf.outer,anat.hdr);
    surf.inner = xformSurfaceWorld2Array(surf.inner,anat.hdr);
  end
end

% just crap out if the params still doesn't exist
if ieNotDefined('params');
  return;
end


% get two additional parameters if they were not passed in before
if ~any(isfield(params, {'threshold', 'flatRes'}));
  paramsInfo = {};
  paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'Whether or not to threshold the flat patch'};
  paramsInfo{end+1} = {'flatRes', 2, 'incdec=[-1 1]', 'Factore by which the resolution of the flat patch is increased'};
  flatParams = mrParamsDialog(paramsInfo,'Flat patch parameters');
  % check for cancel
  if isempty(flatParams);
    return;
  end
else
  flatParams.threshold = params.threshold;
  flatParams.flatRes = params.flatRes;
end

% check to see if we got here from the flatViewer
% we need to translate a few variable names if we did
if ~ieNotDefined('gFlatViewer');
  flat.whichInx  = gFlatViewer.flat.patch2parent(:,2); 
  flat.locsFlat  = gFlatViewer.flat.vtcs;
  flat.curvature = gFlatViewer.curv(flat.whichInx);
  flat.locsInner = gFlatViewer.surfaces.inner.vtcs(flat.whichInx,:);
  flat.locsOuter = gFlatViewer.surfaces.outer.vtcs(flat.whichInx,:);
  flat.hdr       = gFlatViewer.anat.hdr;
elseif ~isempty(flatFile);
  % or maybe we got here from loading a flatFile
  flat.whichInx  = flatFile.patch2parent(:,2); 
  flat.locsFlat  = flatFile.vtcs;
  flat.locsInner = surf.inner.vtcs(flat.whichInx,:);
  flat.locsOuter = surf.outer.vtcs(flat.whichInx,:);
  flat.curvature = surf.curv(flat.whichInx);
  flat.hdr       = anat.hdr;
else
  disp(sprintf('(loadFlatOFF) cannot find paramters needed to make base anatomy'));
end  


% this X-Y swaping only changes the orientation of the image
% isn't a crucial step
flat.locsFlat = [flat.locsFlat(:,2) flat.locsFlat(:,1) flat.locsFlat(:,3)];

flat.minLocsFlat = min(flat.locsFlat);
flat.locsFlat(:,1) = flat.locsFlat(:,1) - flat.minLocsFlat(1) + 1;
flat.locsFlat(:,2) = flat.locsFlat(:,2) - flat.minLocsFlat(2) + 1;

imSize = round(max(flat.locsFlat));

x = flat.locsFlat(:,1);
y = flat.locsFlat(:,2);
xi = [1:(1/flatParams.flatRes):imSize(1)];
yi = [1:(1/flatParams.flatRes):imSize(2)]';

yi = flipud(yi);

warning off
flat.map = griddata(x,y,flat.curvature,xi,yi,'linear');

% grid the 3d coords
for i=1:3
  flat.baseCoordsInner(:,:,i) =  griddata(x,y, flat.locsInner(:,i), xi, yi,'linear');
  flat.baseCoordsOuter(:,:,i) =  griddata(x,y, flat.locsOuter(:,i), xi, yi,'linear');
end

warning on
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
base.name = getLastDir(params.flatFileName);

% Extract permutation matrix to keep track of slice orientation.
base.permutationMatrix = getPermutationMatrix(base.hdr);

base.coordMap.flatDir = params.flatDir;
base.coordMap.flatFileName = params.flatFileName;
base.coordMap.innerFileName = params.innerFileName;
base.coordMap.outerFileName = params.outerFileName;
base.coordMap.curvFileName = params.curvFileName;
base.coordMap.anatFileName = params.anatFileName;
base.coordMap.flatRes = flatParams.flatRes;
base.coordMap.threshold = flatParams.threshold;

% load all the flat maps into the base. We
% need to make all the flat images have
% the same width and height.
if flatParams.threshold
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
base.clip = [0 1.5];

clear global gFlatViewer;

