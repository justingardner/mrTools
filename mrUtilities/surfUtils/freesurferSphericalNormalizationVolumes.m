% function freesurferSphericalNormalizationVolumes(params,<'justGetParams'>)
%
%   goal: Transforms data in source volume(s) from source to destination space using Freesurfer's
%         spherical normalization between fsSourceSubj and fsDestSubj freesrufer subjects, and keeping only data
%         located within the cortical sheet, and saving in Nifti file destVol.
%         The destination space is the space of the surfRelax T1 volume,
%         unless optional argument destVolTemplate is specified.
%         Source and destination volumes are sampled using the space specified by their sform
%         rotation matrix, or their qform if the sform is not set.
%         By default, both hemispheres are transformed
%
%   input parameters:
%         params.sourceVol (mandatory): volume to convert (source)
%         params.fsSourceSubj (mandatory): freesurfer subject ID cooresponding to source volume
%         params.fsSourceSurfSuffix (optional): suffix to add to the surface file names (e.g. fsSourceSubj_left_GMsuffix.off) (default = '')
%         params.fsDestSubj (optional): Freesurfer subject ID corresponding to destination volume (default: 'fsaverage')
%         params.destVol (mandatory): name of destination file(default: surfRelax anatomical scan of destination Freesurfer subject)
%         params.destVolTemplate (optional): template volume for destination space (default: surfRelax anatomical scan of source Freesurfer subject)
%         params.sourceSurfRelaxVolume (optional): surfRelax anatomical scan corresponding to the source volume, in case there are several volumes
%                                  in the surfRelax folder (default: surfRelax anatomical scan of source Freesurfer subject)
%         params.destSurfRelaxVolume (optional): surfRelax anatomical scan corresponding to the destination volumes, in case there are several
%                                  volumes in the surfRelax folder (default: surfRelax anatomical scan of destination Freesurfer subject)
%         params.hemisphere (optional): 'left','right or 'both' (default: 'both')
%         params.interpMethod (optional): interpolation method for ressampling the source volume to the source surface (default from mrGetPref)
%         params.outputBinaryData (optional): whether to binarize resampled data if the source data were binary, useful for ROI masks (default = false)
%         params.recomputeRemapping (optional): whether to recompute the mapping between the source and destination surfaces (default: false)
%         params.cropDestVolume (optional): whether to crop the destination volume to the smallest volume containing the surfaces (default = true)
%         params.dryRun (optional): if true, just checks that the input files exist (default: false)
%
%   author: julien besle (24/07/2020)

function params = freesurferSphericalNormalizationVolumes(params,varargin)

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end

if ieNotDefined('params')
  params = struct;
end

if fieldIsNotDefined(params,'sourceVol')
  params.sourceVol = {''};
end
if fieldIsNotDefined(params,'fsSourceSubj')
  params.fsSourceSubj = '';
end
if fieldIsNotDefined(params,'fsSourceSurfSuffix')
  params.fsSourceSurfSuffix = '';
end
if fieldIsNotDefined(params,'fsDestSubj')
  params.fsDestSubj = 'fsaverage';
end
if fieldIsNotDefined(params,'destVol')
  params.destVol = {''};
end
if fieldIsNotDefined(params,'destVolTemplate')
  params.destVolTemplate = '';
end
if fieldIsNotDefined(params,'sourceSurfRelaxVolume')
  params.sourceSurfRelaxVolume = '';
end
if fieldIsNotDefined(params,'destSurfRelaxVolume')
  params.destSurfRelaxVolume = '';
end
if fieldIsNotDefined(params,'hemisphere')
  params.hemisphere = 'both';
end
if fieldIsNotDefined(params,'interpMethod')
  params.interpMethod = mrGetPref('interpMethod');
end
if fieldIsNotDefined(params,'outputBinaryData')
  params.outputBinaryData = false;
end
if fieldIsNotDefined(params,'recomputeRemapping')
  params.recomputeRemapping = false;
end
if fieldIsNotDefined(params,'cropDestVolume')
  params.cropDestVolume = true;
end
if fieldIsNotDefined(params,'dryRun')
  params.dryRun = false;
end

if justGetParams, return; end

corticalDepthStep = 0.1;
corticalDepths = 0:corticalDepthStep:1;
nDepths = length(corticalDepths);

if ischar(params.sourceVol)
  params.sourceVol = {params.sourceVol};
end
% load source volume data
nSources = length(params.sourceVol);
for iSource = 1:nSources
  if ~exist(params.sourceVol{iSource},'file')
    if ~params.dryRun
      mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Could not find source volume %s',params.sourceVol{iSource}));
      return;
    end
  else
    sourceHdr = mlrImageReadNiftiHeader(params.sourceVol{iSource});
    if isempty(sourceHdr)
      if ~params.dryRun
        return;
      end
    else
      %check that dimensions and sforms are identical
      xformMistmatch = false;
      if iSource==1
        sourceXform = getXform(sourceHdr);
      else
        if ~isequal(sourceXform,getXform(sourceHdr))
          mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Rotation matrices do not match between sources %s and %s',params.sourceVol{1},params.sourceVol{iSource}));
          xformMistmatch = true;
        end
      end
    end
  end
end
if ~params.dryRun && xformMistmatch
  return;
end

if fieldIsNotDefined(params,'destVol')
  params.destVol = cell(size(params.sourceVol));
elseif ischar(params.destVol)
  [~,~,extension] = fileparts(params.destVol);
  if isempty(extension)
    % we assume it is a suffix to add to the source vol name
    destSuffix = params.destVol;
    params.destVol = cell(1,nSources);
    for iSource = 1:nSources
      [path,file,extension] = fileparts(params.sourceVol{iSource});
      params.destVol{iSource}=[path,file,destSuffix,extension];
    end
  else
    params.destVol = {params.destVol};
  end
end
if length(params.destVol) ~= nSources
  mrWarnDlg('(freesurferSphericalNormalizationVolumes) The number of source and destination volumes must match');
  if ~params.dryRun
    return;
  end
end

% first (re)compute mapping between source and destination surfaces (if needed)
[sourcePath, destPath] = remapSurfaces(params.fsSourceSubj,params.fsDestSubj,params.recomputeRemapping,params.dryRun,params.fsSourceSurfSuffix);
if (isempty(sourcePath) || isempty(destPath)) && params.dryRun
  return;
end

% load source surfRelax volume hdr
if fieldIsNotDefined(params,'sourceSurfRelaxVolume')
  params.sourceSurfRelaxVolume = fullfile(sourcePath,'surfRelax',[params.fsSourceSubj '_mprage_pp.nii']);
  if ~exist(params.sourceSurfRelaxVolume,'file')
    params.sourceSurfRelaxVolume = fullfile(sourcePath,'surfRelax',[params.fsSourceSubj '_mprage_pp.img']);
    if ~exist(params.sourceSurfRelaxVolume,'file')
      sourceSurfRelaxVolumes = [dir(fullfile(sourcePath,'surfRelax','*.nii')); dir(fullfile(sourcePath,'surfRelax','*.img'))];
      switch(length(sourceSurfRelaxVolumes))
        case 0
          mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Could not find source surfRelax volume in %s',fullfile(destPath,'surfRelax')));
          if ~params.dryRun
            return;
          end
        case 1
          params.sourceSurfRelaxVolume = sourceSurfRelaxVolumes(1).name;
        otherwise
          mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Multiple volumes in %s, specify a source surfRelax volume',dir(fullfile(sourcePath,'surfRelax'))));
          if ~params.dryRun
            return;
          end
      end
    end
  end
elseif ~exist(params.sourceSurfRelaxVolume,'file')
  mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Could not find source surfRelax volume %s',params.sourceSurfRelaxVolume));
  if ~params.dryRun
    return;
  end
end
if ~params.dryRun
  sourceSurfRelaxHdr = mlrImageReadNiftiHeader(params.sourceSurfRelaxVolume);
  sourceSurfRelaxXform = getXform(sourceSurfRelaxHdr);
end

% load destination surfRelax volume hdr
if fieldIsNotDefined(params,'destSurfRelaxVolume')
  params.destSurfRelaxVolume = fullfile(destPath,'surfRelax',[params.fsDestSubj '_mprage_pp.nii']);
  if ~exist(params.destSurfRelaxVolume,'file')
    params.destSurfRelaxVolume = fullfile(destPath,'surfRelax',[params.fsDestSubj '_mprage_pp.img']);
    if ~exist(params.destSurfRelaxVolume,'file')
      destSurfRelaxVolumes = [dir(fullfile(destPath,'surfRelax','*.nii')); dir(fullfile(destPath,'surfRelax','*.img'))];
      switch(length(destSurfRelaxVolumes))
        case 0
          mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Could not find destination volume in %s',dir(fullfile(destPath,'surfRelax'))));
          if ~params.dryRun
            return;
          end
        case 1
          params.destSurfRelaxVolume = fullfile(destPath,'surfRelax',destSurfRelaxVolumes(1).name);
        otherwise
          mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Multiple volumes in %s, specify a destination surfRelax volume',dir(fullfile(destPath,'surfRelax'))));
          if ~params.dryRun
            return;
          end
      end
    end
  end
elseif ~exist(params.destSurfRelaxVolume,'file')
  mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Could not find destination surfRelax volume %s',params.destSurfRelaxVolume));
  if ~params.dryRun
    return;
  end
end
if ~params.dryRun
  destSurfRelaxHdr = mlrImageReadNiftiHeader(params.destSurfRelaxVolume);
  destSurfRelaxXform = getXform(destSurfRelaxHdr);
end

% load destination volume header
if fieldIsNotDefined(params,'destVolTemplate')
  params.destVolTemplate = params.destSurfRelaxVolume;
elseif ~exist(params.destVolTemplate,'file')
  mrWarnDlg(sprintf('(freesurferSphericalNormalizationVolumes) Cannot find destination volume template %s',params.destVolTemplate));
  if ~params.dryRun
    return;
  end
end
if ~params.dryRun
  destHdr = mlrImageReadNiftiHeader(params.destVolTemplate);
  destXform = getXform(destHdr);
  destHdr.datatype = 16; % make sure data get exported as float32 (single) and that NaNs get saved as NaNs
  uncroppedDestDims = destHdr.dim(2:4)';
end

if params.dryRun
  return;
end

% Load original source and remapped destination surface
switch(params.hemisphere)
  case 'both'
    side = {'left','right'};
  otherwise
    side = {params.hemisphere};
end
surfs = {'WM','GM'};
nSides = length(side);
if params.cropDestVolume
  cropBox = [inf -inf; inf -inf; inf -inf];
end
for iSide=1:nSides
  for iSurf = 1:2
    %get surfaces in OFF format
    sourceSurf{iSurf,iSide} = loadSurfOFF([sourcePath '/surfRelax/' params.fsSourceSubj '_' side{iSide} '_' surfs{iSurf} params.fsSourceSurfSuffix '.off']);
    destSurf{iSurf} = loadSurfOFF([sourcePath '/surfRelax/' params.fsSourceSubj '_' side{iSide} '_' surfs{iSurf} params.fsSourceSurfSuffix '_' params.fsDestSubj '.off']); % same surface mesh as source, but with destination coordinates
    % convert vertices coordinates to surfRelax volume array coordinates
    sourceSurf{iSurf,iSide} = xformSurfaceWorld2Array(sourceSurf{iSurf,iSide},sourceSurfRelaxHdr);
    destSurf{iSurf} = xformSurfaceWorld2Array(destSurf{iSurf},destSurfRelaxHdr);

    % subdivide meshes by adding face centroids  (this avoids missing voxels in the cortical sheet for most regions)
    sourceSurf{iSurf,iSide} = subdivideMesh(sourceSurf{iSurf,iSide},1);
    destSurf{iSurf} = subdivideMesh(destSurf{iSurf},1);
    
    % convert to source and destination array coordinates
    sourceSurf{iSurf,iSide}.vtcs = (sourceXform\sourceSurfRelaxXform*[sourceSurf{iSurf,iSide}.vtcs';ones(1,sourceSurf{iSurf,iSide}.Nvtcs)])';
    destSurf{iSurf}.vtcs = (destXform\destSurfRelaxXform*[destSurf{iSurf}.vtcs';ones(1,destSurf{iSurf}.Nvtcs)])';
    sourceSurf{iSurf,iSide}.vtcs = sourceSurf{iSurf,iSide}.vtcs(:,1:3);
    destSurf{iSurf}.vtcs = destSurf{iSurf}.vtcs(:,1:3);
  end
  if params.cropDestVolume
    % get original (not-remapped) destination surfaces and apply same transformations as above (except subdividing)
    originalDestGMsurf = loadSurfOFF([destPath '/surfRelax/' params.fsDestSubj '_' side{iSide} '_GM.off']);
    originalDestGMsurf = xformSurfaceWorld2Array(originalDestGMsurf,destSurfRelaxHdr);
    originalDestGMsurf.vtcs = (destXform\destSurfRelaxXform*[originalDestGMsurf.vtcs';ones(1,originalDestGMsurf.Nvtcs)])';
    originalDestGMsurf.vtcs = originalDestGMsurf.vtcs(:,1:3);
  end
    
  % compute intermediate depth coordinates for destination mesh
  destCoords = zeros(destSurf{1}.Nvtcs,3,nDepths);
  for iDepth = 1:nDepths
    destCoords(:,:,iDepth) = (1-corticalDepths(iDepth))*destSurf{1}.vtcs + corticalDepths(iDepth)*destSurf{2}.vtcs;
  end
  
  if params.cropDestVolume % find smallest volume including the outer surface
    cropBox(:,1) = min(cropBox(:,1),floor(min(originalDestGMsurf.vtcs))');
    cropBox(:,2) = max(cropBox(:,2),ceil(max(originalDestGMsurf.vtcs))');
  end

  % compute mapping between destination surface and destination volume
  destCoords = permute(destCoords,[1 4 3 2]);
  surf2volumeMap{iSide} = inverseBaseCoordMap(destCoords,uncroppedDestDims);

end
clearvars('destCoords'); %save memory

if params.cropDestVolume % crop destination volume
  destHdr.dim(2:4) = diff(cropBox,[],2)+1;
  cropXform = eye(4);
  cropXform(1:3,4) = -cropBox(:,1)+1;
  destHdr.qform44 = cropXform\destHdr.qform44;
  destHdr.sform44 = cropXform\destHdr.sform44;
end


% compute intermediate depth coordinates for source mesh
for iSide = 1:nSides
  sourceCoords{iSide} = zeros(sourceSurf{1,iSide}.Nvtcs,3,nDepths,nSides);
  for iDepth = 1:nDepths
    sourceCoords{iSide}(:,:,iDepth) = (1-corticalDepths(iDepth))*sourceSurf{1,iSide}.vtcs + corticalDepths(iDepth)*sourceSurf{2,iSide}.vtcs;
  end
  sourceCoords{iSide} = reshape(permute(sourceCoords{iSide},[1 4 3 2]),[],3);
end

% interpolate data to destination volume
for iSource = 1:nSources
  destData = nan(uncroppedDestDims);
  % get source volume data
  [sourceData,sourceHdr] = mlrImageReadNifti(params.sourceVol{iSource});
  dataAreBinary = isequal(unique(sourceData(~isnan(sourceData))),[0 1]');
  if dataAreBinary
    interpMethod = 'linear'; % nearest doesn't work well for binary masks
  else
    interpMethod = params.interpMethod;
  end
  for iSide = 1:nSides %for each hemisphere
    % get surface data from source volume
    surfData = interpn((1:sourceHdr.dim(2))',(1:sourceHdr.dim(3))',(1:sourceHdr.dim(4))',...
               sourceData,sourceCoords{iSide}(:,1),sourceCoords{iSide}(:,2),sourceCoords{iSide}(:,3),interpMethod);
    % transform surface data to destination volume
    thisData = applyInverseBaseCoordMap(surf2volumeMap{iSide},uncroppedDestDims,surfData);
    if dataAreBinary && params.outputBinaryData % re-binarize binary data
      binaryThreshold = 0; %empirical (and conservative) threshold
      thisData(thisData<=binaryThreshold)=0;
      thisData(thisData>binaryThreshold)=1;
    end
    destData(~isnan(thisData)) = thisData(~isnan(thisData)); % we assume left and right surfaces sample exclusive sets of voxels, which is not exactly true at the midline
  end
  % write out the data
  if isempty(params.destVol{iSource})
    [path,file,extension] = fileparts(params.sourceVol{iSource});
    [filename,pathname] = uiputfile(fullfile(path,[file '_' params.fsDestSubj extension]),'Volume save name');
    if ~isnumeric(filename)
      params.destVol{iSource} = fullfile(pathname,filename);
    else
      return;
    end
  end
  if params.cropDestVolume
    destData = destData(cropBox(1,1):cropBox(1,2),cropBox(2,1):cropBox(2,2),cropBox(3,1):cropBox(3,2));
  end
  mlrImageWriteNifti(params.destVol{iSource},destData,destHdr);
end


function surf = subdivideMesh(surf,n)

for i = 1:n
  % calculate face centroids
  nFaces = size(surf.tris,1);
  faceCentroids = squeeze(mean(reshape(surf.vtcs(surf.tris,:),nFaces,3,3),2));
  surf.vtcs = [surf.vtcs;faceCentroids];
  
  % replace old by new faces
  surf.tris = [surf.tris(:) reshape(circshift(surf.tris,1,2),[],1) reshape(repmat(surf.Nvtcs+(1:surf.Ntris)',1,3),[],1)];
  
  % update number of vertices and faces
  surf.Nvtcs = surf.Nvtcs + surf.Ntris;
  surf.Ntris = surf.Ntris * 3;
end


function xform = getXform(hdr)

if fieldIsNotDefined(hdr,'sform44')
  if fieldIsNotDefined(hdr,'qform44')
    keyboard
  else
    xform = hdr.qform44 * shiftOriginXform; % convert from Nifti 0-indexing to Matlab 1-indexing
  end
else
  xform = hdr.sform44 * shiftOriginXform; % convert from Nifti 0-indexing to Matlab 1-indexing
end


