% mlrImportFreeSurfer.m
%
%        $Id$
%      usage: mlrImportFreeSurfer(varargin)
%         by: eli merriam
%       date: 07/11/07
%    purpose: 
%
function params = mlrImportFreeSurfer(varargin)

% check arguments
% if ~any(nargin == [0])
%   help mlrImportFreeSurfer
%   return
% end

if ~exist('mrParamsDialog')
  disp(sprintf('(mlrImportFreeSurfer) You must have mrTools in your path to run this'));
  return
end

if ~ispc
  mriConvert = 'mri_convert';
  [retval retstr] = system('which mri_convert');
  if retval == 1
    warnMessage = '(mlrImportFreeSurfer) Could not find FreeSurfer command mri_convert which is needed to convert the FreeSurfer anatomy file to a nifti file. This is usually in the bin directory under your freesurfer installation. You may need to install freesurfer and add that to your path. See instructions on wiki http://gru.stanford.edu/doku.php/mrTools/howTo#installation';
    if ispc
      warnMessage = [warnMessage '. This is probably because you''re running mrTools on a Windows PC.'];
    end
    mrWarnDlg(warnMessage,'Yes');
    mriConvert = [];
  end
else
  mriConvert = [];
end

% evaluate the arguments
eval(evalargs(varargin));
if ieNotDefined('defaultParams'),defaultParams = 0;end

% check directories
if ~isdir('surf'),disp(sprintf('(mlrImportFreeSurfer) Could not find surf directory'));return,end
if ~isdir('mri'),disp(sprintf('(mlrImportFreeSurfer) Could not find mri directory'));return,end

if ieNotDefined('baseName'), baseName = getLastDir(pwd);end

wmFile    = 'smoothwm';
gmFile    = 'pial';
infFile   = 'inflated';
curvFile  = 'curv';
anatFile  = 'T1.mgz';
hemi      = {'lh', 'rh'};
hemiNames = {'left', 'right'};
outDir    = fullfile(pwd,'surfRelax');
if ~exist('volumeCropSize')
  volumeCropSize = [176 256 256];
end
defaultPixelSize = [1 1 1];
if ~exist('pixelSize')
  pixelSize = defaultPixelSize;
end

paramsInfo = {...
    {'freeSurferDir',pwd,'directory where the freeSurfer files live'}, ...
    {'outDir', outDir,'directory that OFF surfaces will be written to'}, ...
    {'wmFile', wmFile, 'name of the surface defining the white/gray boundary'}, ...
    {'gmFile', gmFile, 'name of the surface defining the gray/pial boundary'}, ...
    {'curvFile', curvFile, 'name of the file specifing the curvature'}, ...
    {'infFile', infFile, 'name of the inflated surface'}, ...
    {'anatFile', anatFile, 'name of the base anatomy file'}, ...
    {'baseName', baseName, 'subject initials'}, ...
    {'volumeCropSize',volumeCropSize, 'Size to crop the volume anatomy to (in voxels).'},...
             };
if isempty(mriConvert) %if mri_convert wasn't found, ask for the resolution
  paramsInfo{end+1} = {'pixelSize',pixelSize,'Resolution of the volume anatomy (in mm). This is normally read from the converted anatomy file, but mri_convert is not available.'};
end

% get the parameters from the user
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo,'mlrImportFreeSurfer',1.5);
end

% return if user hit cancel
if isempty(params)
  return
end

if ~isdir(outDir)
  mkdir(outDir);
end


% convert the volume to mlr volume
anatFile = fullfile(params.freeSurferDir, 'mri', params.anatFile);
niftiExt = mrGetPref('niftiFileExtension');
switch niftiExt
  case '.nii'
    out_type='nii';
  case '.img'
    out_type='nifti1';
end
outFile = fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt));
if mlrIsFile(anatFile)
  commandString = sprintf('%s --out_type %s --out_orientation RAS %s %s',mriConvert,out_type,anatFile,outFile);
  if ~isempty(mriConvert)
    disp(sprintf('(mlrImportFreeSurfer) Converting the volume anatomy to nifti format'))
    % Because of a bug in mri_convert (well actually, not really a bug, freesurfer is not meant to deal with resolutions < 1mm),
    % we first need to convert without cropping to get the full dimensions of the freesurfer T1 volume
    setenv('LD_LIBRARY_PATH', '/usr/pubsw/packages/tiffjpegglut/current/lib:/opt/local/lib:/usr/local/lib:/opt/local/lib')
    system(commandString);
    %read the dimensions
    hdr = cbiReadNiftiHeader(outFile);
    volumeSize = hdr.dim(2:4);
    if all(iseven(volumeSize))
      %now run mri_convert again specifying the centre of the crop volume in the full volume
      % (otherwise the centre of the crop volume will be assumed to be (128,128,128) by mri_convert
      % + the qform/sform will be incorrect (this is the "bug")
      commandString = sprintf('%s --cropsize %i %i %i --crop %i %i %i',commandString,...
        params.volumeCropSize(1),params.volumeCropSize(2),params.volumeCropSize(3),...
        round(volumeSize(1)/2),round(volumeSize(2)/2),round(volumeSize(3)/2) ...
        );
    else
      % note that if the full volume size is an odd number of voxels, the qform/sform will be wrong
      % because the crop should happen at non-integer centre coordinates, which mri_convert does not allow
      % in this case, do not crop
      disp(sprintf('(mlrImportFreeSurfer) Odd number of voxels in original freesurfer anatomical volume (%s), not going to crop...',mat2str(volumeSize)));
      params.volumeCropSize = volumeSize;
    end
    system(commandString);
  else
    if mlrIsFile(outFile)
      fprintf('\n(mlrImportFreesurfer) Getting voxel and volume dimensions from existing %s file\n', strcat(params.baseName, '_', 'mprage_pp', niftiExt));
      hdr = cbiReadNiftiHeader(outFile);
      params.volumeCropSize = hdr.dim(2:4);
      fprintf('Voxel dimensions = %s\n',mat2str(hdr.pixdim(2:4)))
      fprintf('Volume dimensions = %s\n',mat2str(hdr.dim(2:4)))
    else
      disp(sprintf('\n(mlrImportFreeSurfer) To convert the canonical anatomy from Freesurfer to NIFTI format, run:  \n\t mri_convert%s \nin the appropriate terminal\n',commandString));
      disp('Note that if the original freesurfer anatomical volume is not ...x256x256 x 1mm, there will likely be a mismatch between the volume and the surfaces.');
      disp('If you know the dimensions of the freesurfer volume (in voxels), you are strongly encouraged to re-run mlrImportFreesurfer with these values as the cropSize field and to check for any mismatch between the converted volume and the converted surfaces.\n');
      mrWarnDlg(sprintf('(mlrImportFreeSurfer) !!!! Canonical anatomy not created !!!!'));
    end
  end
  
end

if ~mlrIsFile(outFile)
  if fieldIsNotDefined(params,'pixelSize')
    disp('(mlrImportFreeSurfer) Could not determine voxel size. Assuming 1x1x1 mm.');
    pixelSize = defautPixelSize;
  else
    pixelSize = params.pixelSize;
  end
else
  pixelSize=hdr.pixdim(2:4);
end
%   h = cbiReadNiftiHeader(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt)));
%   h.qform44(1,4) = -(h.dim(2)-1)/2;
%   h = cbiSetNiftiQform(h, h.qform44);
%   h = cbiSetNiftiSform(h, h.qform44);
%   cbiWriteNiftiHeader(h, fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt)));

% import the white and gray matter surface, as well as the inflated surface and curvature
disp(sprintf('(mlrImportFreeSurfer) Converting FreeSurfer surfaces to OFF format'))
for i = 1:length(hemi)
  % convert inner surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.wmFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_WM.off'));
  if mlrIsFile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize, pixelSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find inner (white matter) surface %s',getLastDir(surfFile,2)));
  end

  % convert outer surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.gmFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_GM.off'));
  if mlrIsFile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize, pixelSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find outer (pial/gray matter) surface %s',getLastDir(surfFile,2)));
  end

  % convert inflated surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.infFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Inf.off'));
  if mlrIsFile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize, pixelSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find inflated surface %s',getLastDir(surfFile,2)));
  end

  curvFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.curvFile));
  if mlrIsFile(curvFile)
    [curv, fnum] = freesurfer_read_curv(curvFile);
    saveVFF(fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Curv.vff')), -curv)
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find curvature file %s',getLastDir(curvFile)));
  end
end


return;

disp(sprintf('(mlrImportFreeSurfer) Constructing the gray graph'))
setenv('DYLD_LIBRARY_PATH', '/Users/eli/src/TFI/sw/lib/')
for i = 1:length(hemi)
    system(sprintf('surf2graph -verbose -gimage %s -gmin 40 -gmax 120 -expand -1 -gsurface %s -gray %s -savelabels %s %s', ...
                   strcat(params.baseName, '_', 'mprage_pp.hdr'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_GM.off'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '.Gray'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_graylabels.img'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_WM.off')...
                   ));
end


return;
% disp(sprintf('(mlrImportFreeSurfer) Fixing nifti header so that the volume is centered'))
% h = mlrImageReadNiftiHeader(sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));
% h.qform44(1,4) = (h.dim(2)-1)/2;
% h.qoffset_x = -1*(h.dim(2)-1)/2;
% h.sform_code = 0;
% cbiWriteNiftiHeader(h, sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));

