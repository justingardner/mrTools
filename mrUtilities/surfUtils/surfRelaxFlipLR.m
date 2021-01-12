% function thisView = surfRelaxFlipLR(freeSurferID,<force>)
%
%   goal: flip X (left-Right) coordinates of surfRelax surfaces and corresponding volume
%         for use in averaging left and right hemisphere data (see e.g. mlrSphericalNormGroup.m)
%         This should be run after:
%             1. LR flipping a T1-weighted volume
%             2. Creating surfaces from this LR-flipped volume with Freesurfer
%             3. Converting these surfaces to surfRelax format using mlrImportFreeSurfer.m
%         Resulting surfaces (twice LR-flipped) will be spherically-normalized to the reverse hemispheres
%         of fsaverage (left surface normalized to right fsaverage surface and vice-versa),
%         enabling averaging of left and right normalized surface
%         Original surfaces and volumes (once LF-flipped) are saved in subfolder 'surfRelax/Original (once-flipped) files/'
%
%   input: - freeSurferID: ID of the L/R-flipped freesurfer subject/template to flip
%          - force (optional): if true, re-creates files even they already exist and saves old files in subfolder 'surfRelax/Old twice-flipped files'
%
%   author: julien besle (31/07/2020)

function surfRelaxFlipLR(freeSurferID,force)

if isunix || ismac
  freesurferSubjdir = getenv('SUBJECTS_DIR');
end
if ieNotDefined('force')
  force = false;
end

if ieNotDefined('freesurferSubjdir')
  freesurferSubjdir = mrGetPref('volumeDirectory');
  if isempty(freesurferSubjdir)
    mrWarnDlg('(surfRelaxFlipLR) Cannot find the location of Freesurfer subject directory. Check your MR preferences using mrGetPref, or set Freesurfer''s environment variable SUBJECTS_DIR (Linux or Mac)');
    return
  end
  fprintf('(surfRelaxFlipLR) Assuming that the Freesurfer subject directory is %s\n',freesurferSubjdir);
end

if isempty(dir(fullfile(freesurferSubjdir,freeSurferID)))
  mrWarnDlg(['(surfRelaxFlipLR) Freesurfer subject ' freeSurferID ' does not exist']);
  return
end

surfRelaxPath = fullfile(freesurferSubjdir,freeSurferID,'surfRelax');
if isempty(dir(surfRelaxPath))
  mrWarnDlg(['(surfRelaxFlipLR) surfRelax folder does not exist in Freesurfer subject ' freeSurferID '. You must first run mlrImportFreesurfer.']);
  return
end

side = {'left','right'};
surfs = {'GM','WM','Inf'};
niftiExt = mrGetPref('niftiFileExtension');
for iSide = 1:2
  for iSurf = 1:length(surfs)
    surfaceFile{iSide,iSurf} = sprintf('%s_%s_%s.off',freeSurferID,side{iSide},surfs{iSurf});
  end
  curvatureFile{iSide} = sprintf('%s_%s_Curv.vff',freeSurferID,side{iSide});
end
volumeFile = sprintf('%s_mprage_pp%s',freeSurferID,niftiExt);

cwd = pwd;
cd(surfRelaxPath);

onceFlippedFolder = 'Original (once-flipped) files';
twiceFlippedFolder = 'Old twice-flipped files';
if exist(onceFlippedFolder,'dir')
  if force
    fprintf('(surfRelaxFlipLR) Moving previous twice-flipped files to folder ''%s''\n',twiceFlippedFolder);
    mkdir(surfRelaxPath,twiceFlippedFolder);
    for iSide = 1:2
      for iSurf = 1:length(surfs)
        movefile(surfaceFile{iSide,iSurf},twiceFlippedFolder);
      end
      movefile(curvatureFile{iSide},twiceFlippedFolder);
    end
    movefile(volumeFile,twiceFlippedFolder);
  else
    mrWarnDlg(sprintf('(surfRelaxFlipLR) Freesurfer subject %s has already been flipped. Use optional argument <force> to recompute',freeSurferID));
    return;
  end
else
  fprintf('(surfRelaxFlipLR) Moving original (once-flipped) files to folder ''%s''\n',onceFlippedFolder);
  mkdir(surfRelaxPath,onceFlippedFolder);
  for iSide = 1:2
    for iSurf = 1:length(surfs)
      movefile(surfaceFile{iSide,iSurf},onceFlippedFolder);
    end
    movefile(curvatureFile{iSide},twiceFlippedFolder);
  end
  movefile(volumeFile,onceFlippedFolder);
end

% flip volume and surfaces
fprintf('(surfRelaxFlipLR) Left-right flipping surfRelax volume (%s)\n',volumeFile);
[volume,hdr] = mlrImageLoad(fullfile(onceFlippedFolder,volumeFile));
volume = flip(volume,1);
mlrImageSave(volumeFile,volume,hdr);
fprintf('(surfRelaxFlipLR) Left-right flipping surfRelax surfaces\n');
for iSide = 1:2
  for iSurf = 1:length(surfs)
    surf = loadSurfOFF(fullfile(onceFlippedFolder,surfaceFile{iSide,iSurf}));
    % convert from surfRelax to array coordinates
    world2array = mlrXFormFromHeader(mlrImageReadNiftiHeader(fullfile(onceFlippedFolder,volumeFile)),'world2array');
    surf.vtcs = world2array*[surf.vtcs'; ones(1,surf.Nvtcs)];
    surf.vtcs(1,:) = hdr.dim(1) + 1 - surf.vtcs(1,:);    % flip X coordinates
    % convert back to surfRelax coordinates
    surf.vtcs = world2array\surf.vtcs;
    surf.vtcs = surf.vtcs(1:3,:)';
    writeOFF(surf, surfaceFile{iSide,iSurf});
  end
  copyfile(fullfile(onceFlippedFolder,curvatureFile{iSide}),curvatureFile{iSide});
end

cd(cwd);
