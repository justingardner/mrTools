% importFStoSR.m
%
%        $Id$
%      usage: importFStoSR(varargin)
%         by: eli merriam
%       date: 07/11/07
%    purpose: 
%
function retval = importFStoSR(varargin)

% check arguments
% if ~any(nargin == [0])
%   help importFStoSR
%   return
% end

if ~exist('mrParamsDialog')
    disp(sprintf('(mglDoRetinotopy) You must have mrTools in your path to run this'));
    return
end

% evaluate the arguments
eval(evalargs(varargin));

if ieNotDefined('baseName'), baseName = [];end

wmFile    = 'smoothwm';
gmFile    = 'pial';
infFile   = 'inflated';
curvFile  = 'curv';
anatFile  = 'T1.mgz';
hemi      = {'lh', 'rh'};
hemiNames = {'left', 'right'};


paramsInfo = {...
    {'fsDir',    fullfile(getenv('SUBJECTS_DIR'), baseName),    'directory where the freeSurfer files live'}, ...
    {'outDir',   pwd,      'directory that OFF surfaces will be writeen to'}, ...
    {'wmFile',   wmFile,   'name of the surface defining the white/gray boundry'}, ...
    {'gmFile',   gmFile,   'name of the surface defining the gray/pial boundry'}, ...
    {'curvFile', curvFile, 'name of the file specifing the curvature'}, ...
    {'infFile',  infFile,  'name of the inflated surface'}, ...
    {'anatFile', anatFile, 'name of the base anatomy file'}, ...
    {'baseName', baseName,       'subject initials'}, ...
             };

% get the parameters from the user
params = mrParamsDialog(paramsInfo);

% return if use rr hit cancel
if isempty(params)
  return
end

% import the white and gray matter surface, as well as the inflated surface and curvature
disp(sprintf('(importFStoSR) Converting FreeSurfer surfaces to OFF format'))
for i = 1:length(hemi)
    fs2off(fullfile(params.fsDir, 'surf', strcat(hemi{i}, '.', params.wmFile)), ...
           fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_WM.off')));

    fs2off(fullfile(params.fsDir, 'surf', strcat(hemi{i}, '.', params.gmFile)), ...
           fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_GM.off')));

    fs2off(fullfile(params.fsDir, 'surf', strcat(hemi{i}, '.', params.infFile)), ...
           fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Inf.off')));
    
    [curv, fnum] = freesurfer_read_curv(fullfile(params.fsDir, 'surf', strcat(hemi{i}, '.', params.curvFile)));
    eval(sprintf('save %s %s -V6;', ...
                 fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Curv.mat')), 'curv'));
end


disp(sprintf('(importFStoSR) Converting the volume anatomy to nifti format'))
setenv('LD_LIBRARY_PATH', '/usr/pubsw/packages/tiffjpegglut/current/lib:/opt/local/lib:/usr/local/lib:/opt/local/lib')
system(sprintf('mri_convert --out_type nifti1 --out_orientation RAS --cropsize 176 256 256 %s %s', ...
               fullfile(params.fsDir, 'mri', params.anatFile), ...
               fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))))

h = cbiReadNiftiHeader(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr')));
h.qform44(1,4) = -(h.dim(2)-1)/2;
h = cbiSetNiftiQform(h, h.qform44);
h = cbiSetNiftiSform(h, h.qform44);
cbiWriteNiftiHeader(h, fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr')));


return;

disp(sprintf('(importFStoSR) Constructing the gray graph'))
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
% disp(sprintf('(importFStoSR) Fixing nifti header so that the volume is centered'))
% h = cbiReadNiftiHeader(sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));
% h.qform44(1,4) = (h.dim(2)-1)/2;
% h.qoffset_x = -1*(h.dim(2)-1)/2;
% h.sform_code = 0;
% cbiWriteNiftiHeader(h, sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));