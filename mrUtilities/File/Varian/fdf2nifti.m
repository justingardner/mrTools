% fdf2nifti.m
%
%        $Id:$ 
%      usage: [h d] = fdf2nifti(fdffilename,<verbose>,<headerOnly>)
%         by: justin gardner
%       date: 08/08/11
%    purpose: converts fdf file to a nifti file with no return arguments, save file as hdr/img pair
%
function [d hdr] = fdf2nifti(fdffilename,verbose,headerOnly)

d = [];hdr = [];
% check arguments
if ~any(nargin == [1 2 3])
  help fdf2nifti
  return
end

% set verbose default
if nargin < 2,verbose = 0;headerOnly = false;end
if nargin < 3,headerOnly = false;end

% see if we have a wildcard (which means we want to process a bunch of files
if ~isempty(strfind(fdffilename,'*')) && (nargout == 0)
  dirlist = dir(fdffilename);
  for i = 1:length(dirlist)
    fdf2nifti(dirlist(i).name,verbose);
  end
  return
end

% read the file
[d h] = readfdf(fdffilename,verbose,headerOnly);
if isempty(h),return,end

% Y and Z coordinates seems to be flipped and swapped for fdf
swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];
flipYZ = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1];

if 0
  % get the nifti header from the procpar
  [hdr info] = fid2niftihdr(setext(fdffilename,'fid'),verbose);

  % the fdf has X/Y swapped
  qform44 = hdr.qform44*swapXY;
  hdr = cbiSetNiftiQform(hdr,qform44);
end


% get orientation of center of slice around magnet coordinates
orientation = eye(4);
orientation(1,1:3) = h(1).orientation(1:3);
orientation(2,1:3) = h(1).orientation(4:6);
orientation(3,1:3) = h(1).orientation(7:9);
if verbose,disp(sprintf('(fdf2nifti) Orientation\n%s',mlrnum2str(orientation,'compact=0')));end

% get voxel size
voxelSize = diag([10*h(1).span./h(1).matrix h(1).roi(3)*10]);
voxelSize(4,1:4) = [0 0 0 1];
if verbose,disp(sprintf('(fdf2nifti) Voxel size: %s',mlrnum2str(diag(voxelSize)')));end

% get offset to center of slice
locationOffset = [1 0 0 h(1).location(1)*10;0 1 0 h(1).location(2)*10; 0 0 1 h(1).location(3)*10 ; 0 0 0 1];
if verbose,disp(sprintf('(fdf2nifti) Location offset: %s',mlrnum2str(locationOffset(1:3,4)')));end

% fix matrix for 2D images to have a 3rd dim
if length(h(1).matrix) < 3, h(1).matrix(3) = 1;end

sliceCenterOffset = eye(4);
sliceCenterOffset(1,4) = -h(1).matrix(1)/2;
sliceCenterOffset(2,4) = h(1).matrix(2)/2;
sliceCenterOffset(3,4) = h(1).matrix(3)/2;
if verbose,disp(sprintf('(fdf2nifti) Slice center offset: %s',mlrnum2str(sliceCenterOffset(1:3,4)')));end

% then get offset to first voxel
%originOffsetFromSliceCenter = [h(1).matrix(1)/2 h(1).matrix(2)/2 0 1]';
%originOffsetFromSliceCenter = orientation*originOffsetFromSliceCenter

qform44 = orientation*voxelSize*sliceCenterOffset*locationOffset*flipYZ*swapXY;
if verbose,disp(sprintf('(fdf2nifti) Qform44\n%s',mlrnum2str(qform44,'compact=0')));end

% we still need to debug all of this
disp(sprintf('(fdf2nifti) !!! Not sure the qform44 computed from the fdf header is correct (we need to debug this) !!!!'));

% create header
hdr = cbiCreateNiftiHeader; 
hdr.pixdim(2:5) = diag(voxelSize);
hdr.dim(1) = 4;
if length(h(1).matrix) == 2
  hdr.dim(2:3) = h(1).matrix(1:2);
  if isfield(h(1),'slices')
    hdr.dim(4) = h(1).slices;
  else
    hdr.dim(4) = 1;
  end
  hdr.dim(5:6) = 1;
else
  hdr.dim(2:4) = h(1).matrix(1:3);
  hdr.dim(5) = 1;
end
hdr = cbiSetNiftiQform(hdr,qform44);
