% fdf2nifti.m
%
%        $Id:$ 
%      usage: [h d] = fdf2nifti(fdffilename,<verbose>)
%         by: justin gardner
%       date: 08/08/11
%    purpose: converts fdf file to a nifti file with no return arguments, save file as hdr/img pair
%
function [d hdr] = fdf2nifti(fdffilename,verbose)

d = [];hdr = [];
% check arguments
if ~any(nargin == [1 2])
  help fdf2nifti
  return
end

% set verbose default
if nargin < 2,verbose = 0;end

% see if we have a wildcard (which means we want to process a bunch of files
if ~isempty(strfind(fdffilename,'*')) && (nargout == 0)
  dirlist = dir(fdffilename);
  for i = 1:length(dirlist)
    fdf2nifti(dirlist(i).name,verbose);
  end
  return
end

% read the file
[d h] = readfdf(fdffilename,verbose);
if isempty(h),return,end

% read procpar
procpar = readprocpar(setext(fdffilename,'fid'));
if procpar.ni == 1
  procpar.ni = procpar.np/2;
  if verbose,disp(sprintf('(fdf2nifti) Setting procpar.ni to: %i',procpar.ni));end
end

% get the nifti header from the procpar
[hdr info] = fid2niftihdr(procpar,verbose);

hdr.qform44

% get orientation of center of slice around magnet coordinates
orientation = eye(4);
orientation(1,1:3) = h(1).orientation(7:9);
orientation(2,1:3) = h(1).orientation(4:6);
orientation(3,1:3) = h(1).orientation(1:3);
orientation

% get voxel size
voxelSize = diag([10*h(1).span./h(1).matrix h(1).roi(3)*10]);
voxelSize(4,1:4) = [0 0 0 1];

% get offset to center of slice
sliceCenterOffset = [1 0 0 h(1).location(1)*10;0 1 0 h(1).location(2)*10; 0 0 1 h(1).location(3)*10 ; 0 0 0 1];
sliceCenterOffset = eye(4);

sliceCenterOffset(1,4) = info.dim(1)/2;
sliceCenterOffset(2,4) = info.dim(2)/2;
sliceCenterOffset(3,4) = info.dim(3)/2;
% then get offset to first voxel
%originOffsetFromSliceCenter = [h(1).matrix(1)/2 h(1).matrix(2)/2 0 1]';
%originOffsetFromSliceCenter = orientation*originOffsetFromSliceCenter

swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];


qform44 = orientation*voxelSize*sliceCenterOffset;
%keyboard
%hdr = cbiSetNiftiQform(hdr,qform44);

hdr.qform44
if nargout == 0
  niftiFilename = setext(fdffilename,'nii');
  if verbose, disp(sprintf('(fdf2nifti) Saving nifti file %s',niftiFilename)); end
  cbiWriteNifti(niftiFilename,d,hdr);
  d = [];
end
