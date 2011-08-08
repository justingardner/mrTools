% fdf2nifti.m
%
%        $Id:$ 
%      usage: [h d] = fdf2nifti(fdffilename,<verbose>)
%         by: justin gardner
%       date: 08/08/11
%    purpose: converts fdf file to a nifti file with no return arguments, save file as hdr/img pair
%
function [d hdr] = fdf2nifti(fdffilename,verbose)

% check arguments
if ~any(nargin == [1 2])
  help fdf2nifti
  return
end

% set verbose default
if nargin < 2,verbose = 0;end

% read the file
[d h] = readfdf(fdffilename,verbose);
if isempty(h),return,end

% get the nifti header from the fid
hdr = fid2niftihdr(setext(fdffilename,'fid'),verbose);

if nargout == 0
  niftiFilename = setext(fdffilename,'nii');
  if verbose, disp(sprintf('(fid2nifti) Saving nifti file %s',niftiFilename)); end
  cbiWriteNifti(niftiFilename,d,hdr);
  d = [];
end

