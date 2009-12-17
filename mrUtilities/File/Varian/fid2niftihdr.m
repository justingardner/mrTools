% fid2niftihdr.m
%
%        $Id:$ 
%      usage: hdr = fid2niftihdr(fidname,<verbose>)
%         by: justin gardner
%       date: 12/09/09
%    purpose: strips out of fid2nifti the functionality to create nifti headers from procpar
%
function [hdr info] = fid2niftihdr(fidname,verbose)

hdr = [];info = [];
if ieNotDefined('verbose'),verbose=1;end

% check arguments
if ~any(nargin == [1 2])
  help fid2niftihdr
  return
end

% set fid extension
fidname = setext(fidname,'fid',0);

% create an empty header
hdr = cbiCreateNiftiHeader;

% get the qform
[qform44 info] = fid2xform(fidname);

% set the qform
hdr = cbiSetNiftiQform(hdr,qform44);

% now set dimensions in header
nDims = length(info.dim);
hdr.dim(1:nDims+1) = [nDims info.dim];

% set voxel dimensions 
hdr.pixdim(2:4) = info.voxsize(1:3);

% set the tr
hdr.pixdim(5) = info.tr*1000;



