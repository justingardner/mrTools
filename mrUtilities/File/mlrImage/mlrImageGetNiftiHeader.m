% mlrImageGetNiftiHeader.m
%
%        $Id:$ 
%      usage: hdr = mlrImageGetNiftiHeader(h)
%         by: justin gardner
%       date: 11/15/11
%    purpose: Pass in an mlrImage header and will return a nifti header
%             To create a new header:
%             h.dim = [64 64 3];
%             hdr = mlrImageGetNiftiHeader(h);
%
function hdr = mlrImageGetNiftiHeader(h)

hdr = [];

% check arguments
if ~any(nargin == [1])
  help mlrImageGetNiftiHeader
  return
end

% check to make sure the header is an mlrImage
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageGetNiftiHeader) Passed in header is not an mlrImage header'));
  return
end

% create the nifti header
hdr = cbiCreateNiftiHeader;
hdr.dim = [h.nDim h.dim];
hdr.dim(end+1:8) = 0;
hdr.dim = hdr.dim(:);
hdr.dim = hdr.dim(1:8);
hdr.pixdim = [h.nDim h.pixdim];
hdr.pixdim(end+1:8) = 0;
hdr.pixdim = hdr.pixdim(:);
hdr.pixdim = hdr.pixdim(1:8);

% set the qform
if ~isempty(h.qform)
  hdr = cbiSetNiftiQform(hdr,h.qform);
  hdr.qform_code = 1;
  % if the header came from a nifti file and
  % there happens to be a qform_code that is
  % greater than one than use that (this means
  % that this a talairach or some other xform)
  % this info is not supported mlrImage, but
  % to keep compatibility we can set it
  if isfield(h,'hdr') && isfield(h.hdr,'qform_code') && (h.hdr.qform_code > 1)
    hdr.qform_code = h.hdr.qform_code;
  end
end

% set the sform
if ~isempty(h.sform)
  hdr = cbiSetNiftiSform(hdr,h.sform);
  hdr.sform_code = 1;
  % same as sform_code
  if isfield(h,'hdr') && isfield(h.hdr,'sform_code') && (h.hdr.sform_code > 1)
    hdr.sform_code = h.hdr.sform_code;
  end
end

