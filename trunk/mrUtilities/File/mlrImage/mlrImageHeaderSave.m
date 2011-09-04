% mlrImageHeaderSave.m
%
%      usage: mlrImageHeaderSave(filename,header)
%         by: justin gardner
%       date: 09/04/11
%    purpose: Saves the image header using filename. This
%             saves as a nifti header with a mat associated
%             file if there is info in the base
%
function [tf hdr] = mlrImageHeaderSave(filename,h)

% default return argument
tf = false;
hdr = [];

% check arguments
if ~any(nargin == [1 2])
  help mlrImageHeaderSave
  return
end

% set empty filename to filename in h
if (nargin < 2) || isempty(filename)
  filename = h.filename;
end

% check ext
ext = getext(filename);
if isempty(ext)
  filename = setext(filename,mrGetPref('niftiFileExtension'));
  ext = getext(filename);
end

if ~any(strcmp({'hdr','img','nii'},ext))
  disp(sprintf('(mlrImageHeaderSave) Cannot save file of ext %s. Try saving with a nifti extension (hdr,img or nii)',ext));
  return
end

% make sure the header is valid
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageHeaderSave) Invalid header'));
  return
end

% create the nifti h
hdr = cbiCreateNiftiHeader;
hdr.dim = [h.nDim h.dim];
hdr.dim(end+1:8) = 0;
hdr.dim = hdr.dim(:);
hdr.pixdim = [h.nDim h.pixdim];
hdr.pixdim(end+1:8) = 0;
hdr.pixdim = hdr.pixdim(:);

% set the qform
if h.qform_code
  hdr = cbiSetNiftiQform(hdr,h.qform44);
  hdr.qform_code = h.qform_code;
end

% set the sform
if h.sform_code
  hdr = cbiSetNiftiSform(hdr,h.sform44);
  hdr.sform_code = h.sform_code;
end

% if we are not asked to pass back the header
% then save it
if nargout < 2
  hdr = cbiWriteNiftiHeader(hdr,filename);
end

% see if we need to save out a matlab extension
if ~isempty(h.base)
  % set fields that might have changed
  base = h.base;
  base.vol2mag = h.vol2mag;
  base.vol2tal = h.vol2tal;
  % and save
  matFilename = setext(filename,'mat');
  save(matFilename,'base');
end

% success
tf = true;
