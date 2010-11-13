% mlrLoadImageHeader.m
%
%        $Id: mlrLoadImageHeader.m 1268 2008-08-19 15:43:36Z justin $ 
%      usage: mlrLoadImageHeader(filename)
%         by: justin gardner
%       date: 08/19/08
%    purpose: loads a mrLoadRet image (this is a nifti/matlab pair)
%
function header = mlrLoadImageHeader(filename)

header = [];

% check arguments
if ~any(nargin == [1])
  help mlrLoadHeader
  return
end

% set filenames
if isempty(getext(filename))
  hdrFilename = setext(filename,mrGetPref('niftiFileExtension'));
else
  hdrFilename = filename;
end
matlabFilename = setext(filename,'mat');

%check for file
if ~isfile(hdrFilename) 
  disp(sprintf('(mlrLoadImageHeader) Could not find file %s',hdrFilename));
  return
end

% read the nifti header
nifti = cbiReadNiftiHeader(hdrFilename);

% check for matlab file
header = [];
if isfile(matlabFilename)
  % load the header
  matHeader = load(matlabFilename);
  % check for field
  if ~isfield(matHeader,'base')
    disp(sprintf('(mlrLoadImageHeader) Ignoring associated mat file %s because is not a MLR header',matHeader));
  else
    matHeader.base.data = [];
    [tf header] = isbase(matHeader.base);
    if ~tf
      disp(sprintf('(mlrLoadImageHeader) Ignoring associated mat file %s because is not a MLR header',matHeader));
    end
  end
end

% no already created header, so make a new one
if isempty(header)
  header.data = [];
  header.hdr = nifti;
  header.name = filename;
  header.permutationMatrix = getPermutationMatrix(nifti);
  [tf header] = isbase(header);
  if ~tf
    disp(sprintf('(mlrLoadImageHeader) Could not init header'));
    return
  end
end

% set the nifti field
header.hdr = nifti;

% set other fields
header.qform_code = header.hdr.qform_code;
header.qform44 = header.hdr.qform44;
header.sform_code = header.hdr.sform_code;
header.sform44 = header.hdr.sform44;
header.pixdim = header.hdr.pixdim;
header.dim = header.hdr.dim;
    




