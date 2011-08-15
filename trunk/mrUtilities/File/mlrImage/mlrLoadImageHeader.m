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

%check for file
if ~isfile(hdrFilename) 
  disp(sprintf('(mlrLoadImageHeader) Could not find file %s',hdrFilename));
  return
end

% set inital fields of header
header.filename = filename;
header.ext = getext(filename);

% set defaults
header.qform_code = 0;
header.qform44 = eye(4);
header.sform_code = 0;
header.sform44 = eye(4);
header.permuationMatrix = eye(4);
header.nDim = 0;
header.dim = [];
header.pixdim = [];

switch lower(getext(filename))
  case {'hdr','nii','img'}
    header = mlrLoadImageHeaderNifti(filename,header);
  case {'sdt','spr','edt','epr'}
    header = mlrLoadImageHeaderSDT(filename,header);
 otherwise
    disp(sprintf('(mlrLoadImageHeader) Unknown image header type: %s',getext(filename)));
    return
end

% load the associated matlab header
header.matHeader = loadAssociatedMatlabHeader(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLoadImageHeaderSDT    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrLoadImageHeaderSDT(filename,header)

% load sdt
d = readsdt(filename,0);
if isempty(d),return,end

% set fileds
if any(strcmp(getext(filename),{'sdt','spr'}))
  header.type = 'sdt';
else
  header.type = 'edt';
end
header.hdr = d;

% set fields
header.nDim = length(d.dim);
header.dim = d.dim;
if isfield(d,'interval')
  header.pixdim  = d.interval(1:3)*10;
elseif isfield(d,'fov')
  header.pixdim = d.fov(1:3)./d.dim(1:3)*10
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLoadImageHeaderNifti    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrLoadImageHeaderNifti(filename,header)

% read the nifti header
nifti = cbiReadNiftiHeader(hdrFilename);

% set some info
header.type = 'nifti';
header.hdr = nifti;

% set other fields
header.qform_code = header.hdr.qform_code;
header.qform44 = header.hdr.qform44;
header.sform_code = header.hdr.sform_code;
header.sform44 = header.hdr.sform44;
header.pixdim = header.hdr.pixdim;
header.dim = header.hdr.dim;
header.permutationMatrix = getPermutationMatrix(nifti);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadAssociatedMatlabHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matHeader = loadAssociatedMatlabHeader(filename)

matHeader = [];
matlabFilename = setext(filename,'mat');
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
