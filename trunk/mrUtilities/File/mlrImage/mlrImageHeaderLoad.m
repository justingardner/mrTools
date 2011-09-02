% mlrImageHeaderLoad.m
%
%        $Id: mlrImageHeaderLoad.m 1268 2008-08-19 15:43:36Z justin $ 
%      usage: mlrImageHeaderLoad(filename)
%         by: justin gardner
%       date: 08/19/08
%    purpose: loads a mrLoadRet image. This can handle various image types
%             including nifti. 
%
%             To load an image header based on a scan/group
%             v = newView;
%             mlrImageHeaderLoad(v,'groupNum=2','scanNum=3');
%
%             To select a canonical (through dialog box) from the volumeDirecotry
%             mlrImageHeaderLoad canonical
function header = mlrImageHeaderLoad(filename,varargin)

header = [];

% check arguments
if nargin < 1
  help mlrLoadHeader
  return
end

% if passed in a sturcutre with the field data
% then grab the data and convert the rest of the
% fields to a private header
if (isstruct(filename) && isfield(filename,'data'))
  header.hdr = rmfield(filename,'data');
  filename = filename.data;
end

% numeric argument means passed in data structure with no header
if isnumeric(filename)
  header.dim = size(filename);
  [tf header] = mlrImageIsHeader(header);
  return
end

% check input arguments
groupNum = [];scanNum = [];verbose = [];
getArgs(varargin,{'groupNum=1','scanNum=1','verbose=0'});

% if the passed in filename is a view, then load the appropriate group and scan
if isview(filename)
  v = filename;
  filename = viewGet(v,'tseriespathstr',scanNum,groupNum);
end


% set filenames
if isempty(getext(filename))
  hdrFilename = setext(filename,mrGetPref('niftiFileExtension'));
  % get from canonical directory
  if any(strcmp({'canonical','volume','volumedirectory','volumedir','voldir'},lower(stripext(filename)))) && ~isfile(filename)
    filename = getPathStrDialog(mrGetPref('volumeDirectory'),'Choose a volume',{'*.hdr;*.nii', 'Nifti Files (*.hdr, *.nii)'},'off');
    if isempty(filename),return,end
  end
else
  hdrFilename = filename;
end

%check for file
if ~isfile(hdrFilename) && ~isdir(hdrFilename)
  disp(sprintf('(mlrImageHeaderLoad) Could not find file %s',hdrFilename));
  return
end

% set inital fields of header
header.filename = filename;
header.ext = getext(filename);

switch lower(getext(filename))
  case {'img'}
    if ~isdir(filename)
      header = mlrImageHeaderLoadNifti(filename,header);
    else
      header = mlrImageHeaderLoadFDF(filename,header,verbose);
    end
  case {'hdr','nii'}
    header = mlrImageHeaderLoadNifti(filename,header);
  case {'sdt','spr','edt','epr'}
    header = mlrImageHeaderLoadSDT(filename,header);
  case {'fid'}
    header = mlrImageHeaderLoadFid(filename,header,verbose);
 otherwise
    disp(sprintf('(mlrImageHeaderLoad) Unknown image header type: %s',getext(filename)));
    return
end

if isempty(header),return,end

% load the associated matlab header
header.base = loadAssociatedMatlabHeader(filename);

% set the vol2mag field - get it from the base 
% struture if it has not already been set
if ~isfield(header,'vol2mag')
  if isfield(header.base,'vol2mag')
    header.vol2mag = header.base.vol2mag;
  else
    header.vol2mag = [];
  end
end
% set the vol2tal field - get it from the base 
% struture if it has not already been set
if ~isfield(header,'vol2tal')
  if isfield(header.base,'vol2tal')
    header.vol2tal = header.base.vol2tal;
  else
    header.vol2tal = [];
  end
end

% now fill in any missing fields from header
[tf header] = mlrImageIsHeader(header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadFDF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadFDF(filename,header,verbose)

% read the nifti header
[d nifti] = fdf2nifti(filename,verbose,true);
if isempty(nifti),header = [];return;end
  
% set some info
header.type = 'nifti';
header.hdr = nifti;

% set other fields
header.qform_code = header.hdr.qform_code;
header.qform44 = header.hdr.qform44;
header.sform_code = header.hdr.sform_code;
header.sform44 = header.hdr.sform44;
header.nDim = header.hdr.dim(1);
header.dim = header.hdr.dim(2:header.nDim+1);
header.pixdim = header.hdr.pixdim(2:header.nDim+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadFid    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadFid(filename,header,verbose)

% read the nifti header
nifti = fid2niftihdr(filename,verbose);
if isempty(nifti),header = [];return;end
  
% set some info
header.type = 'nifti';
header.hdr = nifti;

% set other fields
header.qform_code = header.hdr.qform_code;
header.qform44 = header.hdr.qform44;
header.sform_code = header.hdr.sform_code;
header.sform44 = header.hdr.sform44;
header.nDim = header.hdr.dim(1);
header.dim = header.hdr.dim(2:header.nDim+1);
header.pixdim = header.hdr.pixdim(2:header.nDim+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadSDT    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadSDT(filename,header)

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
%    mlrImageHeaderLoadNifti    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadNifti(filename,header)

% read the nifti header
nifti = cbiReadNiftiHeader(filename);

% set some info
header.type = 'nifti';
header.hdr = nifti;

% set other fields
header.qform_code = header.hdr.qform_code;
header.qform44 = header.hdr.qform44;
header.sform_code = header.hdr.sform_code;
header.sform44 = header.hdr.sform44;
header.nDim = header.hdr.dim(1);
header.dim = header.hdr.dim(2:header.nDim+1);
header.pixdim = header.hdr.pixdim(2:header.nDim+2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadAssociatedMatlabHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function base = loadAssociatedMatlabHeader(filename)

base = [];
matlabFilename = setext(filename,'mat');
if isfile(matlabFilename)
  % load the header
  matHeader = load(matlabFilename);
  % check for field
  if ~isfield(matHeader,'base')
    disp(sprintf('(mlrImageHeaderLoad) Ignoring associated mat file %s because is not a MLR header',filename));
  else
    matHeader.base.data = [];
    [tf header] = isbase(matHeader.base);
    if ~tf
      disp(sprintf('(mlrImageHeaderLoad) Ignoring associated mat file %s because is not a MLR header',filename));
    else
      base = matHeader.base;
    end
  end
end
