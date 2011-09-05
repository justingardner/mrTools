% mlrImageHeaderLoad.m
%
%        $Id: mlrImageHeaderLoad.m 1268 2008-08-19 15:43:36Z justin $ 
%      usage: mlrImageHeaderLoad(filename)
%         by: justin gardner
%       date: 08/19/08
%    purpose: loads a mlr image. This can handle various image types
%             including nifti. 
%
%             You can load using a filename, view/scanNum/groupNum,
%             select from a dialog box in the current or canonical
%             directory, or from a struct -> see mlrImageParseArgs
%             for details
%
function retval = mlrImageHeaderLoad(varargin)

header = [];

% check arguments
if nargin < 1
  help mlrLoadHeader
  return
end

% parse arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
verbose = [];
getArgs(otherArgs,{'verbose=0'});

% number of images headers to load. Note that for
% a single image, then we just return the header
% for multiple images, we will return a cell array
% of headers
nImages = length(imageArgs);allHeaders = {};

% cycle though each image argument
for iImage = 1:nImages
  header = [];
  % if passed in a sturcutre with the field data
  % then grab the data and convert the rest of the
  % fields to a private header
  if isstruct(imageArgs{iImage})
    if isfield(imageArgs{iImage},'h')
      % see if it has an mlrImage header
      if mlrImageIsHeader(imageArgs{iImage}.h)
	header = imageArgs{iImage}.h;
      else
	% see if this is a nifti header
	header = mlrImageHeaderLoadPassedInNiftiHeader(imageArgs{iImage}.h);
      end
    end
    % if we didn't get it from above, then build a header
    % based on the data field
    if isempty(header) && isfield(imageArgs{iImage},'data')
      header.type = 'data';
      header.ext = '';
      header.hdr = rmfield(imageArgs{iImage},'data');
      header.dim = size(imageArgs{iImage}.data);
      [tf header] = mlrImageIsHeader(header);
    end
    % set in all headers if necessary
    if nImages > 1
      [tf allHeaders{end+1}] = mlrImageIsHeader(header);
    end
    continue;
  elseif isstr(imageArgs{iImage})
    filename = imageArgs{iImage};
    %check for file
    if ~isfile(filename) && ~isdir(filename)
      disp(sprintf('(mlrImageHeaderLoad) Could not find file %s',filename));
      if nImages > 1,allHeaders{end+1} = [];end
      continue
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

    if isempty(header)
      if nImages > 1,allHeaders{end+1} = [];end
      continue
    end

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
  else
    disp(sprintf('(mlrImageHeaderLoad) Unknown argument'));
  end

  % keep all headers for returning
  if nImages > 1,[tf allHeaders{end+1}] = mlrImageIsHeader(header);end
end

% return a cell array if multiple images headers are to be loaded
% or just the one header otherwise
if nImages > 1
  retval = allHeaders;
else
  [tf retval] = mlrImageIsHeader(header);
end
  
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadFDF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadFDF(filename,header,verbose)

% read the nifti header
[d nifti] = fdf2nifti(filename,verbose,true);
if isempty(nifti),header = [];return;end
  
% set some info
header.type = 'fdf';
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
header.type = 'fid';
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

% check for orientation info
if isfield(d,'qform_code')
  header.qform_code = d.qform_code;
end
if isfield(d,'qform44')
  header.qform44 = reshape(d.qform44,4,4);
end
if isfield(d,'sform_code')
  header.sform_code = d.sform_code;
end
if isfield(d,'sform44')
  header.sform44 = reshape(d.sform44,4,4);
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrImageHeaderLoadPassedInNiftiHeader   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadPassedInNiftiHeader(h)

header = [];
% first check some fields
checkFields = {'qform_code','qform44','sform_code','sform44','dim','pixdim'};
foundField = false;
for iField = 1:length(checkFields)
  if isfield(h,checkFields{iField})
    header.(checkFields{iField}) = h.(checkFields{iField});
    foundField = true;
  end
end
% if we found any of the fileds, then just call it a nifti header
if foundField
  header.type = 'nifti';
  header.hdr = h;
end

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
