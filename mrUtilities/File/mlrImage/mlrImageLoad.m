% mlrImageLoad.m
%
%        $Id:$ 
%      usage: mlrImageLoad(filename)
%         by: justin gardner
%       date: 08/15/11
%    purpose: Loads an mlr image (handles nifti, Agilant fid, and
%             other formats)
%  
%             You can load using a filename, view/scanNum/groupNum,
%             select from a dialog box in the current or canonical
%             directory, or from a struct -> see mlrImageParseArgs
%             for details
%
%             For details on the header see:  mlrImageIsHeader
%
%             To reorder the image matrix (and qform/sform) to
%             be in standard LPI orientation, pass arg:
%             'orient=LPI'
%
function [dataRetval headerRetval] = mlrImageLoad(varargin)

% default return values
dataRetval = [];headerRetval = [];

% check arguments
if nargin < 1
  help mlrImageLoad
  return
end

% parse arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);

% check input arguments
verbose=0;orient=[];xMin=1;xMax=inf;yMin=1;yMax=inf;zMin=1;zMax=inf;
validArgs = {'verbose','orient','xMin','xMax','yMin','yMax','zMin','zMax'};
getArgs(otherArgs,validArgs);

% number of images to load. Note that for
% a single image, then we just return the data and header.
% for multiple images, we will return a cell array
% of headers
nImages = length(imageArgs);

for iImage = 1:nImages
  % get the current filename
  filename = imageArgs{iImage};
  data = [];header = [];
  
  % load the header first
  header = mlrImageHeaderLoad(filename);
  if isempty(header)
    if nImages > 1
      headerRetval{iImage} = [];
      dataRetval{iImage} = [];
    else
      headerRetval = [];
      dataRetval = [];
    end
    continue;
  end

  % if this is a struct with filename and loadArgs field
  % then we split those out
  loadArgs = {};altArgs = {};
  if isfield(filename,'filename')
    if isfield(filename,'loadArgs')
      loadArgs = filename.loadArgs;
    end
    if isfield(filename,'altArgs')
      altArgs = filename.altArgs;
    end
    filename = filename.filename;
  end
  
  % set alt args if there are any
  getArgs(altArgs,validArgs);
  
  % ok, now load data
  % if this is a structure with the data field then just return that
  if isstruct(filename) && isfield(filename,'data')
    data = filename.data;
    if isempty(data)
      disp(sprintf('(mlrImageLoad) Empty image'));
      return
    end
  elseif isstr(filename)
    % load the data
    switch lower(getext(filename))
     case {'img'}
      if isdir(filename)
	data = fdf2nifti(filename,verbose);
      else
	data = cbiReadNifti(filename);
      end
     case {'hdr','nii'}
      data = cbiReadNifti(filename);
     case {'sdt','spr','edt','epr'}
      data = readsdt(filename);
      if isfield(data,'data')
	data = data.data;
      else
	data = [];
      end
     case {'fid'}
      data = fid2nifti(filename,verbose,'loadArgs',loadArgs);
    end
  end

  % now make sure dimensions match in header
  header.dim = size(data)';
  header.nDim = length(header.dim);

  % make sure the header is correct
  [tf header] = mlrImageIsHeader(header);

  % shift position if called for
  [data header] = adjustDims(data,header,xMin,xMax,yMin,yMax,zMin,zMax);

  % fix orientation if called for
  if ~isempty(orient)
    [data header] = mlrImageOrient(orient,data,header);
  end

  % package up for returning
  if nImages > 1
    headerRetval{iImage} = header;
    dataRetval{iImage} = data;
  else
    headerRetval = header;
    dataRetval = data;
  end
end

%%%%%%%%%%%%%%%%%%%%
%%   adjustDims   %%
%%%%%%%%%%%%%%%%%%%%
function [data h] = adjustDims(data,h,xMin,xMax,yMin,yMax,zMin,zMax)

% if nothing to do, just return
if (xMin == 1) && (xMax == inf) && (yMin == 1) && (yMax == inf) && (zMin == 1) && (zMax == inf)
  return
end

% get current dimensions
dims = size(data);

% make the dimensions valid
xMin = round(max(1,xMin));
xMax = round(min(xMax,dims(1)));
yMin = round(max(1,yMin));
yMax = round(min(yMax,dims(2)));
zMin = round(max(1,zMin));
zMax = round(min(zMax,dims(3)));

% adjust data size
data = data(xMin:xMax,yMin:yMax,zMin:zMax,:);

% find out how much the new xMin, yMin, zMin have translated the image
t = h.qform * [xMin yMin zMin 1]' - h.qform * [1 1 1 1]';
t = [zeros(3,3) t(1:3,1); 0 0 0 0];
h.qform = t + h.qform;

% reset the dims
h.dim = size(data);



