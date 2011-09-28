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
verbose=0;orient=[];xMin=1;xMax=inf;yMin=1;yMax=inf;zMin=1;zMax=inf;volNum = [];
returnRaw=0;
swapXY=0;swapXZ=0;swapYZ=0;flipX=0;flipY=0;flipZ=0;shiftX=0;shiftY=0;shiftZ=0;
rotateXY=0;rotateXZ=0;rotateYZ=0;interpMethod='linear';applyToHeader=1;applyToData=1;
validArgs = {'verbose','orient','xMin','xMax','yMin','yMax','zMin','zMax','volNum','swapXY','swapXZ','swapYZ','flipX','flipY','flipZ','shiftX','shiftY','shiftZ','rotateXY','rotateXZ','rotateYZ','interpMethod','applyToHeader','applyToData','returnRaw'};
getArgs(otherArgs,validArgs);

% check volNum argument
if ~isempty(volNum) && (length(volNum)>2)
  disp(sprintf('(mlrImageLoad) volNum must be either a single volume or a [min max] array. Ignoring.'));
  volNum = [];
else
  volNum = sort(volNum);
end

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

  % setup for volNum processing
  if ~isempty(volNum) && (header.nDim >= 4)
    volNum(volNum < 1) = 1;
    volNum(volNum > header.dim(4)) = header.dim(4);
  else
    volNum = [];
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
	if isempty(volNum)
	  data = cbiReadNifti(filename);
	else
	  data = cbiReadNifti(filename,{[],[],[],volNum});
	end
      end
      volNum = [];
     case {'hdr','nii'}
      if isempty(volNum)
	data = cbiReadNifti(filename);
      else
	data = cbiReadNifti(filename,{[],[],[],volNum});
      end
      volNum = [];
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

  % if we are allowing xformations
  if ~returnRaw
    % fix orientation if called for
    if ~isempty(orient)
      [data header] = mlrImageOrient(orient,data,header);
    end

    % do any shift or xforms on header
    [data header] = mlrImageXform(data,header,'swapXY',swapXY,'swapXZ',swapXZ,'swapYZ',swapYZ,'flipX',flipX,'flipY',flipY,'flipZ',flipZ,'shiftX',shiftX,'shiftY',shiftY,'shiftZ',shiftZ,'rotateXY',rotateXY,'rotateXZ',rotateXZ,'rotateYZ',rotateYZ,'xMin',xMin,'xMax',xMax,'yMin',yMin,'yMax',yMax,'zMin',zMin,'zMax',zMax,'applyToData',applyToData,'applyToHeader',applyToHeader,'interpMethod',interpMethod,'verbose',verbose);
  end

  % do volNum processing if not already done
  if ~isempty(volNum) && (header.nDim >= 4)
    volNum(volNum < 1) = 1;
    volNum(volNum > header.dim(4)) = header.dim(4);
    if length(volNum) == 1
      data = data(:,:,:,volNum,:);
    else
      data = data(:,:,:,volNum(1):volNum(2),:);
    end
    header.dim(4) = size(data,4);
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

