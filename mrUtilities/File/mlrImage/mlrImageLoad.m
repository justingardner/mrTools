% mlrImageLoad.m
%
%        $Id:$ 
%      usage: mlrImageLoad(filename)
%         by: justin gardner
%       date: 08/15/11
%    purpose: Loads an mlr image (this is usually a nifti file, but we can expand this
%             to load any kind of image).
%
%             To load an image based on a scan/group
%             v = newView;
%             mlrImageLoad(v,'groupNum=2','scanNum=3');
%
function [data header] = mlrImageLoad(filename,varargin)

% default return values
data = [];header = [];

% check arguments
if nargin < 1
  help mlrImageLoad
  return
end

% set filenames
if isstr(filename) && isempty(getext(filename))
  filename = setext(filename,mrGetPref('niftiFileExtension'));
end

% load the header first
header = mlrImageHeaderLoad(filename);
if isempty(header),return,end

% if this is a data structure then extract data
if isstruct(filename) && isfield(filename,'data')
  filename = filename.data;
end
if isnumeric(filename)
  data = filename;
  return
end

% check input arguments
groupNum = [];scanNum = 1;
getArgs(varargin,{'groupNum=1','scanNum=1'});

% if the passed in filename is a view, then load the appropriate group and scan
if isview(filename)
  v = filename;
  filename = viewGet(v,'tseriespathstr',scanNum,groupNum);
end

% load the data
switch lower(getext(filename))
 case {'hdr','img','nii'}
  data = cbiReadNifti(filename);
 case {'sdt','spr','edt','epr'}
  data = readsdt(filename);
  data = data.data;
 case {'fid'}
  data = fid2nifti(filename);
end

