% mlrLoadImage.m
%
%        $Id:$ 
%      usage: mlrLoadImage(filename)
%         by: justin gardner
%       date: 08/15/11
%    purpose: Loads an mlr image (this is usually a nifti file, but we can expand this
%             to load any kind of image).
%
function [data header] = mlrLoadImage(filename)

% check arguments
if ~any(nargin == [1])
  help mlrLoadImage
  return
end

% set filenames
if isempty(getext(filename))
  filename = setext(filename,mrGetPref('niftiFileExtension'));
end

% load the header first
header = mlrLoadImageHeader(filename);

% load the data
switch lower(getext(filename))
 case {'hdr','img','nii'}
  data = cbiReadNifti(filename);
 case {'sdt','spr','edt','epr'}
  data = readsdt(filename);
  data = data.data;
end

