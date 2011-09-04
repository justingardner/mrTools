% mlrImageGetAllFilenames.m
%
%        $Id:$ 
%      usage: imageFilenames = mlrImageGetAllFilenames(<dirname>,<'mustNotHaveDotInFilename=0'>)
%         by: justin gardner
%       date: 07/23/09
%    purpose: Function to return all the valid image filenames in a directory. This looks for all nifti
%             files and will only return the .hdr filename for dual files. It will also validate the 
%             headers and not return filenames with invalid nifti headers. It will return either hdr or nii files.
%
%             dirname is the name of the directory to operate on. defaults to current directory
%
function imageFilenames = mlrImageGetAllFilenames(dirname,varargin)

imageFilenames = {};

% check arguments
if ~any(nargin == [0 1 2 3 4])
  help mlrImageGetAllFilenames
  return
end

mustNotHaveDotInFilename = [];
getArgs(varargin,{'mustNotHaveDotInFilename=0'});

if ieNotDefined('dirname'),dirname = pwd;end

% get the direcory
d = dir(dirname);

% list of valid extensions
validExtensions = {'hdr','nii'};

for i = 1:length(d)
  % check for valid extension
  if any(strcmp(getext(d(i).name),validExtensions))
    % check for valid filename
    if mlrImageIsImage(fullfile(dirname,d(i).name))
      imageFilenames{end+1} = d(i).name;
    end
  end
end

% remove any filenames that have dots in the middle of them, since this causes weird problems later
% since you can't depend on the dot marking extensions
if mustNotHaveDotInFilename
  imageFilenamesWithoutDot = {};
  for i = 1:length(imageFilenames)
    if isempty(strfind(stripext(imageFilenames{i}),'.'))
      imageFilenamesWithoutDot{end+1} = imageFilenames{i};
    else
      mrWarnDlg(sprintf('(mlrImageGetAllFilenames) Ignoring file %s because it has a . in the filename that does not mark the file extension. If you want to use this file, consider renaming to %s',imageFilenames{i},setext(fixBadChars(stripext(imageFilenames{i}),{'.','_'}),'hdr')));
    end
  end
  imageFilenames = imageFilenamesWithoutDot;
end



