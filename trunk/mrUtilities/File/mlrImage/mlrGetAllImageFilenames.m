% mlrGetAllImageFilenames.m
%
%        $Id:$ 
%      usage: imageFilenames = mlrGetAllImageFilenames(<dirname>)
%         by: justin gardner
%       date: 07/23/09
%    purpose: Function to return all the valid image filenames in a directory. This looks for all nifti
%             files and will only return the .hdr filename for dual files. It will also validate the 
%             headers and not return filenames with invalid nifti headers. It will return either hdr or nii files.
%
%             dirname is the name of the directory to operate on. defaults to current directory
%
function imageFilenames = mlrGetAllImageFilenames(dirname)

imageFilenames = {};

% check arguments
if ~any(nargin == [0 1])
  help mlrGetAllImageFilenames
  return
end

if ieNotDefined('dirname'),dirname = pwd;end

% get the direcory
d = dir(dirname);

% list of valid extensions
validExtensions = {'hdr','nii'};

for i = 1:length(d)
  % check for valid extension
  if any(strcmp(getext(d(i).name),validExtensions))
    % check for valid filename
    if mlrIsImageFile(fullfile(dirname,d(i).name))
      imageFilenames{end+1} = d(i).name;
    end
  end
end




