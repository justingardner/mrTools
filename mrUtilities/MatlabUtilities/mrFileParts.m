% this function replaces base Matlab's fileparts and 
% correctly handle the .nii.gz extension

function [path,filename,ext] = mrFileParts(fullPath)

[path,filename,ext] = fileparts(fullPath);

if strcmp(ext,'.gz')
  [~,filename, ext2] = fileparts(filename);
  if strcmp(ext2,'.nii')
    ext = '.nii.gz';
  else
    filename = [filename ext2];
  end
end
