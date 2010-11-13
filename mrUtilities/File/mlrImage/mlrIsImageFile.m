% mlrIsImageFile.m
%
%        $Id$ 
%      usage: retval = mlrIsImageFile(filename)
%         by: justin gardner
%       date: 07/23/09
%    purpose: Test to see if the file is a valid image, if the extension is provided, will
%             test that filename explicitly. If no extension is given, will append the
%             mrGetPref('niftiFileExtension') to the filename
% 
%             returns 0 if not a valid image file, 1 otherwise
%
function retval = mlrIsImageFile(filename)

retval = 0;

% check arguments
if ~any(nargin == [1])
  help mlrIsImageFile
  return
end

% put on nifti file extension, if there is no extension given
if isempty(getext(filename))
  filename = setext(filename,mrGetPref('niftiFileExtension'));
end

% check if the file exists
if ~isfile(filename)
  disp(sprintf('(mlrIsImageFile) File %s does not exist',filename));
  return
end

% no see if we can open the header
hdr = cbiReadNiftiHeader(filename);
if isempty(hdr),return,end

% if we got here, it must be an image file
retval = true;
