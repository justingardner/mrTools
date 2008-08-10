% setext.m
%
%        $Id$ 
%      usage: filename = setext(filename,ext)
%         by: justin gardner
%       date: 08/09/08
%    purpose: Makes sure that file the filename has the specified extension
%
function filename = setext(filename,ext)

% check arguments
if ~any(nargin == [2])
  help setext
  return
end

% strip any leading '.'
if length(ext) && isequal(ext(1),'.')
  ext = ext(2:end);
end

% set the extension
filename = sprintf('%s.%s',stripext(filename),ext);

