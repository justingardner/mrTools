% stripext.m
%
%      usage: stripext(filename)
%         by: justin gardner
%       date: 02/08/05
%    purpose: remove extension if it exists
%
function retval = stripext(filename,delimiter)

if ~any(nargin == [1 2])
  help stripext;
  return
end
% dot delimits end
if exist('delimiter')~=1,delimiter='.';,end

retval = filename;
dotloc = findstr(filename,delimiter);
if length(dotloc) > 0
  retval = filename(1:dotloc(length(dotloc))-1);
end
