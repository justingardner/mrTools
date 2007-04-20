% fixBadChars.m
%
%      usage: str = fixBadChars(str)
%         by: justin gardner
%       date: 04/20/07
%    purpose: takes a string and replaces bad characters not
%    allowed in variable names like space or * with variable name acceptable characters
%
function str = fixBadChars(str)

% check arguments
if ~any(nargin == [1])
  help fixBadChars
  return
end

% this is the list of what characters will map to what
swapchars = {{'-','_'},{' ','_'},{'*','s'},{'+','p'},{'%','p'},{'[','B'},{']','B'},{'(','P'},{')','P'}};

% now swap any occurrences of these characters
for i = 1:length(swapchars)
  swaplocs = strfind(str,swapchars{i}{1});
  str(swaplocs) = swapchars{i}{2};
end

% check for non character at beginning
if ~isempty(regexp(str,'^[^a-zA-Z]'))
  str = sprintf('x%s',str);
end
