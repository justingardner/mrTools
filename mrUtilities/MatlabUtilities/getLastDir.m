% getLastDir(pathStr)
%
%      usage: getLastDir(pathStr)
%         by: justin gardner
%       date: 04/05/07
%    purpose: 
%
function lastDir = getLastDir(pathStr)

% check arguments
if ~any(nargin == [1])
  help getLastDir
  return
end

% remove trailling fileseparator if it is there
if length(pathStr) && (pathStr(end) == filesep)
  pathStr = pathStr(1:end-1);
end

% and return last dir
[pathStr lastDir ext] = fileparts(pathStr);

% paste back on extension
lastDir = [lastDir ext];

