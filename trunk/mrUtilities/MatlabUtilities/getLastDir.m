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

[pathStr lastDir] = fileparts(pathStr);


