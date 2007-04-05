% cellArray.m
%
%      usage: var = cellArray(var)
%         by: justin gardner
%       date: 04/05/07
%    purpose: when passed a single structure it returns it as a
%    cell array of length 1, if var is already a cell array just
%    passes it back
%
function var = cellArray(var)

% check arguments
if ~any(nargin == [1])
  help cellArray
  return
end

if ~iscell(var)
  tmp = var;
  clear var;
  var{1} = tmp;
end

