% getNumSlicesAtATime.m
%
%        $Id$
%      usage: getNumSlicesAtATime(numVolumes,dims)
%         by: justin gardner
%       date: 05/18/07
%    purpose: returns how many slices to process at a tiem,
%             dependent on the preference maxBlocksize
%             numVolumes is the number of volumes in the scan
%             dims is the dimensions of the scan
%
function numSlicesAtATime = getNumSlicesAtATime(numVolumes,dims)

% check arguments
if ~any(nargin == [2])
  help getNumSlicesAtATime
  return
end

maxBlocksize = mrGetPref('maxBlocksize');
if ieNotDefined('maxBlocksize')
  maxBlocksize = 250000000;
end
numSlicesAtATime = max(1,floor(maxBlocksize/(8*numVolumes*prod(dims(1:2)))));
