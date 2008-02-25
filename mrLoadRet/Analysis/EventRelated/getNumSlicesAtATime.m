% getNumSlicesAtATime.m
%
%        $Id$
%      usage: getNumSlicesAtATime(numVolumes,dims)
%         by: justin gardner
%       date: 05/18/07
%    purpose: returns how many slices to process at a time,
%             dependent on the preference maxBlocksize
%             numVolumes is the number of volumes in the scan
%             dims is the dimensions of the scan
%            
%             The number of slices will always be an integer value
%             >=1. To get the unrounded value, check the second
%             return argument
%
%
function [numSlicesAtATime rawNumSlices] = getNumSlicesAtATime(numVolumes,dims,precision)

% check arguments
if ~any(nargin == [2 3])
  help getNumSlicesAtATime
  return
end

% get the precision (# of bytes)
if ieNotDefined('precision')
  precision = 'double';
end
switch (precision)
 case {'double'}
  precision = 8;
 case{'single'}
  precision = 4;
 otherwise
  precision = 8;
end

maxBlocksize = mrGetPref('maxBlocksize');
if ieNotDefined('maxBlocksize')
  maxBlocksize = 250000000;
end
rawNumSlices = maxBlocksize/(precision*numVolumes*prod(dims(1:2)));
numSlicesAtATime = max(1,floor(rawNumSlices));

