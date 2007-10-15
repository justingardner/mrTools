% mrSub2ind.m
%
%      usage: mrSub2ind()
%         by: justin gardner
%       date: 10/15/07
%    purpose: sub2ind that doesn't choke on coordinates
%             outside of dims (just inserts nans for those)
%
function linear = mrSub2ind(dims,x,y,z)

% check arguments
if ~any(nargin == [3 4])
  help mrSub2ind
  return
end

if nargin == 4
  badCoords = find((x < 1) | (x > dims(1)) | ...
		   (y < 1) | (y > dims(2)) | ...
		   (z < 1) | (z > dims(3)));
  x(badCoords) = nan;
  y(badCoords) = nan;
  z(badCoords) = nan;

  linear = sub2ind(dims,x,y,z);
else
  badCoords = find((x < 1) | (x > dims(1)) | ...
		   (y < 1) | (y > dims(2)));
  x(badCoords) = nan;
  y(badCoords) = nan;

  linear = sub2ind(dims,x,y);
end

