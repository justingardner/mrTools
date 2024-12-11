% mrSub2ind.m
%
%        $Id$	
%      usage: mrSub2ind(dims,x,y,z)
%         by: justin gardner
%       date: 10/15/07
%    purpose: sub2ind that doesn't choke on coordinates outside of dims.
%
function linear = mrSub2ind(dims,x,y,z)

% check arguments
if ~any(nargin == [3 4])
  help mrSub2ind
  return
end

if isempty(x),linear = [];return;end

if nargin == 4
  if verLessThan('matlab','24.1')
    badCoords = find((x < 1) | (x > dims(1)) | ...
		   (y < 1) | (y > dims(2)) | ...
		   (z < 1) | (z > dims(3)));
    x(badCoords) = nan;
    y(badCoords) = nan;
    z(badCoords) = nan;
    
    linear = sub2ind(dims,x,y,z);
  else   % As of Matlab 2024, nan coordinates in sub2ind give an error, these coordinates need to be removed completely
    badCoords = (x < 1) | (x > dims(1)) | isnan(x) | ...
		   (y < 1) | (y > dims(2)) | isnan(y) | ...
		   (z < 1) | (z > dims(3)) | isnan(z);
    linear = nan(size(x));
    linear(~badCoords) = sub2ind(dims,x(~badCoords),y(~badCoords),z(~badCoords));
  end
else
  if verLessThan('matlab','24.1')
    badCoords = find((x < 1) | (x > dims(1)) | ...
		     (y < 1) | (y > dims(2)));
    x(badCoords) = nan;
    y(badCoords) = nan;
    
    linear = sub2ind(dims,x,y);
  else   % As of Matlab 2024, nan coordinates in sub2ind give an error, these coordinates need to be removed completely
    badCoords = (x < 1) | (x > dims(1)) | isnan(x) | ...
		   (y < 1) | (y > dims(2)) | isnan(y);
    linear = nan(size(x));
    linear(~badCoords) = sub2ind(dims,x(~badCoords),y(~badCoords));
  end
end

