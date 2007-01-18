% xaxis
%
%      usage: xaxis()
%         by: justin gardner
%       date: 02/17/03
%    purpose: sets the x axis of a figure
%
function xaxis(xmin,xmax)

if ((nargin == 1) && (length(xmin) == 2))
  xmax = xmin(2);
  xmin = xmin(1);
elseif (nargin ~= 2)
  help xaxis;
  return
end

a = axis;
axis([xmin xmax a(3) a(4)]);
