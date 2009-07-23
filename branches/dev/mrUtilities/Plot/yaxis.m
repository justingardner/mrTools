% yaxis
%
%      usage: yaxis()
%         by: justin gardner
%       date: 02/17/03
%    purpose: sets the y axis of a figure
%
function yaxis(ymin,ymax)

if ((nargin == 1) && (length(ymin) == 2))
  ymax = ymin(2);
  ymin = ymin(1);
elseif (nargin ~= 2)
  help yaxis;
  return
end

a = axis;

% default to what is already there
if isempty(ymin),ymin=a(3);end
if isempty(ymax),ymax=a(4);end

axis([a(1) a(2) ymin ymax]);
