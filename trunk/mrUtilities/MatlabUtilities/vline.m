% vline.m
%
%   usage: vline(hpos,linetype)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a vertical line on the current axis
%
function vline(hpos,linetype)

ax = axis;

if (nargin == 1)
  for i = 1:length(hpos)
    hold on
    plot([hpos(i) hpos(i)],[ax(3) ax(4)],'k:');
  end
elseif (nargin == 2)
  for i = 1:length(hpos)
    hold on
    plot([hpos(i) hpos(i)],[ax(3) ax(4)],linetype);
  end
else
  help vline
end

