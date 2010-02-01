% vline.m
%
%   usage: vline(hpos,linetype)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a vertical line on the current axis
%          returns handles 
%          hpos can be an array of positions
function h = vline(hpos,linetype)

ax = axis;
miny = ax(3);maxy = ax(4);
if isequal(get(gca,'YScale'),'log')
  miny = min(get(gca,'Ytick'));
end
h = [];
if (nargin == 1)
  for i = 1:length(hpos)
    hold on
    h(i) = plot([hpos(i) hpos(i)],[miny maxy],'k:');
  end
elseif (nargin == 2)
  for i = 1:length(hpos)
    hold on
    h(i) = plot([hpos(i) hpos(i)],[miny maxy],linetype);
  end
else
  help vline
end

