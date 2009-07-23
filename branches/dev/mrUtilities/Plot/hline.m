% hline.m
%
%   usage: hline(vpos,linetype)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a horizontal line on the current axis
%          returns handles 
%          vpos can be an array of positions
function h = hline(vpos,linetype)

ax = axis;
h = [];
if (nargin == 1)
  for i = 1:length(vpos)
    hold on
    h(i) = plot([ax(1) ax(2)],[vpos(i) vpos(i)],'k:');
  end
elseif (nargin == 2)
  for i = 1:length(vpos)
    hold on
    h(i) = plot([ax(1) ax(2)],[vpos(i) vpos(i)],linetype);
  end
else
  help hline
end

