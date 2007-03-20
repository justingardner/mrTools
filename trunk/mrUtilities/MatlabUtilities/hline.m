% hline.m
%
%   usage: hline(hpos,linetype)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a horizontal line on the current axis
%
function retval = hline(vpos,linetype)

ax = axis;
retval = [];
if (nargin == 1)
  for i = 1:length(vpos)
    hold on
    retval = plot([ax(1) ax(2)],[vpos(i) vpos(i)],'k:');
  end
elseif (nargin == 2)
  for i = 1:length(vpos)
    hold on
    retval = plot([ax(1) ax(2)],[vpos(i) vpos(i)],linetype);
  end
else
  help hline
end

