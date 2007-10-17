% color2RGB.m
%
%      usage: color = color2RGB(color)
%         by: justin gardner
%       date: 10/17/07
%    purpose: Converts a color name to an RGB value
%
function color = color2RGB(color)

% check arguments
if ~any(nargin == [1])
  help color2RGB
  return
end

if isstr(color)
  switch (color)
   case {'yellow','y'}, color = [1 1 0];
   case {'magenta','m'}, color = [1 0 1];
   case {'cyan','c'}, color = [0 1 1];
   case {'red','r'}, color = [1 0 0];
   case {'green','g'}, color = [0 1 0];
   case {'blue','b'}, color = [0 0 1];
   case {'white','w'}, color = [1 1 1];
   case {'black','k'}, color = [0 0 0];
   otherwise, color = [1 1 1];
  end % end switch statement
end
