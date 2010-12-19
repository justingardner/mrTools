% correctMonitorPosition.m
%
%      usage: monitorPositions = correctMonitorPosition(monitorPositions)
%         by: julien besle
%       date: 19/12/2010
%        $Id$
%    purpose: correct MonitorPosition property or root object
%             returned by get(0,'monitorPositions')
%             so that additional monitor positions are in the form [left bottom width height]
%               with the horizontal axis increasing from left to right
%               and the vertical axis increasing from bottom to top
%             consistent with the usual position vector form 
%            (as used by the 'position' property of figures)

function monitorPositions = correctMonitorPosition(monitorPositions)

if size(monitorPositions,1)>1
  if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64') 
    %windows machines return the position of secondary monitors in a 
    % [left top right bottom] form, where
    %       - the vertical origin is the top of the primary monitor
    %       - the vertical axis increases from top to bottom
    %       - the origin of the primary monitor (the top left) is (1,1) 
    %         (but we don't care about this because it's only one pixel)

  %the first position vector is the primary monitor. don't touch it
  % change the origin of the vertical axis, flip it and swap top and bottom
  monitorPositions(2:end,[4 2]) = (monitorPositions(2:end,[2 4]) - monitorPositions(1,4))*-1 +1;

  else
    keyboard
    %I don't know about other platforms... if they do what's described in the Matlab documentation
    %we shouldn't have to do anything
  end
  
  % compute width and heigth
  monitorPositions(2:end,3:4) = monitorPositions(2:end,3:4)+1 - monitorPositions(2:end,1:2);
end