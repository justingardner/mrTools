% setMLRViewAngle.m
%
%      usage: setMLRViewAngle(v)
%         by: justin gardner
%       date: 12/21/07
%    purpose: set the view angle for surfaces in MLR
%
function setMLRViewAngle(v)

% check arguments
if ~any(nargin == [1])
  help setMLRViewAngle
  return
end

% get the gui axis
fig = viewGet(v,'figNum');
gui = guidata(fig);

% get the tilt
baseTilt = viewGet(v,'baseTilt');

% careful not to make rotation and cam angle align
if baseTilt == 90,baseTilt = 89;end
if baseTilt == 270,baseTilt = 279;end

% rotate the view, using the *function* view
view(gui.axis,viewGet(v,'rotateSurface'),-baseTilt);

% change the camera position to avoid the volume
% flipping back and forth, another starnge matlab thing
if (baseTilt > 90) && baseTilt < 270
  camup(gui.axis,[0 0 -1]);
else
  camup(gui.axis,[0 0 1]);
end
