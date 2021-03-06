% setMLRViewAngle.m
%
%      usage: setMLRViewAngle(v,<axisHandle>)
%         by: justin gardner
%       date: 12/21/07
%    purpose: set the view angle for surfaces in MLR. Defaults to working on MLR figure, but you can pass
%             it an axis handle and it will work on that instead.
%
function setMLRViewAngle(v,axisHandle)

% check arguments
if ~any(nargin == [1 2])
  help setMLRViewAngle
  return
end

% get the gui axis
if nargin < 2
  fig = viewGet(v,'figNum');
  gui = guidata(fig);
  axisHandle = gui.axis;
else
  if ~ishandle(axisHandle)
    disp(sprintf('(setMLRViewAngle) Based in axisHandle is not valid'));
    return
  end
end


% get the tilt
baseTilt = viewGet(v,'baseTilt');

% careful not to make rotation and cam angle align
if ismember(rem(baseTilt,180),[-90 90])
  baseTilt = baseTilt-0.1;
end


% rotate the view, using the *function* view
view(axisHandle,360-viewGet(v,'rotateSurface'),-baseTilt);

% change the camera position to avoid the volume
% flipping back and forth, another strange matlab thing
if cos(baseTilt/180*pi)<0
  camup(axisHandle,[0 0 -1]);
else
  camup(axisHandle,[0 0 1]);
end
