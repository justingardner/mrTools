% combineTransformOverlaysPlugin.m
%
%        $Id: combineTransformOverlaysPlugin.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: combineTransformOverlaysPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%
function retval = combineTransformOverlaysPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help combineTransformOverlaysPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(combineTransformOverlaysPlugin) Need a valid view to install plugin'));
  else
    mlrAdjustGUI(thisView,'add','menu','overlaysMenu','/Analysis','label','Overlays','tag','overlaysMenu');
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Combine/Transform Overlays','/Overlays/','callback',@combineOverlaysCallback,'tag','transformOverlaysMenuItem');
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds an item in Menu ''Overlays'' to combine/transform overlays using virtually any type of function';
 otherwise
   disp(sprintf('(combineTransformOverlaysPlugin) Unknown command %s',action));
end

%------------------------- combineOverlaysCallback Function ------------------------------%
function combineOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
combineOverlays(thisView);
