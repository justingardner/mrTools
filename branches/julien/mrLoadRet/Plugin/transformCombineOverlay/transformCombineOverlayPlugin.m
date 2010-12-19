% transformCombineOverlayPlugin.m
%
%        $Id$ 
%      usage: transformCombineOverlayPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%
function retval = transformCombineOverlayPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help transformCombineOverlayPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(transformCombineOverlayPlugin) Need a valid view to install plugin'));
  else
    %first check if there is already an Overlay Menu
    hFig = viewGet(thisView,'fignum');
    hFigChildren = get(hFig,'Children');
    overlaysMenuInstalled = 0;
    for i = 1:length(hFigChildren)
      % if it is a menu item
      if isequal(get(hFigChildren(i),'Type'),'uimenu')&& strcmp(get(hFigChildren(i),'Label'),'Overlays')
        overlaysMenuInstalled=1;
        break;
      end
    end
    if ~overlaysMenuInstalled
      %if not, install it
      mlrAdjustGUI(thisView,'add','menu','Overlays','/Analysis');
    end
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Combine/Transform Overlays','/Overlays/','callback',@combineOverlaysCallback);
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds an item in Menu ''Overlays'' to combine/transform overlays using virtually any type of function';
 otherwise
   disp(sprintf('(transformCombineOverlayPlugin) Unknown command %s',action));
end

%------------------------- combineOverlaysCallback Function ------------------------------%
function combineOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
combineOverlays(thisView);
