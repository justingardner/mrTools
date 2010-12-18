% flexibleCopyPasteOverlaysPlugin.m
%
%        $Id$ 
%      usage: flexibleCopyPasteOverlaysPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%    purpose: replaces behaviour of copy/paste overlay in Edit menu
%
function retval = flexibleCopyPasteOverlaysPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help flexibleCopyPasteOverlaysPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(flexibleCopyPasteOverlaysPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    mlrAdjustGUI(thisView,'set','/Edit/Overlay/Copy Overlay','Callback',@copyOverlayCallback);
    mlrAdjustGUI(thisView,'set','/Edit/Overlay/Paste Overlay','Callback',@pasteOverlayCallback);

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Allows copy/paste of several overlays and automatic interpolation to the destination space';
 otherwise
   disp(sprintf('(flexibleCopyPasteOverlaysPlugin) Unknown command %s',action));
end

% --------------------------------------------------------------------
function copyOverlayCallback(hObject, ~)
mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
clipboard = copyOverlay(thisView); %calls copyOverlay which asks which overlay and for which scans to copy an returns all the copied overlays
if ~isempty(clipboard)
  MLR.clipboard = clipboard;
end

% --------------------------------------------------------------------
function pasteOverlayCallback(hObject, ~)
mrGlobals;
viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');
pasteOverlay(thisView, MLR.clipboard);
refreshMLRDisplay(viewNum);
