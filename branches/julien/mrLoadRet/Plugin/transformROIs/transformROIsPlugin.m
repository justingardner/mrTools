% transformROIsPlugin.m
%
%        $Id: transformROIsPlugin.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: transformROIsPlugin(action,<thisView>)
%         by: julien besle
%       date: 11/01/2011
%
function retval = transformROIsPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help transformROIsPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(transformROIsPlugin) Need a valid view to install plugin'));
  else
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Transform','/ROIs/Combine','callback',@transformROIsCallback,'label','Transform','tag','transformRoiMenuItem');
    mlrAdjustGUI(thisView,'add','menu','Transform','/ROI/Combine','callback',@transformROIsCallback,'label','Transform','tag','transformRoiMenuItem');
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds an item in Menu ''ROI'' to transform ROIs with pre-defined or custom transformation functions';
 otherwise
   disp(sprintf('(transformROIsPlugin) Unknown command %s',action));
end

%------------------------- transformROIsCallback Function ------------------------------%
function transformROIsCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
transformROIs(thisView);
