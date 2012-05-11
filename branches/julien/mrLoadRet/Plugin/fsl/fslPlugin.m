% fslPlugin.m
%
%        $Id$ 
%      usage: DefaultPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%
function retval = fslPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help fslPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(fslPlugin) Need a valid view to install plugin'));
  else
% % %     %first check if there is already an Overlay Menu
% % %     hFig = viewGet(thisView,'fignum');
% % %     hFigChildren = get(hFig,'Children');
% % %     overlaysMenuInstalled = 0;
% % %     for i = 1:length(hFigChildren)
% % %       % if it is a menu item
% % %       if isequal(get(hFigChildren(i),'Type'),'uimenu')&& strcmp(get(hFigChildren(i),'Label'),'Overlays')
% % %         overlaysMenuInstalled=1;
% % %         break;
% % %       end
% % %     end
% % %     if ~overlaysMenuInstalled
% % %       %if not, install it
      mlrAdjustGUI(thisView,'add','menu','overlaysMenu','/Analysis','label','Overlays','tag','overlaysMenu');
% %     end
    %install overlays menu Item
    mlrAdjustGUI(thisView,'add','menu','Apply FNIRT non-linear warps to overlays','/Overlays/','callback',@applyWarpOverlaysCallback,'tag','applyFnirtOverlayMenuItem');
    %install overlays menu Item
    mlrAdjustGUI(thisView,'add','menu','Apply FNIRT non-linear warps to ROIs','/ROI/Combine','callback',@applyWarpROICallback,'tag','applyFnirtRoiMenuItem');
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds FSL functionnalities: (1) add items in menus ''Overlays'' and ''ROI'' to apply FSL FNIRT warp spline coefficients to overlays and/or ROIs';
 otherwise
   disp(sprintf('(fslPlugin) Unknown command %s',action));
end

%------------------------- applyWarpOverlaysCallback Function ------------------------------%
function applyWarpOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
fslApplyWarpOverlays(thisView);

%------------------------- applyWarpROICallback Function ------------------------------%
function applyWarpROICallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
fslApplyWarpROI(thisView);

