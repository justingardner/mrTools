% applyFslFnirtWarpPlugin.m
%
%        $Id: applyFslFnirtWarpPlugin.m 1925 2010-12-14 23:16:51Z julien $ 
%      usage: DefaultPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%
function retval = applyFslFnirtWarpPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help applyFslFnirtWarpPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(applyFslFnirtWarpPlugin) Need a valid view to install plugin'));
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
   retval = 'Adds items in menus ''Overlays'' and ''ROI'' to apply FSL FNIRT warp spline coefficients to overlays and/or ROIs';
 otherwise
   disp(sprintf('(applyFslFnirtWarpPlugin) Unknown command %s',action));
end

%------------------------- applyWarpOverlaysCallback Function ------------------------------%
function applyWarpOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
applyFSLWarpOverlays(thisView);

%------------------------- applyWarpROICallback Function ------------------------------%
function applyWarpROICallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
applyFSLWarpROI(thisView);

