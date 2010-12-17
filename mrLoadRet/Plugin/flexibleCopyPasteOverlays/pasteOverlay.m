function pasteOverlay(thisView, overlay)
% pasteOverlay.m
%
%        $Id$
%      usage: pasteOverlayBaseCoords(thisView, clipboard)
%         by: julien besle
%       date: 25/01/2010
%    purpose: pastes one or many MLR overlay(s) from the clipboard, and converts them to the current scan space

if ~isanalysis(viewGet(thisView,'analysis'))
    mrWarnDlg('(pasteOverlay) Overlays must be pasted into an analysis. Use Edit -> Analysis -> New Analysis.')
    return
end

nScans = viewGet(thisView,'nScans');
scanList = 1:nScans;
if isempty(scanList)
   mrWarnDlg('(pasteOverlay) Could not paste overlay because the group is empty (no scan)' );
   return;
end
while length(scanList)~= length(overlay(1).data)
   scanList=selectScans(thisView);
end
curGroupName = viewGet(thisView,'groupName');
curGroupNum =  viewGet(thisView,'groupNum',curGroupName);

for iOverlay = 1:length(overlay)
   [check thisOverlay] = isoverlay(overlay(iOverlay));
   if ~check
       mrErrorDlg('(pasteOverlay) Cannot paste. Clipboard does not contain a valid overlay. Use Edit -> Overlay -> Copy Overlay.')
   end
   fromGroupNum = viewGet(thisView,'groupNum',thisOverlay.groupName);
   scan2scan = viewGet(thisView,'scan2scan',1,curGroupNum,1,fromGroupNum);
   scandims = viewGet(thisView, 'scandims', 1);
   overlay(iOverlay).data = cell(1,nScans);
   if ~all(all((scan2scan - eye(4))<1e-6)) %check if we're in the scan space
      %transform values in current scan space
      [Ycoords,Xcoords,Zcoords] = meshgrid(1:scandims(2),1:scandims(1),1:scandims(3));
      for iScan = 1:length(scanList)
         fprintf(1,['Overlay ' num2str(iOverlay) ', Scan ' num2str(iScan) '\n']);
         overlay(iOverlay).data{scanList(iScan)} = getNewSpaceOverlay(thisOverlay.data{iScan}, scan2scan,Xcoords,Ycoords,Zcoords);
      end
   else
      overlay(iOverlay).data(scanList) = thisOverlay.data;
   end
   overlay(iOverlay).groupName =  curGroupName;
   overlay(iOverlay).name = [overlay(iOverlay).name ' (from ' thisOverlay.groupName ')'];
   %remove any specific reconcile function and params just in case, this overlay has to emancipate 
   overlay(iOverlay).reconcileFunction = 'defaultReconcileParams';
   overlay(iOverlay).params = [];

   thisView = viewSet(thisView,'newOverlay',overlay(iOverlay));
end


