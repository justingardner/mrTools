% copyOverlay.m
%
%      usage: overlay = copyOverlayBaseCoords(thisView)
%         by: julien besle, based on mrExprt2SR by eli merriam
%       date: 03/20/07
%    purpose: copies many MLR overlays to the clipboard
%        $Id$

function overlay = copyOverlay(thisView)

keepAsking = 1;
overlay = [];

while keepAsking
   scanList = 1:viewGet(thisView,'nScans');
   if length(scanList)>1
      scanList = selectInList(thisView,'scans');
   end
   if isempty(scanList)
      return;
   else
      keepAsking = 0;

      overlayList = 1:viewGet(thisView,'nOverlays');
      if length(overlayList)>1
         overlayList = selectInList(view,'overlays');
         if isempty(overlayList)
            if viewGet(thisView,'nScans')==1
               return
            else
               keepAsking=1;
            end
         end
      end

   end
end

overlay = [];
for iOverlay = 1:length(overlayList)
   overlay = copyFields(viewGet(thisView,'overlay',overlayList(iOverlay)),overlay,iOverlay);
   overlay(iOverlay).data = overlay(iOverlay).data(scanList);
end

return






