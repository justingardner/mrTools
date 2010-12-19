% copyOverlay.m
%
%      usage: overlay = copyOverlayBaseCoords(thisView)
%         by: julien besle, based on mrExprt2SR by eli merriam
%       date: 03/20/07
%    purpose: copies many MLR overlays to the clipboard
%        $Id: copyOverlay.m 1922 2010-12-14 08:39:43Z julien $

function overlay = copyOverlay(thisView)

keepAsking = 1;
overlay = [];

while keepAsking
   scanList = 1:viewGet(thisView,'nScans');
   if length(scanList)>1
      scanList = selectScans(thisView);
   end
   if isempty(scanList)
      return;
   else
      keepAsking = 0;

      overlayList = 1:viewGet(thisView,'nOverlays');
      if length(overlayList)>1
         overlayList = selectOverlays(thisView);
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






