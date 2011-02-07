%
% maskOverlay.m
%
%   usage: [maskOverlay, overlayData] = maskOverlay(thisView, scanList, overlayList)
%      by: jb, taken out of computeOverlay.m
%    date: 11/05/2010
% purpose: constructing a mask excluding clipped/NaN values by looping over all overlays of the current analysis
%           optionally outputs the non-masked overlays corresponding to overlayList
%


function [maskOverlayData, overlayData] = maskOverlay(thisView,overlayList,scanList)


if ieNotDefined('overlayList')
   overlayList = [];
end
if ieNotDefined('scanList')
   scanList = viewGet(thisView,'nScans');
end
maskOverlayData = cell(length(scanList));
if ~isempty(overlayList)
   overlayData = cell(length(scanList),length(overlayList));
end

numOverlays = viewGet(thisView,'numberofOverlays');
% Loop through overlays, filling in NaNs according to clip values.
for iScan = 1:length(scanList)
   for iOverlay = 1:numOverlays
      thisOverlay = viewGet(thisView,'overlayData',scanList(iScan),iOverlay);
      clip = viewGet(thisView,'overlayClip',iOverlay);
      % Find pixels that are within clip
      if diff(clip) > 0
        pts = (thisOverlay >= clip(1) & thisOverlay <= clip(2));
      elseif diff(clip) < 0
        pts = (thisOverlay >= clip(1) | thisOverlay <= clip(2));
      else
        pts = false(size(thisOverlay));
      end
      [dump,index] = ismember(iOverlay,overlayList);
      if ~index
         pts = pts | isnan(thisOverlay);
      else
         % do not clip out for any points that are set to nan
         % this can happen if the current overlay does not
         % exist for this scan
         overlayData{iScan,index} = thisOverlay;
      end

      % now make the maskOverlayData
      if iOverlay ==1
       maskOverlayData{iScan} = pts;
      else
       maskOverlayData{iScan} = maskOverlayData{iScan}& pts;
      end
   end
end

if length(scanList)==1 && length(overlayList)==1
   maskOverlayData = maskOverlayData{1};
   overlayData = overlayData{1,1};
end