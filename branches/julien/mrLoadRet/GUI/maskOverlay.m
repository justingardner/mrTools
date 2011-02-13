%
% maskOverlay.m
%
%   usage: [maskOverlay, overlayData] = maskOverlay(thisView, overlayList, scanList, overlayData)
%      by: jb, taken out of computeOverlay.m
%    date: 11/05/2010
% purpose: constructing a mask excluding clipped values by looping over all overlays of the current analysis
%           optionally outputs the non-masked overlays corresponding to overlayList
%           If a scan list is specified, overlays are retrieved from the view for these scans
%           Instead, a 4/5D array of data can be provided, the 4th and 5th dimensions being the overlays and scans respectively
%             in which case, scanList is ignored

function [maskOverlayData, overlayData] = maskOverlay(thisView,overlayList,scanList,overlayData)
  
if ieNotDefined('overlayData')
  if ieNotDefined('scanList')
    scanList = 1:viewGet(thisView,'nScans');
  end
  numOverlays = viewGet(thisView,'numberofOverlays');
  for iScan = 1:length(scanList)
    for iOverlay = 1:numOverlays
      overlayData{iOverlay,iScan} = viewGet(thisView,'overlayData',scanList(iScan),iOverlay);
    end
  end
else
  numOverlays = size(overlayData,4);
  overlayData= permute((num2cell(overlayData,[1 2 3])),[4 5 1 2 3]);
  scanList = 1;
end


% Loop through overlays, filling in NaNs according to clip values.
maskOverlayData = cell(1,length(scanList));
for iScan = 1:length(scanList)
  for iOverlay = 1:numOverlays
    thisOverlay = overlayData{iOverlay,iScan};
    clip = viewGet(thisView,'overlayClip',iOverlay);
    % Find pixels that are within clip
    if diff(clip) > 0
      pts = (thisOverlay >= clip(1) & thisOverlay <= clip(2));
    elseif diff(clip) < 0
      pts = (thisOverlay >= clip(1) | thisOverlay <= clip(2));
    else
      pts = false(size(thisOverlay));
    end
%     % do not clip out for any points that are set to nan    % JB: Can't figure out what this is for...
%     % this can happen if the current overlay does not       %
%     % exist for this scan                                   %
%     pts = pts | isnan(thisOverlay);                         %

    % now make the maskOverlayData
    if iOverlay ==1
      maskOverlayData{iScan} = pts;
    else
      maskOverlayData{iScan} = maskOverlayData{iScan}& pts;
    end
  end
end

if ieNotDefined('overlayList')
   overlayData = {};
else
   overlayData = overlayData(overlayList,:);
end

% if length(scanList)==1 && length(overlayList)==1
%    maskOverlayData = maskOverlayData{1};
%    overlayData = overlayData{1,1};
% end