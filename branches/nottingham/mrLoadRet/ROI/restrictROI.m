function view = restrictROI(view,ROInum,scan)

% function view = restrictROI(view,[ROInum],[scan])
%
% ROInum can be either a name, a number, or empty (for current ROI).
% Default: current ROI
%
% scan can be either a number or empty (for current scan)
% Default: current scan
%
% djh 9/2005

if ieNotDefined('ROInum')
  ROInum = viewGet(view,'currentroinum');
end
if ischar(ROInum)
  ROInum = viewGet(view,'roinum',ROInum);
end

if ieNotDefined('scan')
  scan = viewGet(view,'curscan');
end

ROIcoords = viewGet(view,'roiCoords',ROInum);
% Save prevCoords for undo
view = viewSet(view,'prevROIcoords',ROIcoords);

% Transform ROI roiScanCoords to overlay
roiScanCoords = round( viewGet(view,'scan2roi',ROInum,scan) \ ROIcoords); 

coordsInfo.base2overlay = eye(4);
coordsInfo.baseCoordsHomogeneous = roiScanCoords;
coordsInfo.baseDims = [size(ROIcoords,2) 1 1];

%find which voxels are not clipped in the current overlay(s) (in overlay space)
roiMask = maskOverlay(view,viewGet(view,'curOverlay'),scan,coordsInfo);
% keep the corresponding voxels in ROI space
ROIcoords = ROIcoords(:,any(roiMask{1},4));


view = viewSet(view,'roiCoords',ROIcoords,ROInum);

return

for overlayNum = 1:viewGet(view,'numberOfOverlays')
  overlayData = viewGet(view,'overlayData',scan,overlayNum);
  overlayClip = viewGet(view,'overlayClip',overlayNum);
  if ~isempty(overlayData) & ~isempty(ROIcoords)
    % Transform ROI roiScanCoords to overlay
    xform = inv(viewGet(view,'scan2roi',ROInum,scan));
    roiScanCoords = round(xform * ROIcoords); 
    overlaySize = size(overlayData);
    % Find roiScanCoords that are within the overlay size
    for v = 1:3
      indices = find((roiScanCoords(v,:) > 0) & (roiScanCoords(v,:) <= overlaySize(v)));
      roiScanCoords = roiScanCoords(:,indices);
      ROIcoords = ROIcoords(:,indices);
    end
    % Find overlay voxels that are within clip
    indices = sub2ind(overlaySize,roiScanCoords(1,:),roiScanCoords(2,:),roiScanCoords(3,:));
    data = overlayData(indices);
    if diff(overlayClip) > 0
      indices = find(data >= overlayClip(1) & data <= overlayClip(2));
    else
      indices = find(data >= overlayClip(1) | data <= overlayClip(2));
    end
    ROIcoords = ROIcoords(:,indices);
  else
    ROIcoords = [];
  end
  view = viewSet(view,'roiCoords',ROIcoords,ROInum);
end

return;
