function roiHistogram(view,overlayNum,scan,x,y,z,roi)
%
% roiHistogram(view,view,overlayNum,scan,x,y,z)
%
%        $Id$
% jb 10/12/2009
% plots the histogram of all non-NaN overlay values across voxels of all
% visible ROIs
% (non-NaN meaning between any boundaries that have been chosen from any
% overlay) 
%


% Error if no current ROI
if isempty(roi)
	disp('Voxel outside an ROI, using all visible ROIs');
   n_rois =  viewGet(view,'numberofROIs');
   if n_rois
      roi{1} = viewGet(view,'roi',1);
      if n_rois >1
         roi{1}.name = 'All ROIs';
         for i_roi = 2:n_rois
            temp_roi = viewGet(view,'roi',i_roi);
            roi{1}.coords = union(roi{1}.coords',temp_roi.coords','rows')';
            %here should check that ROIs are compatible
         end
      end
   else
      mrWarnDlg('(roiHistogram) No ROI currently loaded.');
   end
end

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{roinum}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

% get the analysis structure
analysis = viewGet(view,'analysis');
overlay_data = analysis.overlays(overlayNum).data{scan};

%First, get a mask of non-zero voxel representing the current overlay display
%This is taken from computeOverlay.m
numOverlays = viewGet(view,'numberofOverlays');
% Loop through overlays, filling in NaNs according to clip values.
for i_overlay = 1:numOverlays
    im = analysis.overlays(i_overlay).data{scan};
    clip = viewGet(view,'overlayClip',i_overlay);
    % Find pixels that are within clip
    if diff(clip) > 0
      pts = (im >= clip(1) & im <= clip(2));
    else
      pts = (im >= clip(1) | im <= clip(2));
    end
    % do not clip out for any points that are set to nan
    % this can happen if the current overlay does not
    % exist for this scan
    if i_overlay ~= overlayNum
      pts = pts | isnan(im);
    end
    % now make the mask
    if i_overlay ==1
       mask = pts;
    else
      mask = mask & pts;
    end
end
[mask_coordinates(:,1) mask_coordinates(:,2) mask_coordinates(:,3)] = ind2sub(size(mask),find(mask));

% select the window to plot into
fignum = selectGraphWin;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','roiHistogram');



for i_roi = 1:length(roi)
  subplot(1,length(roi),i_roi);
  coords = intersect(roi{i_roi}.scanCoords',mask_coordinates,'rows');
  overlay_values = overlay_data(sub2ind(size(overlay_data),coords(:,1), coords(:,2), coords(:,3)));
  
  hist(overlay_values,50);
  if strcmp(analysis.overlays(overlayNum).name,'ph')
     %rose(overlay_values,50);
     scale = axis;
     scale(1:2) = [0 2*pi];
     axis(scale);
     set(gca,'XTick',0:pi/2:2*pi)
     set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
  end
  xlabel(analysis.overlays(overlayNum).name);
  ylabel('Number of voxels');

  title([roi{i_roi}.name ' (' num2str(size(coords,1)) '/' num2str(size(roi{i_roi}.scanCoords,2)) ' voxels)']);
end





