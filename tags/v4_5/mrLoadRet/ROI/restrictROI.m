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
if isstr(ROInum)
	ROInum = viewGet(view,'roinum',ROInum);
end

if ieNotDefined('scan')
	scan = viewGet(view,'curscan');
end

ROI = viewGet(view,'roi',ROInum);
ROIcoords = viewGet(view,'roiCoords',ROInum);
ROIxform = viewGet(view,'roiXform',ROInum);
ROIvoxelSize = viewGet(view,'roiVoxelSize',ROInum);

for overlayNum = 1:viewGet(view,'numberOfOverlays')
	overlayData = viewGet(view,'overlayData',scan,overlayNum);
	overlayClip = viewGet(view,'overlayClip',overlayNum);
	overlayXform = viewGet(view,'overlayXform',scan);
    overlayVoxelSize = viewGet(view,'scanVoxelSize',scan);
	if ~isempty(overlayData)
        % Transform ROI coords to overlay
        coords = xformROIcoords(ROIcoords,inv(overlayXform)*ROIxform,overlayVoxelSize,ROIvoxelSize);
		overlaySize = size(overlayData);
		% Find coords that are within the overlay size
		for v = 1:3
			indices = find((coords(v,:) > 0) & (coords(v,:) <= overlaySize(v)));
			coords = coords(:,indices);
		end
		% Find overlay voxels that are within clip
		indices = sub2ind(overlaySize,coords(1,:),coords(2,:),coords(3,:));
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