function view = drawROI(view,descriptor,sgn)

% view = drawROI(view,[descriptor],[sgn])
%
% Adds/deletes to/from the current ROI based on user interaction.
%
% descriptor: option for how the new coordinates are to be specified.
% Current options are:
%    'rectangle'[default]
%
% sgn: If sgn~=0 [default, adds user-specified coordinates to selected ROI
% in current slice. If sgn==0, removes those coordinates from the ROI.
%
% djh, 8/2005 (modified from mrLoadRet-3.1)

% Error if no current ROI
curROInum = viewGet(view,'currentROI');
if isempty(curROInum)
	mrErrorDlg('No ROI currently selected.');
end

if ieNotDefined('descriptor')
	descriptor = 'rectangle'
end
if ieNotDefined('sgn')
	sgn = 1;
end

% Select main axes of view figure for user input
fig = viewGet(view,'figNum');
gui = guidata(fig);
set(fig,'CurrentAxes',gui.axis);

% baseCoords contains the mapping from pixels in the displayed slice to
% voxels in the current base volume.
baseCoords = viewGet(view,'cursliceBaseCoords');
baseSliceDims = [size(baseCoords,1),size(baseCoords,2)];
if isempty(baseCoords)
	mrErrorDlg('Load base anatomy before drawing an ROI');
end

switch descriptor
	
	case 'rectangle'
		% Get region from user.
		region = round(ginput(2));
		
		% Note: ginput hands them back in x, y order (1st col is x and 2nd col is
		% y). But we use them in the opposite order (y then x), so flip 'em.
		region = fliplr(region);
		% Check if outside image
		if (min(region(:,1))< 1 | max(region(:,1))>baseSliceDims(1) | ...
				min(region(:,2))< 1 | max(region(:,2))>baseSliceDims(2))
			mrWarnDlg('Must choose rect entirely within image boundaries');
			return;
		end
		% Make sure 2nd value is larger than the 1st.
		for i=1:2
			if region(2,i) < region(1,i)
				region(:,i)=flipud(region(:,i));
			end
		end
		
		% Extract coordinates in base reference frame
		baseX = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],1);
		baseY = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],2);
		baseZ = baseCoords([region(1,1):region(2,1)],[region(1,2):region(2,2)],3);
		coords = [baseX(:)'; baseY(:)'; baseZ(:)'; ones(1,prod(size(baseX)))];
        
    case 'polygon'
        % Get polygon region using matlab's roipoly function
        polyIm = roipoly;
        
        % Compute image coordinates
        polyImIndices = find(polyIm);
        [x,y] = ind2sub(size(polyIm),polyImIndices);
        
        % Extract coordinates in base reference frame
		baseX = baseCoords(x,y,1);
		baseY = baseCoords(x,y,2);
		baseZ = baseCoords(x,y,3);
		coords = [baseX(:)'; baseY(:)'; baseZ(:)'; ones(1,prod(size(baseX)))];
		
	otherwise
		mrErrorDlg(['Invalid descriptor: ',descriptor]);
end

% Modify ROI
baseNum = viewGet(view,'currentBase');
xform = viewGet(view,'basexform',baseNum);
voxelSize = viewGet(view,'baseVoxelSize',baseNum);
view = modifyROI(view,coords,xform,voxelSize,sgn);


