function tseriesROI = tseriesROI(view, groupNum, roiList, scanList, varargin)
%
% tseries = tseriesROI(view, groupNum, roiList, scanList, [param1], [value1], [param2], [value2])
% 
% Extracts tseries for each voxel in an ROI and for each scan, excluding junkframes.
%
% groupNum: default current group
% roiList: vector of ROI numbers (default: [1:nROIs])
% scanList: vector of scan numbers (default: [1:nscans])
% other arguments as in percentTSeries  
%
% tseries: cell array (nROIs x nScans), each element of which is nFrames x nVoxels matrix
%
% djh 9/2005

if ieNotDefined('groupNum')
	groupNum = viewGet(view,'currentGroup');
end
view = viewSet(view,'currentGroup',groupNum);
if ieNotDefined('roiList')
	roiList = [1:viewGet(view,'numberofROIs')];
end
if ieNotDefined('scanList')
	scanList = [1:viewGet(view,'nscans')];
end

% Initialize cell array
tseriesROI = cell(length(roiList),length(scanList));

% Loop through ROIs, scans and slices, extracting tseries and putting them
% into a cell array.

for iROI = 1:length(roiList)
	roi = roiList(iROI);
	roiCoords = viewGet(view,'roiCoords',roi);
	roiXform = viewGet(view,'roiXform',roi);
	roiVoxelSize = viewGet(view,'roiVoxelSize',roi);

	for iscan = 1:length(scanList)
		scan = scanList(iscan);
		scanXform = viewGet(view,'scanXform',scan);
		scanVoxelSize = viewGet(view,'scanVoxelSize',scan);
		scanDims = viewGet(view,'scanDims',scan);
		sliceDims = scanDims([1,2]);
        junkframes = viewGet(view,'junkFrames',scan);
        nframes = viewGet(view,'nFrames',scan);

		% Transform ROI it to the coordinates of the scan
		coords = round(xformROIcoords(roiCoords,inv(scanXform)*roiXform,roiVoxelSize,scanVoxelSize));
		if isempty(coords)
			slices = [];
		else
			slices = coords(3,:);
			slices = unique(slices);
		end
		
		disp(['Loading tseries from scan ',num2str(scan)]);
		for islice = 1:length(slices)
			slice = slices(islice);
			roiIndices = find(coords(3,:) == slice);
			sliceCoords = coords([1,2],roiIndices);
			% Find coords that are within the FOV
			for v = 1:2
				indices = find((sliceCoords(v,:) > 0) & (sliceCoords(v,:) <= sliceDims(v)));
				sliceCoords = sliceCoords(:,indices);
			end
			sliceIndices = unique(sub2ind(sliceDims,sliceCoords(1,:),sliceCoords(2,:)));
			tseries = loadTSeries(view,scan,slice);
			subtSeries = tseries(junkframes+1:junkframes+nframes,sliceIndices);
			ptSeries = percentTSeries(subtSeries,varargin{:});
			tseriesROI{iROI,iscan} = [tseriesROI{iROI,iscan}, ptSeries];
		end
	end
end

return

% Test
tseries = tseriesROI(MLR.views{1},[],[1:2],[1:3]);
