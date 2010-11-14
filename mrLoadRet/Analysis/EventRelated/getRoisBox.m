% getRoisBox - returns the coordinates of a cube in scan space exactly enclosing all the voxels in loaded ROIs
%
%        $Id$
%      usage: [d, roiVoxelIndices  ] = getRoisBox(thisView,scanNum,margin)
%         by: julien besle
%       date: 2010-04-26
%     inputs: 
%    outputs: scanBoxCoords: scan coordinates of the box
%                  whichRoi: 4D array in which each 3D array indicates the location of each roi in the box
%              marginVoxels: location of margin voxels in the box
%             roiScanCoords: each cell contains the coordinates of each roi in scan space

function [scanBoxCoords, whichRoi, marginVoxels, roiScanCoords] = getRoisBox(thisView,scanNum,margin,roiNums)

if ieNotDefined('margin')
   margin = 0;
end
if ieNotDefined('roiNums')
   roiNums=1:length(thisView.ROIs);
end
nRois = length(roiNums);

scanDims = viewGet(thisView,'scandims',scanNum);

%find the scan coordinates of all rois
minX = inf;
maxX = 0;
minY = inf;
maxY = 0;
minZ = inf;
maxZ = 0;
for iRoi = 1:nRois
   roiScanCoords{iRoi} = getROICoordinates(thisView,roiNums(iRoi),scanNum);
   %remove voxels that are outside the scan
   if ~isempty(roiScanCoords{iRoi})
     roiScanCoords{iRoi} = roiScanCoords{iRoi}(:,all(roiScanCoords{iRoi}(1:3,:)>0));
     %make sure we don't have the same voxel twice for this roi
     roiScanCoords{iRoi} = unique(roiScanCoords{iRoi}','rows')';
     %find the coordinates of the subset box including all the voxels from all the ROIs
     minX = max(min(minX,min(roiScanCoords{iRoi}(1,:))-margin),1);
     maxX = min(max(maxX,max(roiScanCoords{iRoi}(1,:))+margin),scanDims(1));
     minY = max(min(minY,min(roiScanCoords{iRoi}(2,:))-margin),1);
     maxY = min(max(maxY,max(roiScanCoords{iRoi}(2,:))+margin),scanDims(2));
     minZ = max(min(minZ,min(roiScanCoords{iRoi}(3,:))-margin),1);
     maxZ = min(max(maxZ,max(roiScanCoords{iRoi}(3,:))+margin),scanDims(3));
   end
end
scanBoxCoords = [minX maxX;minY maxY;minZ maxZ];
scanBoxDims = diff(scanBoxCoords')+1;

%for each roi, find the indices of the voxels in the box
whichRoi=zeros([scanBoxDims nRois]);
for iRoi = 1:nRois
   whichRoiTemp = zeros(scanBoxDims);
   if ~isempty(roiScanCoords{iRoi})
     XcoordsInBox = roiScanCoords{iRoi}(1,:) - minX + 1;
     YcoordsInBox = roiScanCoords{iRoi}(2,:) - minY + 1;
     ZcoordsInBox = roiScanCoords{iRoi}(3,:) - minZ + 1;
     whichRoiTemp(sub2ind(scanBoxDims,XcoordsInBox',YcoordsInBox',ZcoordsInBox'))=1;
   end
   whichRoi(:,:,:,iRoi) = whichRoiTemp;
end

%if margin, find voxels around Rois
marginVoxels = any(whichRoi,4); %where in the box there is any data
marginVoxels = convn(marginVoxels,ones(2*margin+1),'same')>0;
marginVoxels = logical(marginVoxels-any(whichRoi,4));


