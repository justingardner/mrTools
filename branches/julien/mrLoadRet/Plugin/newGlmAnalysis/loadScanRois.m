% loadScanRois - load timeseries corresponding to the voxels of one or several loaded Rois
%
%        $Id$
%      usage: [d, roiVoxelIndices  ] = loadScanRois(thisView,scanNum,roiList)
%         by: julien besle
%       date: 2010-04-26
%     inputs: 
%    outputs: [
%
%    purpose: load timeseries corresponding to the voxels of one or several Rois (wrapper around loadScan.m)
%

function [d, roiVoxelIndices  ] = loadScanRois(thisView,scanNum,roiList)

%get the smallest box in the scan including all the rois in the list
[subsetBox, whichRoi]  = getRoisBox(thisView,scanNum,0,roiList);  
isInARoi = find(any(whichRoi,4));

d = loadScan(thisView,scanNum,[],subsetBox(3,:),mrGetPref('defaultPrecision'),subsetBox(1,:),subsetBox(2,:));
%only keep the data that's in the ROIs
d.data = reshape(d.data,[prod(d.dim(1:3)),1,1,d.dim(4)]);
d.data = d.data(isInARoi,:,:,:);
d.dim = size(d.data);

%for each roi, find the indices of the voxels in the data
for iRoi = 1:length(roiList)
   [dump, roiVoxelIndices{iRoi}] = intersect(isInARoi,find(whichRoi(:,:,:,iRoi)));
end

%checking for NaNs in the tseries
isnanVoxels = any(isnan(d.data),4);
if nnz(isnanVoxels)
   mrWarnDlg('(loadScanRois) removing voxels containing NaNs in timeseries data')
   subsetDims = diff(subsetBox,1,2)+1;
   newIndices = zeros(size(d.data,1),1);
   newIndices(~isnanVoxels) = 1:sum(~isnanVoxels);
   for iRoi = 1:length(roiList)
      newRoiVoxelIndices = newIndices(roiVoxelIndices{iRoi});
      if ~all(newRoiVoxelIndices)
         fprintf(1,['ROI ' num2str(iRoi) ' (' viewGet(thisView,'roiName',roiList(iRoi)) '):\n']);
         [isnanXCoords,isnanYCoords,isnanZCoords] = ind2sub(subsetDims',roiVoxelIndices{iRoi}(~newRoiVoxelIndices));
         disp(subsetBox(:,1)'-1+[isnanXCoords,isnanYCoords,isnanZCoords]);
      end
      roiVoxelIndices{iRoi} = nonzeros(newRoiVoxelIndices);
   end
   d.data = d.data(~isnanVoxels,:,:,:);
end
      

d.dim = size(d.data);


