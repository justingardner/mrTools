function findMaxRoi(thisView,overlayNum,scanNum,x,y,z,roi)
%
% findMaxContiguousVoxels(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id$
% jb 04/12/2009
% finds the max and min value, as well as coordinates of all non-zero voxels in  ROI
% (non-zero meaning between any boundaries that have been chosen from any
% overlay)  
%
if isempty(roi)
   mrWarnDlg('(findMaxRoi) Please click in an ROI')
   return
end


baseNum = viewGet(thisView,'currentBase');

%First, get a logical mask of the current overlay display, as well as the overlay data
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);

overlayName = viewGet(thisView,'overlayName');

%mask the overlay
overlayData(~mask) = NaN;


base2scan = viewGet(thisView,'base2scan',scanNum,[],baseNum);
base2tal =  viewGet(thisView,'base2tal',baseNum);
if ~isempty(base2tal)
   scan2tal = base2tal/base2scan;
end
%scan2mag = viewGet(thisView,'scanXform',scanNum);
scanvoxelsize = viewGet(thisView,'scanvoxelsize');

disp(['Overlay: ' overlayName])

for iRoi = 1:length(roi)
   
  
%the following does not work, I suspect because of this shiftOriginXform thing
%    roi2mag = roi{iRoi}.xform;
%    roi2scan = scan2mag\roi2mag;
%do that instead:
   roi2scan = inv(viewGet(thisView,'scan2roi',viewGet(thisView,'roinum',roi{iRoi}.name)));
   scanRoiCoords = xformROIcoords(roi{iRoi}.coords,roi2scan,roi{iRoi}.voxelSize,scanvoxelsize);

   if ~isempty(scanRoiCoords)   
      scanRoiCoordsIndex = sub2ind(size(overlayData), scanRoiCoords(1,:)',  scanRoiCoords(2,:)',  scanRoiCoords(3,:)' );
      roiOverlayData = overlayData(scanRoiCoordsIndex);

      [max_value max_index] = max(roiOverlayData);
      [min_value min_index] = min(roiOverlayData);

      max_coordinates = scanRoiCoords(:,max_index);
      min_coordinates = scanRoiCoords(:,min_index);

      max_base_coordinates = (base2scan\max_coordinates);
      min_base_coordinates = (base2scan\min_coordinates);

      fprintf(1,['\tROI ' roi{iRoi}.name '(' num2str(numel(find(~isnan(roiOverlayData)))) '/' num2str(numel(roiOverlayData)) ' scan voxels)\n']);
      fprintf(1,['\t\tmax value :' num2str(max_value) '\n']);
      fprintf(1,['\t\tmax scan coordinates :' num2str(max_coordinates(1:3)') '\n']);
      fprintf(1,['\t\tmax base coordinates :' num2str(round(max_base_coordinates(1:3)')) '\n']);
      if ~isempty(base2tal)
         max_tal_coords = round(scan2tal*max_base_coordinates);
         fprintf(1,['\t\tmax Talairach coordinates: ' num2str(max_tal_coords(1:3)') '\n']);
      end

      fprintf(1,['\t\tmin value :' num2str(min_value) '\n']);
      fprintf(1,['\t\tmin scan coordinates :' num2str(min_coordinates(1:3)') '\n']);
      fprintf(1,['\t\tmin base coordinates :' num2str(round(min_base_coordinates(1:3)')) '\n']);
      if ~isempty(base2tal)
         min_tal_coords = round(scan2tal*min_base_coordinates);
         fprintf(1,['\t\tmin Talairach coordinates: ' num2str(min_tal_coords(1:3)') '\n']);
      end
      fprintf(1,'\n');

   else
      fprintf(1,[roi{iRoi}.name ' is empty\n\n']);
   end
end






