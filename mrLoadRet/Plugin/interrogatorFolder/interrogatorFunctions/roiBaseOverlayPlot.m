function roiBaseOverlayPlot(thisView,overlayNum,scanNum,x,y,z,roi)
%
%      
%        $Id$
%    authors: rs, jb 26/01/2010
%       goal: plot the distribution of current overlay values against current base anatomy values
%              uses box plot is the base anatomy only has two values (0 and 1)
%              use as an interrogate function in mrLoadRet
%

if isempty(thisView.ROIs)
   mrWarnDlg('(roiBaseOverlayPlot) Please load at least one ROI')
   return
end

% get underlying image intensity
baseNum = viewGet(thisView,'currentBase');
base = viewGet(thisView,'baseanatomy');
baseName = viewGet(thisView,'basename');
basevoxelsize = viewGet(thisView,'basevoxelsize');

%First, get a logical mask of the current overlay display, as well as the overlay data
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);
overlay_name = viewGet(thisView,'overlayName');

%mask the overlay
overlayData(~mask) = NaN;

%transform values in base space
overlayData = getBaseSpaceOverlay(thisView, overlayData, scanNum, baseNum);
dataSize = size(overlayData);


figure('name',['Overlay: ' overlay_name ' - Base anatomy: ' baseName]);
overlay_max = 0;
overlay_min = 0;
nTissueVoxels = zeros(length(thisView.ROIs),1);
nVeinsVoxels = zeros(length(thisView.ROIs),1);
nAllVoxels = zeros(length(thisView.ROIs),1);
meanTissue = zeros(length(thisView.ROIs),1);
meanVeins = zeros(length(thisView.ROIs),1);

for i_roi = 1:length(thisView.ROIs)
   h_subplot(i_roi) = subplot(1,length(thisView.ROIs),i_roi);
   hold on
   r = thisView.ROIs(i_roi);

   xform = viewGet(thisView, 'base2roi',i_roi);
   roivoxelsize = viewGet(thisView,'roivoxelsize',i_roi);
   scancoords = xformROIcoords(r.coords,inv(xform),roivoxelsize,basevoxelsize);
   % we need to remove any coordinate that might fall outside the base anatomy
   if ~isempty(scancoords)
      outside_voxels = find(scancoords(1,:)<1 | scancoords(1,:)>dataSize(1) |...
                        scancoords(2,:)<1 | scancoords(2,:)>dataSize(2) |...
                        scancoords(3,:)<1 | scancoords(3,:)>dataSize(3) );
      scancoords(:,outside_voxels) = [];

      %RoiData
      scanRoiCoordsIndex = sub2ind(dataSize, scancoords(1,:)',  scancoords(2,:)',  scancoords(3,:)' );
      roiOverlayData = overlayData(scanRoiCoordsIndex);
      baseData = base.data(scanRoiCoordsIndex);

      nonNanValues = find(~isnan(roiOverlayData));
      roiOverlayData = roiOverlayData(nonNanValues);
      baseData = baseData(nonNanValues);

      %see if base anatomy data are logical or continuous
      if length(unique(baseData))<3
         if length(unique(baseData))==1
            vein_value = 1;
            tissue_value = 0;
         else
            vein_value = max(baseData);
            tissue_value = min(baseData);
         end
         scatter_plot = 0;

         %calculate number of voxel which are veins 
         vein_index = find(baseData==vein_value);
         tissue_index = find(baseData==tissue_value);
         labels = cell(size(roiOverlayData));
         labels(tissue_index) = {'Tissue'};
         labels(vein_index) = {'Veins'};
         
         if ~isempty(vein_index) && ~isempty(tissue_index)
            boxplot(roiOverlayData,labels,'grouporder',{'Veins','Tissue'},'notch','on')
         end
         nTissueVoxels(i_roi) = length(tissue_index);
         nVeinsVoxels(i_roi) = length(vein_index);
         nAllVoxels(i_roi) = size(roiOverlayData,1);

         title(sprintf('%s (%d voxels, %.2f %% veins)', r.name, nTissueVoxels(i_roi), nVeinsVoxels(i_roi)./nAllVoxels(i_roi)*100));

         %calculate mean and range of overlay for the veins and for the tissue
         %voxels
         fprintf(1,[r.name '\n']);
         fprintf(1,['Veins: ' num2str(length(vein_index)) '/' num2str(nAllVoxels(i_roi)) ' voxels\t']);
         [string,meanVeins(i_roi)] = meanAndRange(roiOverlayData(vein_index));
          strings{1} = ['Veins: ' string];
         fprintf(1,['Tissue: ' num2str(length(tissue_index)) '/' num2str(nAllVoxels(i_roi)) ' voxels\t']);
         [string,meanTissue(i_roi)] = meanAndRange(roiOverlayData(tissue_index));
          strings{2} = ['Tissue: ' string];
          xlabel(strings);

      else
         scatter_plot = 1;
         scatter(baseData,roiOverlayData,2);
         title(sprintf('%s (%d voxels)', r.name, size(roiOverlayData,1)));

         Xlabel(baseName);
         Ylabel(overlay_name);

         base_max = max(overlay_max,max(baseData));
         base_min = min(overlay_min,min(baseData));
      end

      overlay_max = max(overlay_max,max(roiOverlayData));
      overlay_min = min(overlay_min,min(roiOverlayData));
   end
end

         nTissueVoxels
         nVeinsVoxels
         nAllVoxels
         meanTissue
         meanVeins

for i_roi = 1:length(thisView.ROIs)
   if ~isempty(nTissueVoxels(i_roi)) && ~isempty(nVeinsVoxels(i_roi))
      subplot(h_subplot(i_roi))
      scale = axis;
      if scatter_plot
         scale = [base_min base_max overlay_min overlay_max];
      else
         scale(3:4) = [overlay_min overlay_max];
      end
      axis(scale);
   end
end
