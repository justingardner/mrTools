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
baseVoxelSize = viewGet(thisView,'baseVoxelSize');

%First, get a logical mask of the current overlay display, as well as the overlay data
nOverlays = length(overlayNum);
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);
mask = mask{1};
overlayData = overlayData{1};

%mask the overlay
overlayData(~mask) = NaN;

%get ROIs
roiList = viewGet(thisView,'visibleRois');
cRoi = 0;
for iRoi = roiList
  cRoi = cRoi+1;
  rois{cRoi} = viewGet(thisView,'roi',iRoi);
  base2roi(:,:,cRoi) = viewGet(thisView, 'base2roi',iRoi);
end
nRois = length(rois);

figure('name',['Base anatomy: ' baseName]);
nTissueVoxels = zeros(nRois,nOverlays);
nVeinsVoxels = zeros(nRois,nOverlays);
nAllVoxels = zeros(nRois,nOverlays);
meanTissue = zeros(nRois,nOverlays);
meanVeins = zeros(nRois,nOverlays);

cOverlay=0;
for iOverlay = overlayNum
  overlayMax = 0;
  overlayMin = 0;
  cOverlay=cOverlay+1;
  %transform values in base space
  thisOverlayData = getBaseSpaceOverlay(thisView, overlayData(:,:,:,cOverlay), scanNum, baseNum);
  dataSize = size(thisOverlayData);
  overlayName = viewGet(thisView,'overlayName',iOverlay);

  for iRoi = 1:nRois
     hSubplot(cOverlay,iRoi) = subplot(nOverlays,nRois,(cOverlay-1)*nRois+iRoi);
     hold on

     scancoords = xformROIcoords(rois{iRoi}.coords,inv(base2roi(:,:,iRoi)),rois{iRoi}.voxelSize,baseVoxelSize);
     % we need to remove any coordinate that might fall outside the base anatomy
     if ~isempty(scancoords)
        outside_voxels = find(scancoords(1,:)<1 | scancoords(1,:)>dataSize(1) |...
                          scancoords(2,:)<1 | scancoords(2,:)>dataSize(2) |...
                          scancoords(3,:)<1 | scancoords(3,:)>dataSize(3) );
        scancoords(:,outside_voxels) = [];
     end

     if ~isempty(scancoords)
        %RoiData
        scanRoiCoordsIndex = sub2ind(dataSize, scancoords(1,:)',  scancoords(2,:)',  scancoords(3,:)' );
        roiOverlayData = thisOverlayData(scanRoiCoordsIndex);
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
           nTissueVoxels(iRoi,cOverlay) = length(tissue_index);
           nVeinsVoxels(iRoi,cOverlay) = length(vein_index);
           nAllVoxels(iRoi,cOverlay) = size(roiOverlayData,1);

           title(sprintf('%s (%d voxels, %.2f %% veins)', rois{iRoi}.name, nTissueVoxels(iRoi,cOverlay), nVeinsVoxels(iRoi,cOverlay)./nAllVoxels(iRoi,cOverlay)*100),'interpreter','none');

           %calculate mean and range of overlay for the veins and for the tissue
           %voxels
           fprintf(1,[rois{iRoi}.name '\n']);
           fprintf(1,['Veins: ' num2str(length(vein_index)) '/' num2str(nAllVoxels(iRoi,cOverlay)) ' voxels\t']);
           [string,meanVeins(iRoi,cOverlay)] = meanAndRange(roiOverlayData(vein_index));
            strings{1} = ['Veins: ' string];
           fprintf(1,['Tissue: ' num2str(length(tissue_index)) '/' num2str(nAllVoxels(iRoi,cOverlay)) ' voxels\t']);
           [string,meanTissue(iRoi,cOverlay)] = meanAndRange(roiOverlayData(tissue_index));
            strings{2} = ['Tissue: ' string];
            xlabel(strings);

        else
           scatter_plot = 1;
           scatter(baseData,roiOverlayData,2);
           title(sprintf('%s (%d voxels)', rois{iRoi}.name, size(roiOverlayData,1)),'interpreter','none');

           xlabel(baseName,'interpreter','none');
           ylabel(overlayName,'interpreter','none');

           base_max = max(overlayMax,max(baseData));
           base_min = min(overlayMin,min(baseData));
        end

        overlayMax = max(overlayMax,max(roiOverlayData));
        overlayMin = min(overlayMin,min(roiOverlayData));
     end
  end

  for iRoi = 1:nRois
     if ~isempty(nTissueVoxels(iRoi,cOverlay)) && ~isempty(nVeinsVoxels(iRoi,cOverlay))
        subplot(hSubplot(cOverlay,iRoi))
        scale = axis;
        if scatter_plot
           scale = [base_min base_max overlayMin overlayMax];
        else
           scale(3:4) = [overlayMin overlayMax];
        end
        axis(scale);
     end
  end
end

           nTissueVoxels
           nVeinsVoxels
           nAllVoxels
           meanTissue
           meanVeins
