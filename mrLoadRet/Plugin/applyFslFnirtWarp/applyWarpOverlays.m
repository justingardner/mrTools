function applyWarpOverlays(thisView)
% applyWarpOverlays(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%   applies non-linear FSL registration to overlay in current view and current analysis
%
% jb 23/07/2010
%
% $Id$ 

if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrWarnDlg('(applyWarpROI) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
  return
end

keepAsking = 1;
while keepAsking
   scanList = 1:viewGet(thisView,'nScans');
   if length(scanList)>1
      scanList = selectScans(thisView);
   end
   if isempty(scanList)
      return;
   else
      keepAsking = 0;

      overlayList = 1:viewGet(thisView,'nOverlays');
      if length(overlayList)>1
         overlayList = selectOverlays(thisView);
         if isempty(overlayList)
            if viewGet(thisView,'nScans')==1
               return
            else
               keepAsking=1;
            end
         end
      end

   end
end

[warpCoefFilename warpCoefPathname] = uigetfile('*.img','FNIRT spline coefficient file');
if isnumeric(warpCoefFilename)
   return
end

warning('off','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
interpMethod = mrGetPref('interpMethod');
scanDims = viewGet(thisView,'scandims'); %we assume scans are all in the same space
scanHdr = viewGet(thisView,'niftiHdr');

overlaysNumber = length(overlayList);
for iOverlay = 1:overlaysNumber
   overlays{iOverlay} = viewGet(thisView,'overlay',overlayList(iOverlay));
   if isfield(overlays{iOverlay},'alphaOverlay') && ~isempty(overlays{iOverlay}.alphaOverlay)
      overlaysNumber = overlaysNumber+1;
      alphaOverlayName = overlays{iOverlay}.alphaOverlay;
      overlays{overlaysNumber} = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum',alphaOverlayName));
      overlays{iOverlay}.alphaOverlay = ['warp(' alphaOverlayName ') - ' warpCoefFilename];
      disp(['(applyWarpOverlays) Alpha overlay ''' alphaOverlayName ''' will also be warped']);
   end
end

for iOverlay = 1:overlaysNumber
   newOverlay = overlays{iOverlay};
   overlayIndices = false(1,length(scanList));
   for iScan = 1:length(scanList)
      overlayIndices(iScan) = ~isempty(overlays{iOverlay}.data(scanList));
   end
   if any(overlayIndices)
      newOverlay.data = cell(1,length(overlays{iOverlay}.data));
      %transform the non empty overlay in a 4D matrix
      data = permute(reshape(cell2mat(overlays{iOverlay}.data(scanList(overlayIndices))),[scanDims([1 2]) sum(overlayIndices) scanDims(3)]), [1 2 4 3]);
      warpedData = applyFSLWarp(data, [warpCoefPathname warpCoefFilename], 'tempInput.img', scanHdr, interpMethod);
      newOverlay.data(scanList(overlayIndices)) = mat2cell(warpedData,scanDims(1),scanDims(2),scanDims(3),ones(1,sum(overlayIndices)));
      newOverlay.name = ['warp(' overlays{iOverlay}.name ') - ' warpCoefFilename];
      %reset the clip values. This is due to the fact that warped data can go beyond the min and max of the original data (even when using nearest-neighboor, probably due to a bug of FSL
      newOverlay.clip = [min(warpedData(:)) max(warpedData(:))];
      thisView = viewSet(thisView,'newoverlay',newOverlay);

   end
end
warning('on','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');

refreshMLRDisplay(thisView.viewNum);



