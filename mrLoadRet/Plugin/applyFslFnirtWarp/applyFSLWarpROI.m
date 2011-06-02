% applyWarpROI.m
%
%        $Id: applyWarpROI.m 2119 2011-05-09 22:46:58Z julien $ 
%      usage: applyWarpROI(thisView)
%         by: julien besle
%       date: 06/10/2010
%    purpose: applies non-linear FSL registration to chosen ROI(s)
%
function thisView = applyFSLWarpROI(thisView)

if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrWarnDlg('(applyFSLWarpROI) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
  return
end

useApplyWarp = 1; %The solution I chose is to create a volume of zeros and put ones where the ROI is, 
%then run this through applywarp using the spline coefficient  and nearest neighbour interpolation
% and finally find the coordinates of the transformed ROI (where the ones are)
% if useApplyWarp=0, coordinates are converted directly using warp fields created using fnirtfileutils
% but this does not preserve the volume of the ROI and takes more time.

currentROIName = viewGet(thisView,'roiname');
needToRefresh=0;

keepAsking = 1;
while keepAsking
  baseNum = [1 2];
  while length(baseNum)>1
    baseNum = selectInList(thisView,'bases','Select the input volume of the FNIRT warp coeffs (or equivalent)');
    if length(baseNum)>1
      mrWarnDlg('(applyFSLWarpROI) Please select only one base');
    end
  end
  if isempty(baseNum)
    return;
  else
    while keepAsking
      roiList = selectInList(thisView,'rois','Select ROI(s) to warp');
      if isempty(roiList)
         break;
      else
         keepAsking = 0;
         [warpCoefFilename warpCoefPathname] = uigetfile('*.img','FNIRT Warp spline coefficient file');
         if isnumeric(warpCoefFilename)
            keepAsking=1;
         end
      end
    end
  end
end
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%% using warping coefficient and applywarp

baseXform = viewGet(thisView,'basexform',baseNum);
baseDims = viewGet(thisView,'basedims',baseNum); 
baseHdr = viewGet(thisView,'basehdr',baseNum);
baseVoxelSize = viewGet(thisView,'baseVoxelSize',baseNum);

fprintf('\n');
for iRoi = 1:length(roiList)
  thisRoi = viewGet(thisView,'ROI',roiList(iRoi));
  %we require that ROI and scans be in the same space, but we keep a resolution of the original space
  %this is done by removing the voxel size from the transformation matrix 
  if any(any( (thisRoi.xform - baseXform)>1e-6))
    
    mrWarnDlg(['(applyFSLWarpROI) ROI ' thisRoi.name ' is not in the requested base space, converting ...']);
    
    scalingFactor = ceil(baseVoxelSize./thisRoi.voxelSize);
    %transformation matrix between roi space and fnirt input volume space (base)
    roi2base = diag([scalingFactor 1])*(baseXform\thisRoi.xform);
    
    thisRoi.coords = xformROIcoords(thisRoi.coords,roi2base,thisRoi.voxelSize,baseVoxelSize./scalingFactor);
%     thisRoi.xform = baseXform;
%     thisRoi.sformCode = viewGet(thisView,'baseSformCode',baseNum);
%     thisRoi.vol2mag = viewGet(thisView,'baseVol2mag',baseNum);
%     thisRoi.vol2tal = viewGet(thisView,'baseVol2tal',baseNum);
%     thisRoi.voxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
    baseDims = scalingFactor.*baseDims;
    baseHdr.pixdim(2:4) = baseHdr.pixdim(2:4)./scalingFactor';
    baseHdr.sform44 = baseHdr.sform44*diag([1./scalingFactor 1]);
    baseHdr.qform44 = baseHdr.qform44*diag([1./scalingFactor 1]);

  else
    scalingFactor = [1 1 1];
  end
  %for some reason, there might be NaNs in the roi coords, remove those
  thisRoi.coords = thisRoi.coords(:,~any(isnan(thisRoi.coords(1:3,:)),1));
  if isempty(thisRoi.coords)
    mrwarnDlg(['Roi ' thisRoi.name ' is empty, skipping ...']);
    fprintf('\n');
  else
    thisWarpedRoi = thisRoi;
    outsideVoxels = any(thisRoi.coords(1:3,:)< ones(3,size(thisRoi.coords,2)) |...
    thisRoi.coords(1:3,:)> repmat(baseDims',1,size(thisRoi.coords,2)));
    if any(outsideVoxels) 
      mrWarnDlg(['(applyFSLWarpROI) Some voxels in ' thisRoi.name ' are outside the input FNIRT volume, removing...']);
      thisRoi.coords(:,outsideVoxels)=[];
    end
    if useApplyWarp
      warpIndices = sub2ind(baseDims(1:3),thisRoi.coords(1,:),thisRoi.coords(2,:),thisRoi.coords(3,:));
      roiVolume = zeros(baseDims(1:3));
      roiVolume(warpIndices) = iRoi;
      warpedRoiVolume = applyFSLWarp(roiVolume, [warpCoefPathname warpCoefFilename], 'tempInput.img', baseHdr, 'nearest');
      warpedCoordsIndices = find(warpedRoiVolume==iRoi);
      [x,y,z] = ind2sub(baseDims(1:3),warpedCoordsIndices);
      thisWarpedRoi.coords = [x y z ones(length(warpedCoordsIndices),1)]';
      

    else
      %%%%%%%%%%%%%%%%%%%%%%%% using warp field, but that doesn't take the ROI volume into account
      thisWarpedRoi.coords = applyFSLWarpCoords([thisRoi.coords;ones(1,size(thisRoi.coords,2))], baseVoxelSize./scalingFactor, [.5 .5 .5], [warpCoefPathname warpCoefFilename], 'tempInput.img', baseHdr);
      %no debugged since modified function to stay in roi space...
    end
    if any(any( (thisRoi.xform - baseXform)>1e-6))
      %convert back to roi space
      thisWarpedRoi.coords = xformROIcoords(thisWarpedRoi.coords,inv(roi2base),baseVoxelSize./scalingFactor,thisRoi.voxelSize);
    end
    thisView = viewSet(thisView,'newROI',thisWarpedRoi);
    disp(['(applyFSLWarpROI) Warped ROI ' thisRoi.name]);
    fprintf('\n');
    needToRefresh = 1;
  end
end



if needToRefresh
 % select the same current roi
 thisView = viewSet(thisView,'curROI',viewGet(thisView,'roinum',currentROIName));
 refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
   
