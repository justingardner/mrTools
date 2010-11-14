% applyWarpROI.m
%
%        $Id$ 
%      usage: applyWarpROI(thisView)
%         by: julien besle
%       date: 06/10/2010
%    purpose: applies non-linear FSL registration to chosen ROI(s)
%
function thisView = applyWarpROI(thisView)

useApplyWarp = 1; %The solution I chose is to create a volume of zeros and put ones where the ROI is, 
%then run this through applywarp using the spline coefficient  and nearest neighbour interpolation
% and finally fidn the coordinates of the transformed ROI (where the ones are)
% if useApplyWarp=0, the program asks for a warp field file (option -fout in fnirt) and apply the transformation
% but the sign of the fields in X and Y are not those expected and it does not work as well.

currentROIName = viewGet(thisView,'roiname');
needToRefresh=0;


if useApplyWarp
   %%%%%%%%%%%%%%%%%%%%%%%%%%%% using warping coefficient and applywarp
  keepAsking = 1;
  while keepAsking
    baseNum = [1 2];
    while length(baseNum)>1
      baseNum = selectBases(thisView,'Select the base space of the FNIRT warp coeffs');
      if length(baseNum)>1
        mrWarnDlg('Please select only one base');
      end
    end
    if isempty(baseNum)
      return;
    else
      while keepAsking
        roiList = selectROIs(thisView,'Select ROI(s) to warp');
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

  baseXform = viewGet(thisView,'basexform',baseNum);
  baseDims = viewGet(thisView,'basedims',baseNum); 
  baseHdr = viewGet(thisView,'basehdr',baseNum);

  fprintf('\n');
  for iRoi = 1:length(roiList)
    thisRoi = viewGet(thisView,'ROI',roiList(iRoi));
    %we require that ROI and scans be in the same space
    if any(any( (thisRoi.xform - baseXform)>1e-6))
      mrwarnDlg(['(applyWarpROI) Roi ' thisRoi.name ' is not in the requested base space, converting ...']);
      thisRoi.coords = getROICoordinates(thisView,roiList(iRoi),0,[],baseNum);
      thisRoi.xform = baseXform;
      thisRoi.sformCode = viewGet(thisView,'baseSformCode',baseNum);
      thisRoi.vol2mag = viewGet(thisView,'baseVol2mag',baseNum);
      thisRoi.vol2tal = viewGet(thisView,'baseVol2tal',baseNum);
      thisRoi.voxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
    end
    if isempty(thisRoi.coords)
      mrwarnDlg(['Roi ' thisRoi.name ' is empty, skipping ...']);
      fprintf('\n');
    else
      thisWarpedRoi = thisRoi;
      outsideVoxels = any(thisRoi.coords(1:3,:)< ones(3,size(thisRoi.coords,2)) |...
        thisRoi.coords(1:3,:)> repmat(baseDims',1,size(thisRoi.coords,2)));
      if any(outsideVoxels) 
        mrWarnDlg(['Some voxels in ' thisRoi.name ' are outside the base volume, removing...']);
        thisRoi.coords(:,outsideVoxels)=[];
      end
      warpIndices = sub2ind(baseDims(1:3),thisRoi.coords(1,:),thisRoi.coords(2,:),thisRoi.coords(3,:));
      roiVolume = zeros(baseDims(1:3));
      roiVolume(warpIndices) = iRoi;
      warpedRoiVolume = applyFSLWarp(roiVolume, [warpCoefPathname warpCoefFilename], 'tempInput.img', baseHdr, 'nearest');
      warpedCoordsIndices = find(warpedRoiVolume==iRoi);
      [x,y,z] = ind2sub(baseDims(1:3),warpedCoordsIndices);
      thisWarpedRoi.coords = [x y z ones(length(warpedCoordsIndices),1)]';
      thisView = viewSet(thisView,'newROI',thisWarpedRoi);
      disp(['warped ' thisRoi.name]);
      fprintf('\n');
      needToRefresh = 1;
    end
  end

else

   %%%%%%%%%%%%%%%%%%%%%%%% using warp field
   % I'm just keeping this code because I might want to come back to this solution, but stopped updating it
   keepAsking = 1;
   while keepAsking

      roiList = selectROIs(thisView,'Select ROI(s) to warp');
      if isempty(roiList)
         return;
      else
         keepAsking = 0;

         [warpCoefFilename warpCoefPathname] = uigetfile('*.img','FNIRT Warp field file');
         if isnumeric(warpCoefFilename)
            keepAsking=1;
         end
      end
   end

   [warpData, warpHdr] = cbiReadNifti([warpCoefPathname warpCoefFilename]);
   fieldSize = size(warpData);

   % transform the warpData in a 2D matrix where first dimension is a single index referring to the 3D volume
   % dx is stored in warpData(:,:,:,1); dy n (:,:,:,2); dz in (:,:,:,3)
   warpData = reshape(warpData, prod(fieldSize(1:3)),3);

   permuteMatrix = [-1 0 0;0 -1 0;0 0 1]; %for some reason, x an y axis are inverted in the warp fields...)
   for iRoi = 1:length(roiList)
      thisRoi = viewGet(thisView,'ROI',roiList(iRoi));
      if any(any( (thisRoi.xform - baseXform)>1e-6))
         mrwarnDlg(['Roi ' thisRoi.name ' is not in scan space, skipping ...']);
         fprintf('\n');
      elseif isempty(thisRoi.coords)
         mrwarnDlg(['Roi ' thisRoi.name ' is empty, skipping ...']);
         fprintf('\n');
      else
         thisConvertedRoi = thisRoi;
         %we assume that the ROI is in a space that's compatible with the reference image of the non-linear transformation
         warpIndices = sub2ind(fieldSize(1:3),thisRoi.coords(1,:),thisRoi.coords(2,:),thisRoi.coords(3,:));
         thisConvertedRoi.coords(1:3,:) = thisRoi.coords(1:3,:) + permuteMatrix* warpData(warpIndices,:)';  
         thisView = viewSet(thisView,'newROI',thisConvertedRoi);
         needToRefresh = 1;
      end
   end

end

if needToRefresh
 % select the same current roi
 thisView = viewSet(thisView,'curROI',viewGet(thisView,'roinum',currentROIName));
 refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
   
