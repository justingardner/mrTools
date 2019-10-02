% makeROIsExactlyContiguous.m
%
%      usage: transformedRois = makeROIsExactlyContiguous(rois,margin,<kernelType>)
%         by: julien besle
%       date: 11/01/2011
%
%    purpose: attribute voxels shared by two ROIs to the closest of these two ROIs, with 
%             the two resulting two ROIs then becoming exactly contiguous.
%             This function will run with more than two ROIs and will apply to all pairs of ROIs
%             This assumes that no voxel is shared by more than two ROIs. If this is the case,  
%             then results will depend on the order in which the ROIs are passed

function rois = makeROIsExactlyContiguous(rois)

if ~ismember(nargin,[1])
  help expandROI;
  return
end

if numel(rois) < 2
  mrWarnDlg('(makeROIsExactlyContiguous) You need to provide at least 2 ROIs ');
  return
end

% Check that transformation matrices are identical for all ROIs
for iRoi = 2:length(rois)
  if any(any( (rois(1).xform - rois(iRoi).xform) > 10e-6))
    mrWarnDlg('(makeROIsExactlyContiguous) All ROIs must be converted to the same space (set the roiSpace option to something other than ''Native'').');
    return
  end
end


for iRoi = 1:length(rois)
  for jRoi = iRoi+1:length(rois)
    coords1 = rois(iRoi).coords';
    coords2 = rois(jRoi).coords';
    [commonCoordinates, indexROI1, indexROI2] = intersect(coords1,coords2,'rows');
    if ~isempty(commonCoordinates)
      %remove common coordinates from ROIs 1 and 2
      coords1 = setdiff(coords1,commonCoordinates,'rows');
      coords2 = setdiff(coords2,commonCoordinates,'rows');
      %attribute common coordinates to one or the other ROI depending on distance
      belongsToROI1 = false(size(commonCoordinates,1),1);
      for iCoords = 1:size(commonCoordinates,1)
        %compute distance between these coordinates and all coordinates unique to either both ROI
        distanceCoords1 = sqrt(sum((repmat(commonCoordinates(iCoords,1:3),size(coords1,1),1) - coords1(:,1:3)).^2,2));
        distanceCoords2 = sqrt(sum((repmat(commonCoordinates(iCoords,1:3),size(coords2,1),1) - coords2(:,1:3)).^2,2));
        %identify closest ROI
        if min(distanceCoords1) < min(distanceCoords2)
          belongsToROI1(iCoords) = true;
        end
      end
      % delete coords that belong to the other ROI
      rois(iRoi).coords(:,indexROI1(~belongsToROI1'))=[];
      rois(jRoi).coords(:,indexROI2(belongsToROI1'))=[];
    end
  end
end


