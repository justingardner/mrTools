% makeROIsExactlyContiguous.m
%
%      usage: transformedRois = makeROIsExactlyContiguous(rois)
%         by: julien besle
%       date: 11/01/2011
%
%    purpose: make two or more ROIs mutually exclusive
%
function rois = makeROIsExactlyContiguous(rois)

if ~ismember(nargin,[1])
  help makeROIsExactlyContiguous;
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
  for jRoi = 1:length(rois)
    if iRoi ~= jRoi
      % first ensure that all voxel coordinates are rounded and unique
      coords1 = unique(round(rois(iRoi).coords'),'rows');
      coords2 = unique(round(rois(jRoi).coords'),'rows');
      [commonCoordinates, indexROI1, indexROI2] = intersect(coords1,coords2,'rows'); % indexROI1, indexROI2 used to be used below
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
  %       % delete coords that belong to the other ROI
  %       rois(iRoi).coords(:,indexROI1(~belongsToROI1'))=[];
  %       rois(jRoi).coords(:,indexROI2(belongsToROI1'))=[];
          % instead of deleting common voxels, replace all voxels in each ROI by its unique voxels
          % and the common voxels that have been attributed to it
          % (replacing is necessary because coordinates might have been rounded and duplicates removed)
          rois(iRoi).coords = [coords1; commonCoordinates(belongsToROI1,:)]';
          rois(jRoi).coords = [coords2; commonCoordinates(~belongsToROI1,:)]';
      end
    end
  end
end


