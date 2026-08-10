% fillHolesInROI.m
%
%      usage: transformedRoi = fillHolesInROI(roi,<connectivity>)
%         by: julien besle
%       date: 06/04/2021
%
%    purpose: fills holes in ROI using a simple algorithm 
%      input:   - connectivity: number of neighboring voxels each voxel can have in 2D or 3D (default = 6)
%                               6: contiguous voxel faces in 3D
%                               18: contiguous voxel faces and edges in 3D
%                               26: contiguous voxel faces, edges and corners in 3D
%                               4: contiguous faces in X-Y plane
%                               8: contiguous faces and edges in X-Y plane

function roi = fillHolesInROI(roi,connectivity)

if ~ismember(nargin,[1 2])
  help fillHolesInROI;
  return
end

if ieNotDefined('connectivity')
  connectivity = 6;
end
if ~ismember(connectivity,[4 6 8 18 26])
  mrWarnDlg(['(expandROI) unknown connectivity value ' connectivity]);
  roi=[];
  return
end

boxCoords = [min(roi.coords(1:3,:),[],2)-[1 1 1]'  max(roi.coords(1:3,:),[],2)+[1 1 1]'];

%shift coordinates so that the boxes starts at 1 on all dimensions
voxelShift = -boxCoords(:,1)+1;

boxCoords = boxCoords+repmat(voxelShift,1,2);
roiCoords = roi.coords(1:3,:)+repmat(voxelShift,1,size(roi.coords,2));

volume = zeros(boxCoords(:,2)');
volume(sub2ind(boxCoords(:,2)',roiCoords(1,:),roiCoords(2,:),roiCoords(3,:)))=1;

% if trim
%   volume = 1-volume;
% end

switch(connectivity)
  case 4
    kernel = zeros(3,3,3);
    kernel(:,:,2) = [0 1 0;1 1 1;0 1 0];
  case 8
    kernel = zeros(3,3,3);
    kernel(:,:,2) = ones(3,3,1);
  case 6
    kernel(:,:,1) = [0 0 0;0 1 0;0 0 0];
    kernel(:,:,2) = [0 1 0;1 1 1;0 1 0];
    kernel(:,:,3) = [0 0 0;0 1 0;0 0 0];
  case 18
    kernel(:,:,1) = [0 1 0;1 1 1;0 1 0];
    kernel(:,:,2) = [1 1 1;1 1 1;1 1 1];
    kernel(:,:,3) = [0 1 0;1 1 1;0 1 0];
  case 26
    kernel = ones(3,3,3);
   
end

volume = logical(volume);
holesRemain = true;
while holesRemain
  volume2 = logical(convn(volume,kernel,'same')); % fills voxels neighboring any filled voxel
  volume2 = ~volume2;
  volume2 = logical(convn(volume2,kernel,'same')); % empty voxels neighboring any empty voxel
  volume2 = ~volume2;
  
  if isequal(volume,volume2) %until all holes have been filled
    holesRemain = false;
  else
    volume = volume2;
  end
end

[newCoordsX,newCoordsY,newCoordsZ] = ind2sub(boxCoords(:,2)',find(volume));
roi.coords = [newCoordsX-voxelShift(1) newCoordsY-voxelShift(2) newCoordsZ-voxelShift(3)]';


