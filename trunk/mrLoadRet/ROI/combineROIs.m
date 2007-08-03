function view = combineROIs(view,roi1,roi2,action)

% function view = combineROIs(view,roi1,roi2,action,[name])
%
% Logical combination (union, intersection, xor, set difference) of ROIs.
% Modifies roi1 by combining it with roi2
%
% roi1 and roi2 can be either ROI names, ROI numbers, or empty (for current ROI).
% Default: current ROI
%
% action must be empty or one of the following strings: 
%     'Intersection', 'Union', 'XOR', 'A not B'
% If empty, then 'Intersection' is used as the default.
%
% djh 8/2007

nROIs = viewGet(view,'numberofROIs');
if ieNotDefined('roi1')
  roi1 = viewGet(view,'currentROI');
end
if isstr(roi1)
  roi1 = viewGet(view,'roiNum',roi1);
end
if ~isnumeric(roi1) | (roi1 < 1) | (roi1 > nROIs)
  mrErrorDlg('Invalid ROI');
end
if ieNotDefined('roi2')
  roi2 = viewGet(view,'currentROI');
end
if isstr(roi2)
  roi2 = viewGet(view,'roiNum',roi2);
end
if ~isnumeric(roi2) | (roi2 < 1) | (roi2 > nROIs)
  mrErrorDlg('Invalid ROI');
end
if ieNotDefined('action')
  action = 'Intersection';
end

% Get coordinates
roiXform1 = viewGet(view,'roiXform',roi1);
roiVoxelSize1 = viewGet(view,'roiVoxelSize',roi1);
roiCoords1 = viewGet(view,'roiCoords',roi1);
roiXform2 = viewGet(view,'roiXform',roi2);
roiVoxelSize2 = viewGet(view,'roiVoxelSize',roi2);
roiCoords2 = viewGet(view,'roiCoords',roi2);

% Transform coords or roi2, using xformROIcoords to supersample the
% coordinates
if (size(roiCoords1,1)==3)
  roiCoords1(4,:) = 1;
end
if (size(roiCoords2,1)==3)
  roiCoords2(4,:) = 1;
end
roiCoords2 = round(xformROIcoords(roiCoords2,inv(roiXform1)*roiXform2,roiVoxelSize2,roiVoxelSize1));


% Transpose because matlab functions work on rows, not cols
coords1 = roiCoords1';
coords2 = roiCoords2';

% Combine
switch action
  case 'Intersection'
     newCoords = intersect(coords1,coords2,'rows');
  case 'Union'
     newCoords = union(coords1,coords2,'rows');
  case 'XOR'
     newCoords = setxor(coords,coords2,'rows');
  case 'A not B'
     newCoords = setdiff(coords1,coords2,'rows');
end

% Transpose back 
newCoords = newCoords';

% Select ROI and modify it
view = viewSet(view,'currentROI',roi1);
view = modifyROI(view,roiCoords1,roiXform1,roiVoxelSize1,0);
view = modifyROI(view,newCoords,roiXform1,roiVoxelSize1,1);
