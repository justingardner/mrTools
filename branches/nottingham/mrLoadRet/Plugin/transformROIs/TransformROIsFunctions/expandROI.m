% expandROI.m
%
%        $Id: expandROI.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: transformedRoi = expandROI(roi,margin,<kernelType>)
%         by: julien besle
%       date: 11/01/2011
%
%    purpose: expands roi coordinates in 3D by margin voxels
%      input:   - margin: number of voxels that will be added (removed if negative) around the ROI.
% kernelType:   - (optional) type of kernel to convolve the ROI with: 
%                     sphere (default), cube, disc or square.
%                     the margin is the radius of the sphere/disc or the half-side of the square/cube 
%                     discs and squares are in the X-Y plane  

function roi = expandROI(roi,margin,kernelType)

if ~ismember(nargin,[2 3])
  help expandROI;
  return
end

if ieNotDefined('kernelType')
  kernelType = 'sphere';
end
if ~ismember(kernelType,{'sphere','cube','disc','square'})
  mrWarnDlg(['(expandROI) unknown kernel type ' kernelType]);
  roi=[];
end

trim = margin<0;
margin = abs(margin);

boxCoords = [min(roi.coords,[],2)-margin  max(roi.coords,[],2)+margin];

%shift coordinates so that the boxes starts at 1 on all dimensions
voxelShift = -boxCoords(:,1)+1;

boxCoords = boxCoords+repmat(voxelShift,1,2);
roiCoords = roi.coords+repmat(voxelShift,1,size(roi.coords,2));

volume = zeros(boxCoords(:,2)');
volume(sub2ind(boxCoords(:,2)',roiCoords(1,:),roiCoords(2,:),roiCoords(3,:)))=1;

if trim
  volume = 1-volume;
end

switch(kernelType)
  case 'sphere'
    [sphereX,sphereY,sphereZ] = ndgrid(-margin:margin,-margin:margin,-margin:margin);
    kernel = zeros(2*margin+1,2*margin+1,2*margin+1);
    kernel((sphereX.^2+sphereY.^2+sphereZ.^2)<margin^2)=1;
  case 'cube'
    kernel = ones(2*margin,2*margin,2*margin);
  case 'square'
    kernel = ones(2*margin,2*margin);
  case 'disc'
    [discX,discY] = ndgrid(-margin:margin,-margin:margin);
    kernel = zeros(2*margin+1,2*margin+1);
    kernel((discX.^2+discY.^2)<margin^2)=1;
end

volume = logical(convn(volume,kernel,'same'));
if trim
  volume = ~volume;
end
[newCoordsX,newCoordsY,newCoordsZ] = ind2sub(boxCoords(:,2)',find(volume));
roi.coords = [newCoordsX-voxelShift(1) newCoordsY-voxelShift(2) newCoordsZ-voxelShift(3)]';


