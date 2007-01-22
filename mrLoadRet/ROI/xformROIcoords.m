function newcoords = xformROIcoords(coords,xform,inputVoxSize,outputVoxSize,sampRate)
%
% newcoords = xformROIcoords(coords,xform,inputVoxSize,outputVoxSize,[sampRate])
%
% Transforms ROI coords using Xform, supersampling in each dimension to
% accumulate partial volumes, then keeping only those voxels with partial
% volumes above thresh to maintain ROI volume.
% 
% coords: 4xN matrix of coordinates (x,y,z,1).
% xform: 4x4 homogeneous transform
% inputVoxSize: 3-vector, size of voxels (mm) in coords
% outputVoxSize: 3-vector, size of voxels (mm) in newCoords
% sampRate: 3-vector, supersampling rate for each dimension
%           default is odd number >= 4x ratio of inputVoxSize/outputVoxSize
%
% newcoords: 4xN matrix of (x,y,z,1) 
%
% djh, 8/98.  Modified from ROIcoords/transformROI.m in mrLoadRet-1
% 7/19/02 djh, Modified to maintain equal volumes
% 8/2005 djh, Updated to mrLoadRet-4.0

if ~exist('sampRate','var')
	sampRate = ceil(inputVoxSize ./ outputVoxSize) .* [4,4,4];
	sampRate = 2*floor(sampRate/2) + 1;
end

if isempty(coords)
	newcoords = [];
	return;
end

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin.
shiftXform = shiftOriginXform;
xform = inv(shiftXform) * xform * shiftXform;

% First, just transform the coords. This is insufficient because it will
% leave holes if the input voxels are larger than the output voxels.
%
tmpNewCoords = xform * coords;

% Find bounding (min and max) volume.
%
minPos = [min(tmpNewCoords(1,:));min(tmpNewCoords(2,:));min(tmpNewCoords(3,:));1];
maxPos = [max(tmpNewCoords(1,:));max(tmpNewCoords(2,:));max(tmpNewCoords(3,:));1];
minPos = floor(minPos)-[2,2,2,0]';
maxPos = ceil(maxPos)+[2,2,2,0]';
dims = (maxPos(1:3)-minPos(1:3)+ones(3,1))';

% Initialize accumulator for partial volume calculation, a vector
% of length appropriate to index the bounding volume.
%
accum = zeros(1,prod(dims));

% Calculate offsets that will be added within the loop to do the
% partial voluming.
%
xoffsets=[-.5+1/(2*sampRate(1)):1/sampRate(1):.5-1/(2*sampRate(1))];
yoffsets=[-.5+1/(2*sampRate(2)):1/sampRate(2):.5-1/(2*sampRate(2))];
zoffsets=[-.5+1/(2*sampRate(3)):1/sampRate(3):.5-1/(2*sampRate(3))];
% xoffsets=[0:1/sampRate(1):1-1/sampRate(1)];
% yoffsets=[0:1/sampRate(2):1-1/sampRate(2)];
% zoffsets=[0:1/sampRate(3):1-1/sampRate(3)];

% Divide alpha by prod(sampRate) to get partial volume for the
% supersampled voxels.
%
alpha = repmat(1/prod(sampRate),[1 size(coords,2)]);

% Loop through supersamples, transform them, and accumulate
% partial volume.
%
for ioff=1:length(xoffsets)
	xoff=xoffsets(ioff);
	for yoff=yoffsets
		for zoff=zoffsets
			% Add offset
			tmpNewCoords(1:3,:) = coords(1:3,:) + ...
				repmat([xoff;yoff;zoff],[1,size(coords,2)]);
			% Transform
			tmpNewCoords = xform * tmpNewCoords;
			% Round and subtract minPos
			tmpNewCoords(1:3,:) = round(tmpNewCoords(1:3,:)) - ...
				repmat(minPos(1:3),[1,size(tmpNewCoords,2)]);
			% Convert to indices
			indices = sub2ind(dims,tmpNewCoords(1,:),tmpNewCoords(2,:),tmpNewCoords(3,:));
			% Accumulate partial volume. Need to do it in a loop
			% instead of:
			%    accum(indices) = accum(indices) + alpha;
			% because an index can appear twice in indices and we want
			% to accumulate them both.
			for jj=1:length(indices)
				accum(indices(jj)) = accum(indices(jj)) + alpha(jj);
			end
		end
	end
end

% Build newROIcoords
%
[sortedAccum,indices] = sort(accum);
nonZeroSize = length(find(accum > 0));
newROIsize = round(prod(inputVoxSize)*size(coords,2) / prod(outputVoxSize));
newROIsize = min(nonZeroSize,newROIsize);
indices = indices(length(indices)-newROIsize+1:length(indices));
if ~isempty(indices)
  [newcoords1,newcoords2,newcoords3] = ind2sub(dims,indices);
  newcoords = [newcoords1; newcoords2; newcoords3; ones(size(newcoords1))];
  newcoords(1:3,:) = newcoords(1:3,:) + repmat(minPos(1:3),[1,length(indices)]);
else
  newcoords = [];
end

return;

%%%%%%%%%%%%%%
% Debug/test %
%%%%%%%%%%%%%%

coords1 = [0; 0; 0; 1];
coords2 = [1 2 3 4;
	       1 1 1 1;
	       1 1 1 1;
		   1 1 1 1];
xform = [1 0 0 0;
	     0 1 0 0;
	     0 0 1 0.5;
	     0 0 0 1];
inputVoxSize = [1,1,1];
outputVoxSize = [1,1,1/2];
xformROIcoords(coords1,xform,inputVoxSize,outputVoxSize)
xformROIcoords(coords2,xform,inputVoxSize,outputVoxSize)
