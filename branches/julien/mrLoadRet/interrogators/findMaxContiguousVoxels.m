function findMaxContiguousVoxels(thisView,overlayNum,scanNum,x,y,z,roi)
%
% findMaxContiguousVoxels(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id$
% jb 04/12/2009
% finds the max and min value, as well as coordinates of all non-zero voxels contiguous to the selected non-zero voxel
% (non-zero meaning between any boundaries that have been chosen from any
% overlay)  
%


baseNum = viewGet(thisView,'currentBase');


% get the analysis structure
analysis = viewGet(thisView,'analysis');

%First, get a mask of non-zero voxel representing the current overlay display
mask = maskOverlay(thisView,overlayNum,scanNum);

base2scan = viewGet(thisView,'base2scan',scanNum,[],baseNum);

base2tal =  viewGet(thisView,'base2tal',baseNum);

%Now find contiguous voxels, starting from the selected non-zeros voxel
selected_voxel_index = sub2ind(size(mask),x,y,z);
if mask(x,y,z)
%    while ~isempty(contiguous_next_slice)
%       
%       %the idea here was to go up the slices from the current one and find
%       %the voxels connected to the selected contiguous blob 
%       %until there is none or the last slice is reached
%       %then go down the slices repeating the same procedure
%       %it would go up and down until no more connected voxels are found
%       % but I'll use bwconncomp.m from the image processing toolbox instead
%       
%    end
   cc = bwconncomp(mask,6);
   for i_object = 1:cc.NumObjects
      if ismember(selected_voxel_index,cc.PixelIdxList{i_object})
         contiguous_voxels = cc.PixelIdxList{i_object};
         break;
      end
   end
   
   [max_value max_index] = max(analysis.overlays(overlayNum).data{scanNum}(contiguous_voxels));
   [min_value min_index] = min(analysis.overlays(overlayNum).data{scanNum}(contiguous_voxels));
   
   [max_coordinates(1) max_coordinates(2) max_coordinates(3)] = ind2sub(size(mask),contiguous_voxels(max_index));
   [min_coordinates(1) min_coordinates(2) min_coordinates(3)] = ind2sub(size(mask),contiguous_voxels(min_index));
   
   max_base_coordinates = (inv(base2scan)*[max_coordinates 1]')';
   min_base_coordinates = (inv(base2scan)*[min_coordinates 1]')';

   disp([num2str(length(contiguous_voxels)) ' contiguous voxels']);
   disp(['max value :' num2str(max_value)]);
   disp(['max scan coordinates :' num2str(max_coordinates)]);
   disp(['max base coordinates :' num2str(round(max_base_coordinates(1:3)))]);
   if ~isempty(base2tal)
      max_tal_coords = round(base2tal*max_base_coordinates')';
      disp(['max Talairach coordinates: ' num2str(max_tal_coords(1:3))]);
   end
   
   disp(['min value :' num2str(min_value)]);
   disp(['min scan coordinates :' num2str(min_coordinates)]);
   disp(['min base coordinates :' num2str(round(min_base_coordinates(1:3)))]);
   if ~isempty(base2tal)
      min_tal_coords = round(base2tal*min_base_coordinates')';
      disp(['min Talairach coordinates: ' num2str(min_tal_coords(1:3))]);
   end

else
   mrWarnDlg('(findMaxContiguousVoxel) Please select a non-zero voxel');
end






