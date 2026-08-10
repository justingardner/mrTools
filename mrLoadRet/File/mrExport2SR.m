function [overlayData,hdr] = mrExport2SR(viewNum, pathstr, baseNum)
% mrExport2SR.m
%
%      usage: [] = mrExprt2SR(viewNum, pathstr)
%         by: eli merriam
%       date: 03/20/07
%    purpose: exports a MLR overlay to a Nifti file in base space (compatible with SurfRelax for surfaces)
%             if baseNum is 0, the overlay is exported to scan space
%

% Get view
thisView = viewGet(viewNum,'view');

if ieNotDefined('baseNum')
  baseNum = viewGet(thisView,'currentBase');
end

% Get values from the GUI
scanNum = viewGet(thisView,'curscan');
overlayNum = viewGet(thisView,'currentOverlay');
overlayData = viewGet(thisView,'overlayData',scanNum,overlayNum);

if baseNum
  %transform values in base space
  [overlayData, new_base_voxel_size] = getBaseSpaceOverlay(thisView, overlayData, scanNum, baseNum);
  if isempty(overlayData)
    return
  end
  hdr = viewGet(thisView,'basehdr');
  if any(new_base_voxel_size ~= viewGet(thisView,'basevoxelsize',baseNum))
     hdr.pixdim = [0 new_base_voxel_size 0 0 0 0]';        % all pix dims must be specified here
     hdr.qform44 = diag([new_base_voxel_size 0]);
     hdr.sform44 = hdr.qform44;
  end
else
  hdr = viewGet(thisView,'niftihdr',scanNum);
  hdr.dim(5) = length(overlayNum);
end

hdr.datatype = 16; % make sure data are written as float32 (single)
hdr.scl_slope = 1;
hdr.scl_inter = 0;
if viewGet(thisView,'baseType',baseNum)==2 %for surfaces, leave as it was in the original mrExport2SR
  hdr.is_analyze = 1;
  hdr.endian = 'l';
end

% set the file extension
niftiFileExtension = mrGetPref('niftiFileExtension');
if isempty(niftiFileExtension)
  niftiFileExtension = '.img';
end

if ~ieNotDefined('pathstr') || nargout>0
  %write nifti file
  mlrImageWriteNifti(sprintf('%s%s',stripext(pathstr),niftiFileExtension),overlayData,hdr)
end
