% mlrExportROI.m
%
%        $Id$ 
%      usage: mlrExportROI(v,saveFilename,<hdr>)
%         by: justin gardner
%       date: 07/14/09
%    purpose: Export an ROI to a nifti image. Uses current roi
%             and current base in view to export. Pass in a nifti
%             header as hdr argument if you want to use a different header
%
function mlrExportROI(v,saveFilename,varargin)

% check arguments
if nargin < 2
  help mlrExportROI
  return
end

% optional arguments
getArgs(varargin,{'hdr=[]'});

% get the roi we are being asked to export
roiNum = viewGet(v,'currentroi');
if isempty(roiNum)
  mrWarnDlg('(mlrExportROI) No current ROI to export');
  return
end
  
if ischar(saveFilename)
  saveFileName = {saveFilename};
end
if ~isequal(length(roiNum),length(saveFilename))
  mrWarnDlg('(mlrExportROI) number of file names must be identical to number of ROIs');
  return
end

% get the base nifti header
passedInHeader = false;
if ~isempty(hdr)
  passedInHeader = true;
else
  hdr = viewGet(v,'basehdr');
  if isempty(hdr)
    mrWarnDlg('(mlrExportROI) Could not get base anatomy header');
    return
  end
end

baseCoordMap = viewGet(v,'basecoordmap');
baseType = viewGet(v,'basetype');
if ~isempty(baseCoordMap) && baseType==1  %for flats, use basecoordmap 
  [~,baseCoords,baseCoordsHomogeneous] = getBaseSlice(v,1,3,viewGet(v,'baseRotate'),viewGet(v,'curBase'),baseType);
  % make sure that baseCoords are rounded (they may not be
  % if we are working with a baseCoordMap's flat map
  baseDims = size(baseCoords);
  baseDims = baseDims ([1 2 4]);
  
  baseCoordsHomogeneous = reshape(baseCoordsHomogeneous,4,prod(baseDims));
  baseCoordsHomogeneous = round(baseCoordsHomogeneous);
  baseCoordsLinear = mrSub2ind(baseCoordMap.dims,baseCoordsHomogeneous(1,:),baseCoordsHomogeneous(2,:),baseCoordsHomogeneous(3,:));

  % estimate voxel size (taken from getBaseOverlay, assuming mask is invarioant to rotation, which it should be since it is a flat map)
  oldBaseVoxelSize=viewGet(v,'basevoxelsize',viewGet(v,'curBase'));
  Xcoords0Mask = permute(baseCoords(:,:,1,:)==0,[1 2 4 3]);
  Xcoords0Mask = convn(Xcoords0Mask,ones(5,5,5),'same'); %expand the mask a bit to make sure we don't include any edge voxels
  XcoordsNaN = permute(baseCoords(:,:,1,:),[1 2 4 3]); 
  XcoordsNaN(Xcoords0Mask>0)=NaN;
  YcoordsNaN = permute(baseCoords(:,:,2,:),[1 2 4 3]);
  YcoordsNaN(Xcoords0Mask>0)=NaN;
  ZcoordsNaN = permute(baseCoords(:,:,3,:),[1 2 4 3]);
  ZcoordsNaN(Xcoords0Mask>0)=NaN;
  newBaseVoxelSize(1) = oldBaseVoxelSize(1)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,1).^2 + diff(YcoordsNaN,1,1).^2 + diff(ZcoordsNaN,1,1).^2))));
  newBaseVoxelSize(2) = oldBaseVoxelSize(2)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,2).^2 + diff(YcoordsNaN,1,2).^2 + diff(ZcoordsNaN,1,2).^2))));
  newBaseVoxelSize(3) = oldBaseVoxelSize(3)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,3).^2 + diff(YcoordsNaN,1,3).^2 + diff(ZcoordsNaN,1,3).^2))));
  if any(newBaseVoxelSize ~= oldBaseVoxelSize)
   hdr.pixdim = [0 newBaseVoxelSize 0 0 0 0]';        % all pix dims must be specified here
   hdr.qform44 = diag([newBaseVoxelSize 0]);
   hdr.sform44 = hdr.qform44;
  end
  
else
  baseDims = hdr.dim(2:4)';
end

if ~passedInHeader
  b = viewGet(v,'base');
  % if the orientation has been changed in loadAnat, undo that here.
  if ~isempty(b.originalOrient)
    % create a data structure that has all 0's
    d = zeros(baseDims);
    % convert into mlrImage
    [d h] = mlrImageLoad(d,hdr);
    % convert the orientation back to original
    [d h] = mlrImageOrient(b.originalOrient,d,h);
    % covert back to nifti
    hdr = mlrImageGetNiftiHeader(h);
  end
end

for iRoi = 1:length(roiNum)
  % tell the user what is going on
  disp(sprintf('(mlrExportROI) Exporting ROI to %s with dimensions set to match base %s: [%i %i %i]',saveFilename{iRoi},viewGet(v,'baseName'),baseDims(1),baseDims(2),baseDims(3)));

  % create a data structure that has all 0's
  d = zeros(baseDims);

  % get  roi coordinates in base coordinates
  roiBaseCoords = getROICoordinates(v,roiNum(iRoi),0);

  % check roiBaseCoords
  if isempty(roiBaseCoords)
    mrWarnDlg('(mlrExportROI) This ROI does not have any coordinates in the base');
    return
  end

  % make sure we are inside the base dimensions
  xCheck = (roiBaseCoords(1,:) >= 1) & (roiBaseCoords(1,:) <= hdr.dim(2));
  yCheck = (roiBaseCoords(2,:) >= 1) & (roiBaseCoords(2,:) <= hdr.dim(3));
  sCheck = (roiBaseCoords(3,:) >= 1) & (roiBaseCoords(3,:) <= hdr.dim(4));

  % only use ones that are in bounds
  roiBaseCoords = roiBaseCoords(:,xCheck & yCheck & sCheck);

  % convert to linear coordinates
  roiBaseCoordsLinear = mrSub2ind(hdr.dim(2:4)',roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));

  if ~isempty(baseCoordMap) && baseType==1  %for flats, use basecoordmap to transform ROI from canonical base to multi-depth flat map
    roiBaseCoordsLinear = ismember(baseCoordsLinear,roiBaseCoordsLinear);
  end
  
  % set all the roi coordinates to 1
  d(roiBaseCoordsLinear) = 1;

  % now save the nifti file
  cbiWriteNifti(saveFilename{iRoi},d,hdr);
end
  
  
  
