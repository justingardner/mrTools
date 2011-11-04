% mlrExportROI.m
%
%        $Id:$ 
%      usage: mlrExportROI(v,saveFilename)
%         by: justin gardner
%       date: 07/14/09
%    purpose: Export an ROI to a nifti image
%
function mlrExportROI(v,saveFilename)

% check arguments
if ~any(nargin == [2])
  help mlrExportROI
  return
end

% get the roi we are being asked to export
roiNum = viewGet(getMLRView,'currentroi');
if isempty(roiNum)
  mrWarnDlg('(mlrExportROI) No current ROI to export');
  return
end
  

% get the base nifti header
hdr = viewGet(v,'basehdr');
if isempty(hdr)
  mrWarnDlg('(mlrExportROI) Could not get base anatomy header');
  return
end

% tell the user what is going on
disp(sprintf('(mlrExportROI) Exporting ROI to %s with dimensions set to match base %s: [%i %i %i]',saveFilename,viewGet(v,'baseName'),hdr.dim(2),hdr.dim(3),hdr.dim(4)));

% create a data structure that has all 0's
d = zeros(hdr.dim(2:4)');

% get  roi coordinates in base coordinates
roiBaseCoords = getROICoordinates(v,roiNum,0);

% check roiBaseCoords
if isempty(roiBaseCoords)
  mrWarnDlg('(mlrExportROI) This ROI does not have any coordinates in the base');
end

% make sure we are inside the base dimensions
xCheck = (roiBaseCoords(1,:) >= 1) & (roiBaseCoords(1,:) <= hdr.dim(2));
yCheck = (roiBaseCoords(2,:) >= 1) & (roiBaseCoords(2,:) <= hdr.dim(3));
sCheck = (roiBaseCoords(3,:) >= 1) & (roiBaseCoords(3,:) <= hdr.dim(4));

% only use ones that are in bounds
roiBaseCoords = roiBaseCoords(:,find(xCheck & yCheck & sCheck));

% convert to linear coordinates
roiBaseCoordsLinear = sub2ind(hdr.dim(2:4)',roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));

% set all the roi coordinates to 1
d(roiBaseCoordsLinear) = 1;

% now save the nifti file
cbiWriteNifti(saveFilename,d,hdr);

  
  
  
