% mlrExportROI.m
%
%        $Id$ 
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

% create a data structure that has all 0's
d = zeros(hdr.dim(2:4)');

% get  roi coordinates in base coordinates
roiBaseCoords = getROICoordinates(v,roiNum,0);

% check roiBaseCoords
if isempty(roiBaseCoords)
  mrWarnDlg('(mlrExportROI) This ROI does not have any coordinates in the base');
end

% convert to linear coordinates
roiBaseCoordsLinear = sub2ind(hdr.dim(2:4)',roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));

% set all the roi coordinates to 1
d(roiBaseCoordsLinear) = 1;

% now save the nifti file
cbiWriteNifti(saveFilename,d,hdr);

  
  
  
