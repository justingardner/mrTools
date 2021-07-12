% mlrExportROI.m
%
%        $Id$ 
%      usage: mlrExportROI(v,saveFilename,<'baseNum',baseNum>,<'scanNum',scanNum>,<'groupNum',groupNum>,<'hdr',hdr>,<'exportToFreesurferLabel',true/false>)
%         by: justin gardner
%       date: 07/14/09
%    purpose: Export ROI(s) to a nifti image or Freesurfer label file. Uses
%             current roi(s) and current base in view to export. To use a different
%             base, specify 'baseNum'. To use a scan, specify 'scanNum' and 'groupNum'.
%             Pass in a nifti header as hdr argument if you want to use a different header
%             'baseNum', 'scanNum', 'groupNum' and 'hdr' are ignored if exportToFreesurferLabel is true
%
function mlrExportROI(v,saveFilename,varargin)

% check arguments
if nargin < 2
  help mlrExportROI
  return
end

% optional arguments
getArgs(varargin,{'baseNum=[]','scanNum=[]','groupNum=[]','hdr=[]','exportToFreesurferLabel=0'});

% get the roi we are being asked to export
roiNum = viewGet(v,'currentroi');
if isempty(roiNum)
  mrWarnDlg('(mlrExportROI) No current ROI to export');
  return
end
  
if ischar(saveFilename)
  saveFilename = {saveFilename};
end
if ~isequal(length(roiNum),length(saveFilename))
  mrWarnDlg('(mlrExportROI) number of file names must be identical to number of ROIs');
  return
end

if ~isempty(scanNum) && ~isempty(groupNum)
  exportToScanSpace = true;
  exportVolumeString = 'scan';
  baseCoordMap = [];
  baseType = 0;
else
  exportToScanSpace = false;
  exportVolumeString = 'base';
  if isempty(baseNum)
    baseNum = viewGet(v,'curBase');
  end
  baseType = viewGet(v,'basetype',baseNum);
  baseCoordMap = viewGet(v,'basecoordmap',baseNum);
end

% get the base nifti header
passedInHeader = false;
if ~isempty(hdr)
  passedInHeader = true;
else
  if exportToScanSpace
    hdr = viewGet(v,'niftiHdr',scanNum,groupNum);
  else
    hdr = viewGet(v,'basehdr',baseNum);
  end
  if isempty(hdr)
    mrWarnDlg(sprintf('(mlrExportROI) Could not get %s header',exportVolumeString));
    return
  end
end

if ~ismember(baseType,[2]) && exportToFreesurferLabel
  mrWarnDlg('(mlrExportROI) Load a surface in order to export to freesurfer label format');
  % for surfaces, the required list of surface vertices is in baseCoordMap
  % for flat maps (or volumes), there is no easy access to the corresponding surface vertices, so returning
  return;
end
if ~isempty(baseCoordMap) && (baseType==1 || exportToFreesurferLabel) %for flats, or when exporting to freesurfer label file, use basecoordmap 
  
  if baseType == 1 && ~exportToFreesurferLabel 
      [~,baseCoords,baseCoordsHomogeneous] = getBaseSlice(v,1,3,viewGet(v,'baseRotate',baseNum),baseNum,baseType);
  else  
      baseCoords = permute(baseCoordMap.coords,[1 2 4 5 3]);
      baseCoordsHomogeneous = [permute(reshape(baseCoordMap.coords, ...
        [size(baseCoordMap.coords,1)*size(baseCoordMap.coords,2) size(baseCoordMap.coords,3) size(baseCoordMap.coords,4) size(baseCoordMap.coords,5)]),...
        [3 1 4 2]); ones(1,size(baseCoordMap.coords,1)*size(baseCoordMap.coords,2),size(baseCoordMap.coords,5))];
  end
    
  volDims = size(baseCoords);
  volDims = volDims ([1 2 4]);
  if baseType==1 && ~exportToFreesurferLabel
    mrWarnDlg(sprintf('(mlrExportROI) Exporting ROI(s) to flat space (%d x %d x %d voxels). If you do not want this, load base volume or surface',volDims(1),volDims(2),volDims(3)));
  elseif exportToFreesurferLabel
    mrWarnDlg('(mlrExportROI) Vertex coordinates in label file might be incorrect.');
  end
  % make sure that baseCoords are rounded (they may not be if we are working with a baseCoordMap's flat map)
  baseCoordsHomogeneous = reshape(baseCoordsHomogeneous,4,prod(volDims));
  baseCoordsHomogeneous = round(baseCoordsHomogeneous);
  baseCoordsLinear = mrSub2ind(baseCoordMap.dims,baseCoordsHomogeneous(1,:),baseCoordsHomogeneous(2,:),baseCoordsHomogeneous(3,:));

  if baseType==1 && ~exportToFreesurferLabel
    % estimate voxel size (taken from getBaseOverlay, assuming mask is invariant to rotation, which it should be since it is a flat map)
    oldBaseVoxelSize=viewGet(v,'basevoxelsize',baseNum);
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
  end
else
  volDims = hdr.dim(2:4)';
end

if ~passedInHeader && ~exportToFreesurferLabel
  b = viewGet(v,'base',baseNum);
end

for iRoi = 1:length(roiNum)
  roiName = viewGet(v,'roiName',roiNum(iRoi));
  % tell the user what is going on
  fprintf('(mlrExportROI) Exporting ROI %s to %s with dimensions set to match %s: [%i %i %i]\n',roiName, saveFilename{iRoi},exportVolumeString,volDims(1),volDims(2),volDims(3));

  % get  roi coordinates in base/scan coordinates
  if exportToScanSpace
    roiBaseCoords = getROICoordinates(v,roiNum(iRoi),scanNum,groupNum);
  else
    roiBaseCoords = getROICoordinates(v,roiNum(iRoi),0,[],sprintf('baseNum=%d',baseNum));
  end

  if ~exportToFreesurferLabel
    % create a data structure that has all 0's
    d = zeros(volDims);


    % make sure we are inside the base dimensions
    xCheck = (roiBaseCoords(1,:) >= 1) & (roiBaseCoords(1,:) <= hdr.dim(2));
    yCheck = (roiBaseCoords(2,:) >= 1) & (roiBaseCoords(2,:) <= hdr.dim(3));
    sCheck = (roiBaseCoords(3,:) >= 1) & (roiBaseCoords(3,:) <= hdr.dim(4));

    % only use ones that are in bounds
    roiBaseCoords = roiBaseCoords(:,xCheck & yCheck & sCheck);
  end
  

  if ~isempty(baseCoordMap) && (baseType==1 || exportToFreesurferLabel)  %for flats and surfaces, use basecoordmap to transform ROI from canonical base to multi-depth flat map
    roiBaseCoordsLinear = mrSub2ind(baseCoordMap.dims',roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
    roiBaseCoordsLinear = ismember(baseCoordsLinear,roiBaseCoordsLinear);
  else
    % convert to linear coordinates
    roiBaseCoordsLinear = mrSub2ind(hdr.dim(2:4)',roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
  end

  % check roiBaseCoords
  if isempty(roiBaseCoords) || ~nnz(roiBaseCoordsLinear)
    mrWarnDlg(sprintf('(mlrExportROI) This ROI (%s) does not have any coordinates in the %s',roiName,exportVolumeString));
    
  else
    
    if exportToFreesurferLabel
      % in order to export to label format, select vertices that are within the ROI

      % reshape to vertices * depths
      roiBaseCoordsLinear = reshape(roiBaseCoordsLinear, [size(baseCoordMap.coords,1)*size(baseCoordMap.coords,2) size(baseCoordMap.coords,5)]);
      % Multiple cortical depths are not taken into account: a vertex can be in the ROI at any depth:
      % but let's be conservative and consider only ROI voxels in the central part of cortical ribbon
      nDepths = size(roiBaseCoordsLinear,2);
      roiBaseCoordsLinear = find(any(roiBaseCoordsLinear(:,ceil((nDepths-1)/4)+1:floor(3*(nDepths-1)/4)+1),2));
      %actual coordinates in label file will be midway between inner and outer surface
      vertexCoords = (baseCoordMap.innerVtcs(roiBaseCoordsLinear,:)+baseCoordMap.outerVtcs(roiBaseCoordsLinear,:))/2;
      % change vertex coordinates to freesurfer system: 0 is in the middle of the volume and coordinates are in mm
      % (this does not seem to give the correct coordinates, but coordinates are usually not needed in label files)
      vertexCoords = (vertexCoords - repmat(baseCoordMap.dims/2 ,size(vertexCoords,1),1)) ./ repmat(hdr.pixdim([2 3 4])',size(vertexCoords,1),1);
      % (this assumes that the base volume for this surface is either the original Freesurfer volume, or has been cropped symmetrically, which is usually the case)
      
      % find freesurfer subject name
      if isempty(which('extractBetween'))
        freesurferName = viewGet(v,'subject'); % this is an old version of Matlab and I can't be bothered to find a replacement for extractBetween()
      else
        freesurferName = extractBetween(baseCoordMap.path,'subjects\','\surfRelax');
        if isempty(freesurferName)
          freesurferName= viewGet(v,'subject');
        elseif iscell(freesurferName)   %in newer versions of Matlab, extractBetween may return a cell array
          freesurferName = freesurferName{1};
        end
      end
      labelHeader = sprintf('#!ascii label, ROI exported from subject %s using mrTools (mrLoadRet v%.1f)\n', freesurferName, mrLoadRetVersion);
      % a Freesurfer label file is text file with a list of vertex numbers and coordinates
      mlrWriteFreesurferLabel(saveFilename{iRoi},labelHeader,[roiBaseCoordsLinear-1, vertexCoords, ones(size(vertexCoords,1),1)])
    else
      % set all the roi coordinates to 1
      d(roiBaseCoordsLinear) = 1;

      % if the orientation has been changed in loadAnat, undo that here.
      if ~passedInHeader && ~isempty(b.originalOrient)
        % convert into mlrImage
        [d, h] = mlrImageLoad(d,hdr);
        % convert the orientation back to original
        [d, h] = mlrImageOrient(b.originalOrient,d,h);
        % convert back to nifti
        reorientedHdr = mlrImageGetNiftiHeader(h);
        cbiWriteNifti(saveFilename{iRoi},d,reorientedHdr);
      else
        cbiWriteNifti(saveFilename{iRoi},d,hdr);
      end
    end
  end
end
  
  
  
