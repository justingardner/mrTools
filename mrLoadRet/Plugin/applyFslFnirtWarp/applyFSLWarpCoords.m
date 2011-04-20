% warpedCoords = applyFSLWarpCoords(coords,voxelSize, warpCoefFilename, tempFilename, hdr, verbose)
%
%   applies non-linear FSL registration to coordinates
%
%  input: 
%        coords: coordinates in the space of the input FNIRT volume, in homogeneous format (4*n matrix, with 1s on the last row)
%     voxelSize: resolution in mm of the warp fiel computed using fnirtfileutils
%  tempFilename: name of the temporary file to write the data to the disc
%           hdr: header of a volume in the space of the input FNIRT volume for dimensions and voxel size (transformation matrix is not taken into account)
%
% jb 23/07/2010
%
% $Id: applyFSLWarp.m 2107 2011-04-17 19:49:52Z julien $

function warpedCoords = applyFSLWarpCoords(coords,voxelSize, warpCoefFilename, tempFilename, hdr, verbose)


if ieNotDefined('verbose')
   verbose = true;
end

scalingFactor = round(hdr.pixdim(2:4)'./voxelSize); 
dataSize = hdr.dim(2:4)'.*scalingFactor;
data = ones(dataSize,'single');
hdr.pixdim(2:4) = hdr.pixdim(2:4)./scalingFactor';
hdr.sform44 = hdr.sform44*diag([1./scalingFactor 1]);
hdr.qform44 = hdr.qform44*diag([1./scalingFactor 1]);

cbiWriteNifti(tempFilename, data, hdr,[],[],[],verbose);
clear data

try
command =  sprintf('fnirtfileutils --in=%s --ref=%s --out=%s  --withaff', warpCoefFilename, tempFilename, tempFilename);
  if verbose
    fprintf('(applyFSLWarpCoords) Computing FNIRT warp fields at a resolution of %s mm:\n',mat2str(hdr.pixdim(2:4)));
    disp(['  ' command])
  end;
  [s,w] = unix(command);
  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(applyFSLWarp) There was a problem running FSL applywarp command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the warped fields
warpFields = cbiReadNifti(tempFilename);
warpFields = permute(warpFields,[4 1 2 3]);

coordsIndex = round(repmat([scalingFactor 1]',1,size(coords,2)).*(coords-.5) + .5);
coordsIndex = mrSub2ind(dataSize, coordsIndex(1,:)',coordsIndex(2,:)',coordsIndex(3,:)');
isInScan = ~isnan(coordsIndex); %ignore voxels outside the input volume

warpFields = warpFields(:,coordsIndex(isInScan));
%for some reason, x an y axis are inverted in the warp fields...)
warpFields = [-1 0 0;0 -1 0;0 0 -1]*warpFields;
warpedCoords = NaN(size(coords));
warpedCoords(:,isInScan) = [warpFields;zeros(1,size(warpFields,2))] + coords(:,isInScan);
