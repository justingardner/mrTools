function result = warpAffineInterp3(tseries,frame,A,sliceTimes,badVal,method)
%
% function result = warpAffineInterp3(tseries,frame,A,[sliceTimes],[badVal],[method])
%
% tseries is a time series of volumes (4d array)
% frame is the frame number of the volume that is being warped and
%    interpolated.
% A: 3x4 affine transform matrix or a 4x4 matrix with [0 0 0 1]
%    for the last row.
% badVal: if a transformed point is outside of the volume, badVal is used
%    (default = 0).
% method:  'nearest', 'linear', 'cubic', or 'spline' (default = 'linear')
%
% result: output volume, same size as one frame of tseries
%
% 10/99 ON: added Border parameter to permit nearest neighbor interpolation at the edges
% 7/2007 DJH: modified from warpAffine3 to do slice time correction
%    (interplation in time) as well as warping in space. This involved
%    changing the input parameters to be a time series instead of a single
%    frame.

if ieNotDefined('sliceTimes')
  sliceTimes = zeros(1,size(tseries,3));
end
if ieNotDefined('badVal')
  badVal = 0;
end
if ieNotDefined('method')
  method = 'linear';
end

if (size(A,1)>3)
  A = A(1:3,:);
end

% original size
[Ny Nx Nz Nt] = size(tseries);

% Compute coordinates corresponding to input volume
% and transformed coordinates for result
[ygrid,xgrid,zgrid,tgrid] = ndgrid(1:Ny,1:Nx,1:Nz,frame);
coords = [xgrid(:)'; ygrid(:)'; zgrid(:)'];
homogeneousCoords = [coords; ones(1,size(coords,2))];
warpedCoords = A*homogeneousCoords;
% Slice time correction. Interpolated to the end of each TR. Would have to
% add 1/2 to interpolate to the middle of each TR but that busts the last
% frame of the tseries.
for slice = 1:size(tgrid,3)
  tgrid(:,:,slice) = tgrid(:,:,slice) - sliceTimes(slice);
end

% Interpolate
result = interpn(tseries,...
  warpedCoords(2,:),warpedCoords(1,:),warpedCoords(3,:),tgrid(:)',...
  method, badVal);

% Reshape result
result = reshape(result,[Ny Nx Nz]);

return;

%%% Debug

slice=[0 0 0; 0.1 0.1 0.1; 0.2 0.2 0.2]';
tseries=ones(3,3,4,2);
for z=1:4
  tseries(:,:,z,1)=slice+z;
end
tseries(:,:,:,2)=tseries(:,:,:,1)+4;

A = eye(4);
sliceTimes = [0 .25 .5 .75];

res=warpAffineInterp3(tseries,2,A,sliceTimes,NaN)
res=warpAffineInterp3(tseries,1,A,sliceTimes,NaN)
res=warpAffineInterp3(tseries,1,A,sliceTimes,NaN,'nearest')
