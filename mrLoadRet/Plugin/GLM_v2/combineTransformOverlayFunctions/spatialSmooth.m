%
% function smoothed = spatialSmooth(overlay,FWHM)
%
% spatialSmooth.m: spatially smooth overlay with a 3D gaussian of given FWHM (in voxels)
%
%        $Id: spatialSmooth.m 2733 2013-05-13 11:47:54Z julien $
%           

function smoothed = spatialSmooth(overlay,FWHM)

smoothed = nanconvn(overlay,gaussianKernel(FWHM),'same');
smoothed(isnan(overlay))=NaN;

function kernel = gaussianKernel(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
if length(w)==1
  w = [w w w];
  sigma_d = [sigma_d sigma_d sigma_d];
end
kernelDims = 2*w+1;
kernelCenter = ceil(kernelDims/2);
[X,Y,Z] = meshgrid(1:kernelDims(1),1:kernelDims(2),1:kernelDims(3));
if sigma_d(2) == 0 && sigma_d(3) == 0
  kernel = exp(-((X-kernelCenter(1)).^2/(2*sigma_d(1)^2))); % 1D Gaussian function along X
elseif sigma_d(1) == 0 && sigma_d(3) == 0
  kernel = exp(-((Y-kernelCenter(2)).^2/(2*sigma_d(2)^2))); % 1D Gaussian function along Y
elseif sigma_d(1) == 0 && sigma_d(2) == 0
  kernel = exp(-((Z-kernelCenter(3)).^2/(2*sigma_d(3)^2))); % 1D Gaussian function along Z
elseif sigma_d(3) == 0
  kernel = exp(-((X-kernelCenter(1)).^2/(2*sigma_d(1)^2)+(Y-kernelCenter(2)).^2/(2*sigma_d(2)^2))); % 2D Gaussian function in XY plane
elseif sigma_d(2) == 0
  kernel = exp(-((X-kernelCenter(1)).^2/(2*sigma_d(1)^2)+(Z-kernelCenter(3)).^2/(2*sigma_d(3)^2))); % 2D Gaussian function in XZ plane
elseif sigma_d(1) == 0
  kernel = exp(-((Y-kernelCenter(2)).^2/(2*sigma_d(2)^2)+(Z-kernelCenter(3)).^2/(2*sigma_d(3)^2))); % 2D Gaussian function in YZ plane
else
  kernel = exp(-((X-kernelCenter(1)).^2/(2*sigma_d(1)^2)+(Y-kernelCenter(2)).^2/(2*sigma_d(2)^2)+(Z-kernelCenter(3)).^2/(2*sigma_d(3)^2))); % 3D Gaussian function
end
kernel = kernel./sum(kernel(:));

