% spatialSmooth.m: spatially smooth overlay with a 3D gaussian of given FWHM (in voxels)
%
%        $Id: spatialSmooth.m 2733 2013-05-13 11:47:54Z julien $
%           

function smoothed = spatialSmooth(overlay,FWHM)

smoothed = nanconvn(overlay,gaussianKernel(FWHM),'same');

function kernel = gaussianKernel(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
kernelDims = 2*[w w w]+1;
kernelCenter = ceil(kernelDims/2);
[X,Y,Z] = meshgrid(1:kernelDims(1),1:kernelDims(2),1:kernelDims(3));
kernel = exp(-((X-kernelCenter(1)).^2+(Y-kernelCenter(2)).^2+(Z-kernelCenter(3)).^2)/(2*sigma_d^2)); %Gaussian function
kernel = kernel./sum(kernel(:));

