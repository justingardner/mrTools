function checkAlignment(correction)
% function checkAlignment(vol,inp,xform,correction)
%
% Interpolates the inplanes corresponding to the estimated rotation (rot)
% and translation (trans) and displays slice by slice the original
% inplanes, the interpolated inplanes and a mosaic of both.
%
% Oscar Nestares - 5/99
% djh 9/2003, modified into a function and cleaned up various things
% djh 1/2007, updated to mrAlign-4.5

global ALIGN

if ~exist('correction','var')
	RESP = questdlg('Correct intensity?', 'Correct', 'Yes', 'No', 'Yes');
	if strcmp(RESP, 'Yes')
		correction = 1;
	else
		correction = 0;
	end
end

% get alignment xform
xform = ALIGN.xform * ALIGN.guiXform;

hmsgbox = mrMsgBox('Please wait: interpolating slices, correcting intensity & contrast...');

% interpolate the inplanes
volInterp = interpVolume(ALIGN.volume, xform, ALIGN.inplaneSize, ALIGN.origin);
inp = ALIGN.inplanes;

% intensity & contrast correction
if correction
	inpIC = intensityContrastCorrection(inp, ALIGN.crop);
	volInterpIC = intensityContrastCorrection(volInterp, ALIGN.crop);
else
	inpIC = inp;
	volInterpIC = volInterp;
end

% make mosaic
mosaic = zeros(size(inpIC));
for k = 1:size(inp,3)
   mosaic(:,:,k) = imgMosaic(volInterpIC(:,:,k),inpIC(:,:,k));
end
    
mrCloseDlg(hmsgbox);

% open checkAlignmentGUI
checkAlignmentGUI(inpIC,mosaic,volInterpIC);


function im = imgMosaic(im1, im2, N)
% imgMosaic - generates a checkerboard mosaic from the two input images
%             with blocks of size NxN (default 10)
%
%       im = regMosaic(im1, im2, <N>)
%
% Oscar Nestares - 5/99

[Ny Nx] = size(im1);
if nargin<3
   N=round(max(Ny,Nx)/15);
end

im1(find(isnan(im1))) = 0;
im2(find(isnan(im2))) = 0;

basic = [ones(N) zeros(N); zeros(N) ones(N)];
check = repmat(basic, ceil(Ny/N), ceil(Nx/N));
check = check(1:Ny,1:Nx);

im = check.*im1 + (~check).*im2;

% puts a slightly different contrast in one image than in the other
% im = check.*im1*0.9 + (~check).*im2;

