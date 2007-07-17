function[] = mrExport2SR(viewNum, pathstr)
% mrExport2SR.m
%
%      usage: [] = mrExprt2SR(viewNum, pathstr)
%         by: eli merriam
%       date: 03/20/07
%    purpose: exports a MLR overlay to a Nifti file compatible with SurfRelax
%        $Id$	
%

mrGlobals

% interpMethod and interpExtrapVal are used by calls to interp3 when
% extracting slices from the base and overlay arrays.
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
    interpMethod = 'linear';
end
interpExtrapVal = NaN;

% Get view
view = viewGet(viewNum,'view');

% Get values from the GUI
scan = viewGet(view,'curscan');

% sliceIndex depends both on sliceOrientation and on the orientation
% of baseVolume.
baseNum = viewGet(view,'currentBase');

rotate = 0;
sliceIndex = 1;

basedims = viewGet(view, 'basedims');
overlayIm = zeros(basedims);
overlayNum = viewGet(view,'currentOverlay');
overlayData = viewGet(view,'overlayData',scan,overlayNum);

disppercent(-inf,sprintf('Exporting resampled Nifti file: %s', pathstr));
for baseSlice = 1:basedims(1)
    % Compute base and overlay coordinates for the current slice
    [baseCoords,overlayCoords] = getSliceCoords(view,scan,baseSlice,sliceIndex,rotate);
    view = viewSet(view,'cursliceBaseCoords',baseCoords);
    view = viewSet(view,'cursliceOverlayCoords',overlayCoords);
    
    % Extract slice from current overlay.
    if ~isempty(overlayNum) & ~isempty(overlayCoords) & ~isempty(overlayData)
        overlayIm(baseSlice,:,:) = interp3(overlayData,overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3), interpMethod,interpExtrapVal);
    end
    disppercent(baseSlice/basedims(1));
end

disppercent(inf)
baseVolume = viewGet(viewNum,'baseVolume');
hdr = baseVolume.hdr;
hdr.bitpix = 32;   
hdr.datatype = 16;
hdr.is_analyze = 1;
hdr.scl_slope = 1;
hdr.endian = 'l';
hdr.pixdim = [1 1 1 1 0 0 0 0]';        % all pix dims must be specified here
                                        %cbiWriteNifti(pathstr, overlayIm, hdr);
cbiWriteNifti(pathstr, overlayIm);

return


function [baseCoords,overlayCoords] = getSliceCoords(view,scanNum,sliceNum,sliceIndex,rotate)

baseCoords = [];
overlayCoords = [];

% Use baseVolume.xform to transform the overlay data into the
% selected slice.
baseNum = viewGet(view,'currentBase');
if baseNum
    volSize = viewGet(view,'baseDims',baseNum);
    baseXform = viewGet(view,'baseXform',baseNum);
else
    return
end

% Generate coordinates with meshgrid
x = sliceNum * ones(volSize(2),volSize(3));
[z,y] = meshgrid(1:volSize(3),1:volSize(2));

% Reformat coordinates
sliceDims = size(x);
numPixels = prod(sliceDims);
xvec = reshape(x,1,numPixels);
yvec = reshape(y,1,numPixels);
zvec = reshape(z,1,numPixels);
baseCoordsHomogeneous = [xvec; yvec; zvec; ones(1,numPixels)];
baseCoords = reshape(baseCoordsHomogeneous(1:3,:)',[sliceDims 3]);

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;

% Transform to overlay coordinates
overlayNum = viewGet(view,'currentOverlay');
if overlayNum
    overlayXform = viewGet(view,'overlayXform',scanNum);
    if ~isempty(overlayXform)
        xform = inv(shiftXform) * inv(overlayXform) * baseXform * shiftXform;
        % Transform coordinates
        overlayCoordsHomogeneous = xform * baseCoordsHomogeneous;
        overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[sliceDims 3]);
    end
end




