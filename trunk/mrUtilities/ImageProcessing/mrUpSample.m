function out = mrUpSample(im, nLevels, filt)
% out = mrUpSample(im, nLevels, [filt])
%
% Upsamples the image im by the integer nLevels.
% The upsampled image is blurred and edges are 
% dealt with by zero-padding.
%
% The blurring is done with the filter kernel specified by filt
% (default = [.125 .5 .75 .5 .125]'), which should be a vector 
% (applied separably as a 1D convolution kernel in X and Y), 
% or a matrix (applied as a 2D convolution kernel)
%
% 99.08.16 RFD wrote it, based on upBlur (and helpful
%				comments from DJH)
% 07.08.19, EPM added upsampling across slices
% $Id$	

if nLevels ~= fix(nLevels)
    error('nLevels must be an integer!!!');
end

if ~exist('filt', 'var')
    % default filter for post-upsample convolution
    filt = sqrt(2)*namedFilter('binom5');
end

% use a little recursion to deal with upsample steps > 2
if nLevels > 1
    im = mrUpSample(im, nLevels-1);
end
if (nLevels >= 1)
    if (any(size(im)==1))
        if (size(im,1)==1)
            filt = filt';
        end
        out = upConv(im, filt, 'zero',(size(im)~=1)+1);
    else
        out = im;
        [nrows ncols nslices nframes] = size(out);

        disppercent(-inf,'(mrUpSample) upsampling the rows');
        out = reshape(out,[nrows ncols*nslices*nframes]);
        out = upConv(out, filt, 'zero', [2 1]);
        nrows = nrows*2;
        out = reshape(out,[nrows ncols nslices nframes]);
        disppercent(inf);

        disppercent(-inf,'(mrUpSample) upsampling the columns');
        out = permute(out,[2 1 3 4]);
        out = reshape(out,[ncols nrows*nslices*nframes]);
        out = upConv(out, filt, 'zero', [2 1]);
        ncols = ncols*2;
        out = reshape(out,[ncols nrows nslices nframes]);
        out = permute(out,[2 1 3 4]);
        disppercent(inf);
        
        disppercent(-inf,'(mrUpSample) upsampling the slices');
        out = permute(out,[3 2 1 4]);
        out = reshape(out,[nslices ncols*nrows*nframes]);
        out = upConv(out, filt, 'zero', [2 1]);
        nslices = nslices*2;
        out = reshape(out,[nslices ncols nrows nframes]);
        out = permute(out,[3 2 1 4]);
        disppercent(inf);
    end
else
    out = im;
end

return

