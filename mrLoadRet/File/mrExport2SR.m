function[] = mrExport2SR(viewNum, pathstr)
% mrExport2SR.m
%
%      usage: [] = mrExprt2SR(viewNum, pathstr)
%         by: eli merriam
%       date: 03/20/07
%    purpose: exports a MLR overlay to a Nifti file compatible with SurfRelax
%	$Id$	
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
alpha = viewGet(view,'alpha');
%rotate = viewGet(view,'rotate');

% sliceIndex depends both on sliceOrientation and on the orientation
% of baseVolume.
baseNum = viewGet(view,'currentBase');
sliceIndex = viewGet(view,'baseSliceIndex',baseNum);

% jg: don't think that this code is handling rotate
% correctly. when the overlay gets stored in the overlayIn
% if rotate is sent to something other than 0 then the dimensions
% don't match. also, why should this code even bother to save
% the image rotated? why not always save in the same way that
% the base anatomy is stored since that is probably what is expeted.
% so, always set sliceIndex to 1 and rotate = 0 no matter what
% the user settings are...
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
hdr.datatype = 4;                       % SurfaceViewer requires shorts
hdr.pixdim = [1 1 1 1 1 1 1 1]';        % all pix dims must be specified here
cbiWriteNifti(pathstr, overlayIm, hdr)

return


function [baseCoords,overlayCoords] = getSliceCoords(view,scanNum,sliceNum,sliceIndex,rotate)

baseCoords = [];
overlayCoords = [];

switch view.viewType

  case 'Volume'
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
    switch sliceIndex
      case 1
        x = sliceNum * ones(volSize(2),volSize(3));
        [z,y] = meshgrid(1:volSize(3),1:volSize(2));
      case 2
        y = sliceNum * ones(volSize(1),volSize(3));
        [z,x] = meshgrid(1:volSize(3),1:volSize(1));
      case 3
        z = sliceNum * ones(volSize(1),volSize(2));
        [y,x] = meshgrid(1:volSize(2),1:volSize(1));
    end

    % Rotate
    if (rotate ~= 0)
        xrot = imrotate(x,rotate,'nearest','loose');
        yrot = imrotate(y,rotate,'nearest','loose');
        zrot = imrotate(z,rotate,'nearest','loose');
    else
        xrot = x;
        yrot = y;
        zrot = z;
    end

    % Reformat coordinates
    sliceDims = size(xrot);
    numPixels = prod(sliceDims);
    xvec = reshape(xrot,1,numPixels);
    yvec = reshape(yrot,1,numPixels);
    zvec = reshape(zrot,1,numPixels);
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

  case {'Surface'}
    mrErrorDlg('Surface view not implemented yet');

  case 'Flat'
    mrErrorDlg('Flat view not implemented yet');

end


function displayROIs(view)
%
% displayROIs(view)
%
% Draws the ROIs in the current slice.
%
% viewGet(view,'showROIs') can be:
% 'all'
% 'selected'
% 'all perimeter'
% 'selected perimeter'
% 'hide'

s = viewGet(view,'currentroi');
n = viewGet(view,'numberOfROIs');

% Order in which to draw the ROIs
option = viewGet(view,'showROIs');
switch option
  case{'all','all perimeter'}
    if s
        order = [1:s-1,s+1:n,s];
    else
        order = 1:n;
    end
  case{'selected','selected perimeter'}
    order = s;
  otherwise
    return
end

% Get baseCoords corresponding to the current slice, outside loop through
% ROIs.
baseNum = viewGet(view,'currentBase');
baseXform = viewGet(view,'baseXform',baseNum);
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
baseCoords = viewGet(view,'curSliceBaseCoords');
sliceDims = [size(baseCoords,1),size(baseCoords,2)];
numPixels = prod(sliceDims);
baseX = baseCoords(:,:,1);
baseY = baseCoords(:,:,2);
baseZ = baseCoords(:,:,3);
baseCoordsHomogeneous = [baseX(:)'; baseY(:)'; baseZ(:)'; ones(1,numPixels)];

% Loop through ROIs in order
for r = order
    if (r == s)
        % Selected ROI: set color=white
        color = [1 1 1];
    else
        % Non-selected ROI, get color
        thisCol = viewGet(view,'roiColor',r);
        % If it's a 'text' color, translate it...
        switch (thisCol)
          case 'yellow', color = [1 1 0];
          case 'magenta', color = [1 0 1];
          case 'cyan', color = [0 1 1];
          case 'red', color = [1 0 0];
          case 'green', color = [0 1 0];
          case 'blue', color = [0 0 1];
          case 'white', color = [1 1 1];
          case 'black', color = [0 0 0];
          otherwise, color = [1 1 1];
        end % end switch statement
    end % end loop

    % Get subset of coords corresponding to the current slice
    roiCoords = viewGet(view,'roiCoords',r);
    roiXform = viewGet(view,'roiXform',r);
    roiVoxelSize = viewGet(view,'roiVoxelSize',r);
    
    if ~isempty(roiCoords)

        % Use xformROI to supersample the coordinates
        tmpCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));
        [roiCoords,roiIndices,baseIndices] = intersect(tmpCoords',baseCoordsHomogeneous','rows');

        % Get corresponding image slice coordinates
        [x,y] = ind2sub(size(baseX),baseIndices);

        % Draw it
        fig = viewGet(view,'figNum');
        gui = guidata(fig);
        set(fig,'CurrentAxes',gui.axis);
        lineWidth = 0.5;
        w = 0.5;
        hold on
        switch option

          case{'all','selected'}
            % Draw the lines around each pixel, w=1/2 because you don't
            % want to connect the centers of the pixels, rather you want to
            % draw around each pixel, e.g, from (x-.5,y-.5) to (x+.5,y-.5).
            for i=1:length(x);
                line([y(i)-w,y(i)+w,y(i)+w,y(i)-w,y(i)-w],...
                     [x(i)-w,x(i)-w,x(i)+w,x(i)+w,x(i)-w], ...
                     'Color',color,'LineWidth',lineWidth);
            end

          case{'all perimeter','selected perimeter'}
            % Draw only the perimeter
            for i=1:length(x);
                xMinus = find(x == x(i)-1);
                xEquals = find(x == x(i));
                xPlus = find(x == x(i)+1);
                if isempty(xMinus)
                    line([y(i)-w,y(i)+w],[x(i)-w, x(i)-w],'Color',color,'LineWidth',lineWidth);
                else
                    if ~any(y(i) == y(xMinus))
                        line([y(i)-w,y(i)+w],[x(i)-w, x(i)-w],'Color',color,'LineWidth',lineWidth);
                    end
                end
                if isempty(xPlus)
                    line([y(i)-w,y(i)+w],[x(i)+w, x(i)+w],'Color',color,'LineWidth',lineWidth);
                else
                    if ~any(y(i) == y(xPlus))
                        line([y(i)-w,y(i)+w],[x(i)+w, x(i)+w],'Color',color,'LineWidth',lineWidth);
                    end
                end
                if ~isempty(xEquals)
                    if ~any(y(i) == y(xEquals)-1)
                        line([y(i)+w,y(i)+w],[x(i)-w, x(i)+w],'Color',color,'LineWidth',lineWidth);
                    end
                    if ~any(find(y(i) == y(xEquals)+1))
                        line([y(i)-w,y(i)-w],[x(i)-w, x(i)+w],'Color',color,'LineWidth',lineWidth);
                    end
                end
            end
        end
    end
    hold off

end

return;


