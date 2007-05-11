function [img] = refreshMLRDisplay(viewNum)
%	$Id$

mrGlobals

% Get current view and baseNum.
% Get interp preferences.
% Get slice, scan, alpha, rotate, and sliceIndex from the gui.
%disppercent(-inf,'viewGet');
interpMethod = mrGetPref('interpMethod');
if isempty(interpMethod)
    interpMethod = 'linear';
end
interpExtrapVal = NaN;
view = viewGet(viewNum,'view');
scan = viewGet(view,'curscan');
slice = viewGet(view,'curslice');
alpha = viewGet(view,'alpha');
rotate = viewGet(view,'rotate');
baseNum = viewGet(view,'currentBase');
sliceIndex = viewGet(view,'baseSliceIndex',baseNum);
%disppercent(inf);

% Compute base coordinates and extract baseIm for the current slice
%disppercent(-inf,'extract base image');
[baseIm,baseCoords,baseCoordsHomogeneous] = ...
    getBaseSlice(view,slice,sliceIndex,rotate,baseNum);
view = viewSet(view,'cursliceBaseCoords',baseCoords);
baseDims = size(baseIm);
%disppercent(inf);

% Rescale base volume
%disppercent(-inf,'rescale base');
if ~isempty(baseIm)
    baseCmap = gray(256);
    baseClip = viewGet(view,'baseClip',baseNum);
    baseRGB = rescale2rgb(baseIm,baseCmap,baseClip);
else
    baseRGB = [];
end
%disppercent(inf);

% Extract overlay images and overlay coords, and alphaMap
%disppercent(-inf,'extract overlay images');
curOverlay = viewGet(view,'currentOverlay');
analysisNum = viewGet(view,'currentAnalysis');
[overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(view,scan,slice,sliceIndex,rotate,...
    baseNum,baseCoordsHomogeneous,baseDims,...
    analysisNum,interpMethod,interpExtrapVal);
view = viewSet(view,'cursliceOverlayCoords',overlayCoords);
if ~isempty(overlayImages)
    overlayIm = overlayImages(:,:,curOverlay);
else
    overlayIm = [];
end
%disppercent(inf);

%disppercent(-inf,'alphaMap');
numOverlays = viewGet(view,'numberofOverlays');
mask = ones(size(baseIm));
% Loop through overlays, filling in NaNs according to clip values.
if ~isempty(overlayImages)
    for ov = 1:numOverlays
        im = overlayImages(:,:,ov);
        clip = viewGet(view,'overlayClip',ov);
        % Find pixels that are within clip
        if diff(clip) > 0
            pts = (im >= clip(1) & im <= clip(2));
        else
            pts = (im >= clip(1) | im <= clip(2));
        end
        mask = mask & pts;
    end
end
% Finally, make the alphaMap.
alphaMap = repmat(alpha*mask,[1 1 3]);
%disppercent(inf);

% Rescale current overlay.
%disppercent(-inf,'rescale overlay');
if ~isempty(overlayIm)
    overlayCmap = viewGet(view,'overlayCmap',curOverlay);
    overlayRange = viewGet(view,'overlayRange',curOverlay);
    if strcmp(viewGet(view,'overlayCtype',curOverlay),'setRangeToMax')
        clip = viewGet(view,'overlayClip',curOverlay);
        overlayRange(1) = clip(1);
        overlayRange(2) = min(max(overlayIm(:)),clip(2));
    end
    overlayRGB = rescale2rgb(overlayIm,overlayCmap,overlayRange);
else
    overlayRGB = [];
end
%disppercent(inf);

% figure
% image(overlayRGB)
% colormap(overlayCmap)
% return

% Combine base and overlay
%disppercent(-inf,'combine base and overlay');
if ~isempty(baseRGB) & ~isempty(overlayRGB)
    img = (1-alphaMap).*baseRGB + alphaMap.*overlayRGB;
    cmap = overlayCmap;
    cbarRange = overlayRange;
elseif ~isempty(baseRGB)
    img = baseRGB;
    cmap = baseCmap;
    cbarRange = baseClip;
else
    % If no image at this point then display blank
    img = 0;
    cmap = gray(1);
    cbarRange = [0 1];
end
%disppercent(inf);

% If no image at this point then return
if ieNotDefined('img')
    return
end

% Display the image
%disppercent(-inf,'displayImage');
fig = viewGet(view,'figNum');
gui = guidata(fig);
% set(fig,'CurrentAxes',gui.axis);
image(img,'Parent',gui.axis);
axis(gui.axis,'off');
axis(gui.axis,'image');
%disppercent(inf);

% Display colorbar
%disppercent(-inf,'colorbar');
% set(fig,'CurrentAxes',gui.colorbar);
cbar = rescale2rgb([1:256],cmap,[1,256]);
image(cbar,'Parent',gui.colorbar);
set(gui.colorbar,'YTick',[]);
set(gui.colorbar,'XTick',[1 64 128 192 256]);
set(gui.colorbar,'XTicklabel',num2str(linspace(cbarRange(1),cbarRange(2),5)',3));
%disppercent(inf);

% Display the ROIs
%disppercent(-inf,'displayROI');
displayROIs(view,slice,sliceIndex,rotate,...
    baseNum,baseCoordsHomogeneous,baseDims);
%disppercent(inf);

return


%-------------------------------------------------------------------------
function rgb = rescale2rgb(image,cmap,clip)
%
% function rgb = rescale2rgb(image,cmap,[clipMin,clipMax])
%
% Clips top and bottom tails of image histogram.
% Sets NaNs to the lowest index in the colormap

if ~exist('clip','var')
    % Choose clipping based on histogram
    histThresh = length(image(:))/1000;
    [cnt, val] = hist(image(:),100);
    goodVals = find(cnt>histThresh);
    clipMin = val(min(goodVals));
    clipMax = val(max(goodVals));
    clip = [clipMin,clipMax];
else
    clipMin = clip(1);
    clipMax = clip(2);
end

% Clip
result = image;
result(find(image < clipMin)) = clipMin;
result(find(image > clipMax)) = clipMax;

% Scale
indices = round(255 * (result-clipMin)/(clipMax-clipMin)) + 1;
indices = max(1,min(indices,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(image));
g = zeros(size(image));
b = zeros(size(image));
r(:) = cmap(indices,1);
g(:) = cmap(indices,2);
b(:) = cmap(indices,3);

% Stuff them into rgb
dims = [size(image),3];
rgb = cat(length(dims),r,g,b);


%-------------------------------------------------------------------------
function [baseIm,baseCoords,baseCoordsHomogeneous] = ...
    getBaseSlice(view,sliceNum,sliceIndex,rotate,baseNum)
%
% getBaseSlice: extracts base image and corresponding coordinates

baseCoords = [];
baseCoordsHomogeneous = [];
baseIm = [];

% viewGet
volSize = viewGet(view,'baseDims',baseNum);
baseData = viewGet(view,'baseData',baseNum);

if ~isempty(volSize)

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

    % Rotate coordinates
    if (rotate ~= 0)
        x = imrotate(x,rotate,'nearest','loose');
        y = imrotate(y,rotate,'nearest','loose');
        z = imrotate(z,rotate,'nearest','loose');
    end

    % Reformat base coordinates
    imageDims = size(x);
    numPixels = prod(imageDims);
    xvec = reshape(x,1,numPixels);
    yvec = reshape(y,1,numPixels);
    zvec = reshape(z,1,numPixels);
    baseCoordsHomogeneous = [xvec; yvec; zvec; ones(1,numPixels)];
    baseCoords = reshape(baseCoordsHomogeneous(1:3,:)',[imageDims 3]);
end

% Extract base image
if ~isempty(baseData)
    switch sliceIndex
        case 1
            baseIm = squeeze(baseData(sliceNum,:,:));
        case 2
            baseIm = squeeze(baseData(:,sliceNum,:));
        case 3
            baseIm = squeeze(baseData(:,:,sliceNum));
    end
    if (rotate ~= 0)
        baseIm = imrotate(baseIm,rotate,'nearest','loose');
    end
end

return


%-------------------------------------------------------------------------
function [overlayImages,overlayCoords,overlayCoordsHomogeneous] = ...
    getOverlaySlice(view,...
    scanNum,sliceNum,sliceIndex,rotate,...
    baseNum,baseCoordsHomogeneous,imageDims,...
    analysisNum,interpMethod,interpExtrapVal);
%
% getOverlaySlice: extracts overlay image and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

% viewGet
baseXform = viewGet(view,'baseXform',baseNum);
overlayXform = viewGet(view,'overlayXform',scanNum);
numOverlays = viewGet(view,'numberofoverlays',analysisNum);
interpFnctn = viewGet(view,'overlayInterpFunction',analysisNum);

% Transform base coords corresponding to this slice/image to overlay
% coordinate frame.
% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin.
shiftXform = shiftOriginXform;
if ~isempty(overlayXform) & ~isempty(baseXform) & ~isempty(baseCoordsHomogeneous)
    xform = inv(shiftXform) * inv(overlayXform) * baseXform * shiftXform;
    % Transform coordinates
    overlayCoordsHomogeneous = xform * baseCoordsHomogeneous;
    overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[imageDims 3]);
end

% Extract overlayImages
if ~isempty(interpFnctn)
    overlayImages = feval(interpFnctn,view,scanNum,imageDims,...
        analysisNum,overlayCoords,interpMethod,interpExtrapVal);
else
    overlayImages = zeros([imageDims,numOverlays]);
    for ov = 1:numOverlays
        overlayData = viewGet(view,'overlayData',scanNum,ov,analysisNum);
        if ~isempty(overlayData) & ~isempty(overlayCoords)
            % Extract the slice
            overlayImages(:,:,ov) = interp3(overlayData,...
                overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
                interpMethod,interpExtrapVal);
        end
    end
end


%-------------------------------------------------------------------------
function overlayImages = corAnalInterp(view,scanNum,...
    imageDims,analysisNum,overlayCoords,interpMethod,interpExtrapVal);
%
% corAnalInterp: special case for corAnal. Need to treat amp/ph as complex
% valued.

% Initialize
numOverlays = viewGet(view,'numberofoverlays',analysisNum);
overlayImages = zeros([imageDims,numOverlays]);

% Interpolate complex values for amp and ph
ampNum = viewGet(view,'overlayNum','amp',analysisNum);
phNum = viewGet(view,'overlayNum','ph',analysisNum);
ampData = viewGet(view,'overlaydata',scanNum,ampNum,analysisNum);
phData = viewGet(view,'overlaydata',scanNum,phNum,analysisNum);
zData = ampData .* exp(i*phData);
if ~isempty(zData) & ~isempty(overlayCoords)
    zInterp = interp3(zData,...
        overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
        interpMethod,interpExtrapVal);
    overlayImages(:,:,ampNum) = abs(zInterp);
    overlayImages(:,:,phNum) = angle(zInterp);
end

% Loop through overlays and extract images, using ZInterp for amp and ph
for ov = 1:numOverlays
    overlayData = viewGet(view,'overlayData',scanNum,ov,analysisNum);
    if ~isempty(overlayData) & ~isempty(overlayCoords)
        % Extract the slice
        overlayImages(:,:,ov) = interp3(overlayData,...
            overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3),...
            interpMethod,interpExtrapVal);
    end
end



function [x,y] = getROISlice(view,sliceNum,sliceIndex,rotate,...
    baseNum,baseCoordsHomogeneous,imageDims,roiNum);
%
% getROISlice: extracts ROI coords transformed to the image

x = [];
y = [];

% viewGet
baseXform = viewGet(view,'baseXform',baseNum);
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
roiCoords = viewGet(view,'roiCoords',roiNum);
roiXform = viewGet(view,'roiXform',roiNum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);

if ~isempty(roiCoords) & ~isempty(roiXform) & ~isempty(baseXform)
    % Use xformROI to supersample the coordinates
    tmpCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));
    [roiCoordsSlice,roiIndices,baseIndices] = intersect(tmpCoords',baseCoordsHomogeneous','rows');
    % Get corresponding image slice coordinates
    [x,y] = ind2sub(imageDims,baseIndices);
end


%-------------------------------------------------------------------------
function displayROIs(view,sliceNum,sliceIndex,rotate,...
    baseNum,baseCoordsHomogeneous,imageDims);
%
% displayROIs: draws the ROIs in the current slice.
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
            case {'yellow','y'}, color = [1 1 0];
            case {'magenta','m'}, color = [1 0 1];
            case {'cyan','c'}, color = [0 1 1];
            case {'red','r'}, color = [1 0 0];
            case {'green','g'}, color = [0 1 0];
            case {'blue','b'}, color = [0 0 1];
            case {'white','w'}, color = [1 1 1];
            case {'black','k'}, color = [0 0 0];
            otherwise, color = [1 1 1];
        end % end switch statement
    end % end loop

    % Get ROI coords transformed to the image
    [x,y] = getROISlice(view,sliceNum,sliceIndex,rotate,...
        baseNum,baseCoordsHomogeneous,imageDims,r);
    if ~isempty(x) & ~isempty(y)

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


