function cropRegion = selectCropRegion(volume)
% function cropRegion = selectCropRegion(volume)

cropRegion = [];
if isempty(volume)
  mrWarnDlg('(selectCropRegion) The source volume is empty');
  return
end

nImages = size(volume,3);
aSize = [size(volume,1) size(volume,2)];

h = figure('MenuBar', 'none');
OK = 0;    % Flag to accept chosen crop
while ~OK
    clf
    m = round(sqrt(nImages));
    n = ceil(nImages/m);
    for sliceNum = 1:nImages
        h_slice = subplot(m, n, sliceNum);
        h_image = imagesc(volume(:, :, sliceNum), 'Tag', sprintf(' %d', sliceNum));
        colormap(gray)
        axis off
        axis equal
    end
    brighten(0.6);

    subplot(m, n, 1)
    title('Click on the first slice.')
    sliceNum = 0;
    while sliceNum == 0
        waitforbuttonpress
        tag = get(gco, 'Tag');
        if ~isempty(tag)
            sliceNum = str2num(tag);
        end
    end
    firstSlice = sliceNum;
    
    subplot(m, n, 1)
    title('Click on the last slice.')
    sliceNum = 0;
    while sliceNum == 0
        waitforbuttonpress
        tag = get(gco, 'Tag');
        if ~isempty(tag)
            sliceNum = str2num(tag);
        end
    end
    if sliceNum > firstSlice
        lastSlice = sliceNum;
    else
        lastSlice = firstSlice;
        firstSlice = sliceNum;
    end
    subplot(m, n, firstSlice)
    title('First slice')
    subplot(m, n, lastSlice)
    title('Last slice')

    subplot(m, n, 1)
    title('Click on an image to crop.')
    sliceNum = 0;
    while sliceNum == 0
        waitforbuttonpress
        tag = get(gco, 'Tag');
        if ~isempty(tag)
            sliceNum = str2num(tag);
        end
    end

    clf
    imagesc(volume(:, :, sliceNum));
    colormap(gray)
    brighten(0.6);
    axis off
    axis equal
    title('Crop the image.')

    [x,y] = ginput(2);
    x = sort(round(x));
    y = sort(round(y));

    x = max(1, min(x, size(volume, 2)));
    y = max(1, min(y, size(volume, 1)));
    
    % 
    cropRegion = [y(1), x(1), firstSlice; y(2), x(2), lastSlice];

    % Confirm crop
    tmp = volume([cropRegion(1,1):cropRegion(2,1)], ...
        [cropRegion(1,2):cropRegion(2,2)], ...
        [cropRegion(1,3):cropRegion(2,3)]);

    for sliceNum = 1:size(tmp,3)
        subplot(m, n, sliceNum)

        h_slice(sliceNum) = imagesc(tmp(:, :, sliceNum));
        colormap(gray)
        axis off
        axis equal
    end
    brighten(0.6);

    switch questdlg('Does this look OK?', 'Confirm crop');
        case 'Cancel'
            % Cancel
            cropRegion = [1, 1, 1; size(volume)];
            close(h);
            disp('Crop aborted');
            return
        case 'Yes'
            % Okay
            OK = 1;
        case 'No'
            % No
            disp('Repeating crop');
    end
end  %if ~OK

close(h)

