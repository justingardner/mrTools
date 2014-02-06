function cropRegion = selectCropRegion(volume)
% function cropRegion = selectCropRegion(volume)

cropRegion = [];
if isempty(volume)
  mrWarnDlg('(selectCropRegion) The source volume is empty');
  return
end

nImages = size(volume,3);
aSize = [size(volume,1) size(volume,2)];

%get the screen size of the default monitor
screenSize = getMonitorPositions;
screenSize = screenSize(1,3:4);
hFigure = figure('MenuBar', 'none');

OK = 0;    % Flag to accept chosen crop
while ~OK
    clf
    %optimize number of columns and rows based on screen dimensions and images dimensions
    [m,n]=getArrayDimensions(nImages,screenSize(2)/aSize(1)*aSize(2)/screenSize(1));
    %optimize figure proportions for subplot
    figureDims = get(hFigure,'position');
    figureDims(3)=figureDims(4)*n/m;
    set(hFigure,'position',figureDims);
    for row = 1:m
      for col = 1:n
        sliceNum = (row-1)*n+col;
        if sliceNum<=nImages
          h(sliceNum) =  subplot('position',getSubplotPosition(col,row,ones(1,n),ones(1,m),0.02,0.02));
          imagesc(volume(:, :, sliceNum), 'Tag', sprintf(' %d', sliceNum));
          colormap(gray)
          axis off
          axis equal
        end
      end
    end
    brighten(0.6);

    set(hFigure,'name','Click on the first slice.')
    sliceNum = 0;
    while sliceNum == 0
        waitforbuttonpress
        tag = get(gco, 'Tag');
        if ~isempty(tag)
            sliceNum = str2num(tag);
        end
    end
    firstSlice = sliceNum;
    hFirst = text(aSize(2)/2,aSize(1)/2,{'First','slice'},'parent',h(firstSlice),'HorizontalAlignment','center','color','g','fontweight','bold');
    
    set(hFigure,'name','Click on the last slice.')
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
        text(aSize(2)/2,aSize(1)/2,{'Last','slice'},'parent',h(lastSlice),'HorizontalAlignment','center','color','g','fontweight','bold');
    else
        lastSlice = firstSlice;
        firstSlice = sliceNum;
        set(hFirst,'string',{'Last','slice'});
        text(aSize(2)/2,aSize(1)/2,{'First','slice'},'parent',h(firstSlice),'HorizontalAlignment','center','color','g','fontweight','bold');
    end

    set(hFigure,'name','Click on an image to crop.')
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
    set(hFigure,'name','Crop the image.')

    [x,y] = ginput(2);
    x = sort(round(x));
    y = sort(round(y));

    x = max(1, min(x, size(volume, 2)));
    y = max(1, min(y, size(volume, 1)));
    
    cropRegion = [y(1), x(1), firstSlice; y(2), x(2), lastSlice];

    % Confirm crop
    tmp = volume(cropRegion(1,1):cropRegion(2,1), ...
        cropRegion(1,2):cropRegion(2,2), ...
        cropRegion(1,3):cropRegion(2,3));
    clf
    %optimize number of columns and rows based on screen dimensions and images dimensions
    [m,n]=getArrayDimensions(size(tmp,3),screenSize(2)/size(tmp,1)*size(tmp,2)/screenSize(1));
    %optimize figure proportions for subplot
    figureDims = get(hFigure,'position');
    figureDims(3)=figureDims(4)*n/m;
    set(hFigure,'position',figureDims);
    for row = 1:m
      for col = 1:n
        sliceNum = (row-1)*n+col;
        if sliceNum<=size(tmp,3)
          h(sliceNum) =  subplot('position',getSubplotPosition(col,row,ones(1,n),ones(1,m),0.02,0.02));
          imagesc(tmp(:, :, sliceNum), 'Tag', sprintf(' %d', sliceNum));
          colormap(gray)
          axis off
          axis equal
        end
      end
    end
    brighten(0.6);

    switch questdlg('Does this look OK?', 'Confirm crop');
        case 'Cancel'
            % Cancel
            cropRegion = [1, 1, 1; size(volume)];
            close(hFigure);
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

close(hFigure)

