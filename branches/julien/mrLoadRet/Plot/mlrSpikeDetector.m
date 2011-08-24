% spikedetector.m
%
%        $Id$	
%      usage: spikeinfo = mlrSpikeDetector(scanNum,groupNum,<'criterion=10'>,<'dispfigs=1'>
%         by: justin gardner
%       date: 01/09/06
%    purpose: looks for spiking artifact in data set
%       e.g.: spikedetector(1,3)
% 
%             criterion sets how many std above mean a fourier
%             component has to be to be considered a
%             spike. (default=10)
%
function v = mlrSpikeDetector(v,scanNum,groupNum,varargin)

% check arguments
if nargin < 3
  help mlrSpikeDetector
  return
end

% get optional arguments
eval(evalargs(varargin));

if exist('criterion')~=1,criterion = 4;,end
if exist('dispfigs')~=1,dispfigs = 1;,end
if exist('recompute')~=1,recompute = 0;,end

% try to load spikeinfo from view
spikeinfo = viewGet(v,'spikeinfo',scanNum,groupNum);

if isempty(spikeinfo) || recompute
  % ask user how to recalculate
  paramsInfo{1} = {'criterion',criterion,'incdec=[-1 1]','minmax=[0 inf]','Criterion for spike detection. Spike detection works by computing the mean and standard deviation of each fourier component of the images across time. If a fourier component on any single volume exceeds criterion standard deviations of the mean, it is considered to be a spike. Default value is 4. i.e. A fourier component has to be 4 standard deviations greater from the mean to be considered to be a spike. This is a low threshold, but can be changed later on'};
  paramsInfo{2} = {'useMedian', 0, 'type=checkbox', 'Use the median and interquartile range to calculate the center and spread of the data.  This is useful if the data are very noisy.'};
  paramsInfo = mrParamsDialogSelectScans(v,groupNum,paramsInfo,scanNum);
  params = mrParamsDialog(paramsInfo,'mlrSpikeDetector params');
  if isempty(params)||~any(params.include),return,end
  % go through and recalculate
  scanNums = find(params.include);
  for s = scanNums
    spikeinfo = calcSpikeInfo(v,s,groupNum,params);
    v = viewSet(v,'spikeinfo',spikeinfo,s,groupNum);
  end
  saveSession;
  % now get back the one we want
  if ~ismember(scanNum,scanNums)
    scanNum = scanNums(1);
  end
  spikeinfo = viewGet(v,'spikeinfo',scanNum,groupNum);
end
if dispfigs
  selectGraphWin;
  spikeinfo.v = v;
  spikeinfo = spikePlot(spikeinfo);
  spikePlotController(spikeinfo);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   calcSpikeInfo   %%
%%%%%%%%%%%%%%%%%%%%%%%
function spikeinfo = calcSpikeInfo(v,scanNum,groupNum,params)

% load TSeries
spikeinfo.scanNum = scanNum;
spikeinfo.groupNum = groupNum;
spikeinfo.filename = viewGet(v,'tSeriesFile',scanNum,groupNum);
if isempty(spikeinfo.filename)
  disp(sprintf('(mlrSpikeDetector) Could not find scan %i in group %i',scanNum,groupNum));
  spikeinfo = [];
  return
end
  
% load file
v = viewSet(v,'curGroup',groupNum);
disppercent(-inf,sprintf('(mlrSpikeDetector) Loading time series for scan %i, group %i',scanNum,groupNum));
data = loadTSeries(v,scanNum);
% Dump junk frames
junkFrames = viewGet(v, 'junkframes', scanNum);
nFrames = viewGet(v, 'nFrames', scanNum);
data = data(:,:,:,junkFrames+1:junkFrames+nFrames);

spikeinfo.dim = size(data);
disppercent(inf);

% compute timecourse means for later 
for slicenum = 1:spikeinfo.dim(3)
  % calculate mean timecourse for slice
  sliceMeans(slicenum,:) = nanmean(nanmean(data(:,:,slicenum,:),1),2);
  % normalize to % signal change
  sliceMeans(slicenum,:) = 100*sliceMeans(slicenum,:)/mean(sliceMeans(slicenum,:));
  % compute standard deviation
  sliceStds(slicenum) = std(sliceMeans(slicenum,:));
end

% compute fourier transform of data
% calculating fourier transform of data
disppercent(-inf,'(mlrSpikeDetector) Calculating FFT');
% skip some frames in the beginning to account
% for saturation
if junkFrames < 5
    startframe = min(5,spikeinfo.dim(4));
else
    startframe = 1;
end
data = data(:,:,:,startframe:spikeinfo.dim(4));
for i = 1:size(data,4)
  disppercent(i/size(data,4));
  for j = 1:size(data,3)
    %first need to remove NaNs from data
    %let's replace them by the mean of each image
    thisData =  data(:,:,j,i);
    thisData(isnan(thisData)) = nanmean(thisData(:));
    %compute the spatial fourier transform
    data(:,:,j,i) = abs(fftshift(fft2(thisData)));
  end
end
disppercent(inf);

% get mean and std
if params.useMedian
    disppercent(-inf,'(mlrSpikeDetector) Calculating median and iqr');
    for slicenum = 1:spikeinfo.dim(3)
        disppercent(slicenum/spikeinfo.dim(3));
        meandata(:,:,slicenum) = squeeze(median(data(:,:,slicenum,:),4));
        stddata(:,:,slicenum) = squeeze(iqr(data(:,:,slicenum,:),4));
    end
else
    disppercent(-inf,'(mlrSpikeDetector) Calculating mean and std');
    for slicenum = 1:spikeinfo.dim(3)
        meandata(:,:,slicenum) = squeeze(mean(data(:,:,slicenum,:),4));
        stddata(:,:,slicenum) = squeeze(std(data(:,:,slicenum,:),0,4));
    end
end
disppercent(inf);


% now subtract off mean and see
% if there are any points above std criterion
slice = [];time = [];numspikes = [];spikelocs = {};meanZvalue=[];
disppercent(-inf,'(mlrSpikeDetector) Looking for spikes');
for i = 1:size(data,4)
  disppercent(i/spikeinfo.dim(4));
  data(:,:,:,i) = squeeze(data(:,:,:,i))-meandata;
  % see if any voxels are larger then expected
  for slicenum = 1:spikeinfo.dim(3)
    [spikex spikey] = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
    if ~isempty(spikex)
      slice(end+1) = slicenum;
      time(end+1) = startframe + i -1;
      numspikes(end+1) = length(spikex);
      spikelocs{end+1}.x = spikex;
      spikelocs{end}.y = spikey;
      spikelocs{end}.linear = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
      for iSpike = 1:length(spikex)
        spikelocs{end}.zValue(iSpike) =  data(spikex(iSpike),spikey(iSpike),slicenum,i)/squeeze(stddata(spikex(iSpike),spikey(iSpike),slicenum));
      end
      meanZvalue(end+1) = mean(spikelocs{end}.zValue);
    end
  end
end
disppercent(inf);
if length(slice)
  disp(sprintf('======================================================'));
  disp(sprintf('(mlrSpikeDetector) Found %i spikes at z>%.2f in scan %i, group %i',length(slice),params.criterion,scanNum,groupNum));
  disp(sprintf('======================================================'));
else
  disp(sprintf('(mlrSpikeDetector) No spikes in scan %i, group %i',scanNum,groupNum));
end

% now look for any time point that exceeds criteria
% (in std units) -- this way looks averaged over slice
%[slice time] = find(abs(sliceMeans(:,:)-100) > (sliceStds')*ones(1,spikeinfo.dim(4))*criterion);

% pass them back in spikeinfo
spikeinfo.n = length(slice);
spikeinfo.slice = slice;
spikeinfo.time = time;
spikeinfo.numspikes = numspikes;
spikeinfo.spikelocs = spikelocs;
spikeinfo.criterion = params.criterion;
spikeinfo.sliceMeans = single(sliceMeans);
spikeinfo.meanZvalue = meanZvalue;
spikeinfo.maxZvalue = [zeros(spikeinfo.dim(3),startframe-1) permute(max(max(data./repmat(stddata,[1 1 1 size(data,4)]),[],1),[],2),[3 4 1 2])];

%%%%%%%%%%%%%%%%%%%
%%   spikePlot   %%
%%%%%%%%%%%%%%%%%%%
function hFigure = spikePlot(spikeinfo)

if isempty(spikeinfo)
  return
end

if fieldIsNotDefined(spikeinfo,'currentCriterion')
  spikeinfo.currentCriterion=spikeinfo.criterion;
else
  spikeinfo.currentCriterion=max(spikeinfo.criterion,spikeinfo.currentCriterion);
end
  
% plot the slice means
hFigure = selectGraphWin(0,'replace');
colorOrder = jet(spikeinfo.dim(3));
set(hFigure,'DefaultAxesColorOrder',colorOrder,'name',...
  sprintf('Spike Detector - %s:%i %s (%s)',...
           viewGet(spikeinfo.v,'groupName',spikeinfo.groupNum),spikeinfo.scanNum,...
           viewGet(spikeinfo.v,'description',spikeinfo.scanNum,spikeinfo.groupNum),...
           spikeinfo.filename));

% if spikeinfo.n == 0
%   spikeinfo.nrows = 4;
%   spikeinfo.ncols = 3;
%   plot1 = 1:3;
%   plot2 = [4:12];
%   
% else
  spikeinfo.nrows = 4;
  spikeinfo.ncols = 5;
  plot1 = 1:3;
  plot2 = [6:8,11:13,16:18];
  spikeinfo.imagePlot = subplot(spikeinfo.nrows,spikeinfo.ncols,[4 5 9 10]);
  spikeinfo.fftPlot = subplot(spikeinfo.nrows,spikeinfo.ncols,[14 15 19 20]);
% end

spikeinfo.hTseries = subplot(spikeinfo.nrows,spikeinfo.ncols,plot1);
hold on
for slicenum = 1:spikeinfo.dim(3)
  plot(spikeinfo.sliceMeans(slicenum,:),'.-','color',colorOrder(slicenum,:));
end
set(spikeinfo.hTseries,'xLim',[0 spikeinfo.dim(4)]);
spikeinfo.tSeriesYlim = get(spikeinfo.hTseries,'Ylim');
spikeinfo.hCursorTseries = plot(spikeinfo.hTseries,[1 1],spikeinfo.tSeriesYlim,'k','visible','off');
%title(sprintf('%s:%i %s (%s)\nSpikes found=%i (criterion=%0.1f)',viewGet(spikeinfo.v,'groupName',spikeinfo.groupNum),spikeinfo.scanNum,viewGet(spikeinfo.v,'description',spikeinfo.scanNum,spikeinfo.groupNum),spikeinfo.filename,spikeinfo.n,spikeinfo.criterion),'interpreter','none');

% print out slice labels
xmax = min(get(gca,'XTick'));
ymax = max(get(gca,'yTick'));
for slicenum = 0:spikeinfo.dim(3)
  if slicenum==0
    htext = text(xmax,ymax,'Slices:');
  else
    htext = text(xmax,ymax,sprintf('%i',slicenum));
    set(htext,'Color',colorOrder(slicenum,:));
  end
  position = get(htext,'Position');
  textextent = get(htext,'Extent');
  position(2) = position(2)-textextent(4);
  position(1) = position(1)+10;
  set(htext,'Position',position);
  xmax = xmax+textextent(3);
end
ylabel(spikeinfo.hTseries,{'Mean over each slice','(% signal change)'});

% plot Z value Matrix and Spikes
spikeinfo.hSpikeMatrix = subplot(spikeinfo.nrows,spikeinfo.ncols,plot2);
title(spikeinfo.hSpikeMatrix,'Spike matrix: Use mouse or arrow keys to move the cursor and display corresponding image, FFT and spikes');
hold on
if isfield(spikeinfo,'maxZvalue')
  imagesc(spikeinfo.maxZvalue);
else
  imagesc(zeros(spikeinfo.dim(3),spikeinfo.dim(4)));
end
colormap('gray');
caxis([0 8]);
set(spikeinfo.hSpikeMatrix,'xLim',[0.5 spikeinfo.dim(4)+.5],'yLim',[0.5 spikeinfo.dim(3)+.5]);
xlabel(sprintf('Frame number'));
ylabel('Slice number');
h = colorbar('location','southoutside');
set(get(h,'XLabel'),'string','Max FFT Z value in Slice/Frame')

if ~spikeinfo.n
    spikeinfo.currentSpike =0;
elseif fieldIsNotDefined(spikeinfo,'currentSpike')
    spikeinfo.currentSpike =1;
end
if fieldIsNotDefined(spikeinfo,'currentFrame') || fieldIsNotDefined(spikeinfo,'currentSlice')
  if spikeinfo.currentSpike
    spikeinfo.currentFrame =spikeinfo.time(spikeinfo.currentSpike);
    spikeinfo.currentSlice =spikeinfo.slice(spikeinfo.currentSpike);
  else
    spikeinfo.currentFrame =round(spikeinfo.dim(4)/2);
    spikeinfo.currentSlice =round(spikeinfo.dim(3)/2);
  end    
end

guidata(hFigure,spikeinfo);

% display all spikes
setCriterion(hFigure,spikeinfo.currentCriterion);

set(hFigure,'WindowButtonDownFcn',{@MouseDownCallback});
set(hFigure,'WindowButtonMotionFcn',{@MouseMotionCallback});
set(hFigure,'KeyPressFcn',{@KeypressCallback});
set(hFigure,'interruptible','off')
set(hFigure,'BusyAction','cancel');

% set(hFigure,'WindowButtonUpFcn',);
%drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MouseDownCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseDownCallback(hFigure,eventData)

spikeinfo = guidata(hFigure);

spikeMatrixCoords = round(get(spikeinfo.hSpikeMatrix,'CurrentPoint'));
if all(spikeMatrixCoords(1,[1 2])>0 & spikeMatrixCoords(1,[1 2])<=spikeinfo.dim([4 3]))
  spikeinfo.currentSlice = spikeMatrixCoords(1,2);
  spikeinfo.currentFrame = spikeMatrixCoords(1,1);
  guidata(hFigure,spikeinfo);
  spikePlotImage(hFigure,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MouseMotionCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseMotionCallback(hFigure,eventData)

spikeinfo = guidata(hFigure);

spikeMatrixCoords = round(get(spikeinfo.hSpikeMatrix,'CurrentPoint'));
tSeriesCoords = round(get(spikeinfo.hTseries,'CurrentPoint'));
if ishandle(spikeinfo.hCursorTseries)
  delete(spikeinfo.hCursorTseries);
end
if all(spikeMatrixCoords(1,[1 2])>0 & spikeMatrixCoords(1,[1 2])<=spikeinfo.dim([4 3])) ||...
  all(tSeriesCoords(1,[1 2])>[0 spikeinfo.tSeriesYlim(1)] & tSeriesCoords(1,[1 2])<=[spikeinfo.dim(4) spikeinfo.tSeriesYlim(2)])
  spikeinfo.hCursorTseries = plot(spikeinfo.hTseries,repmat(tSeriesCoords(1,1),1,2),spikeinfo.tSeriesYlim,'k');
end

guidata(hFigure,spikeinfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   KeypressCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeypressCallback(hFigure,eventData)

spikeinfo = guidata(hFigure);
switch(double(get(hFigure,'CurrentCharacter')))
  case 28 %left arrow
    spikeinfo.currentFrame = max(spikeinfo.currentFrame-1,1);
  case 29 %right arrow
    spikeinfo.currentFrame = min(spikeinfo.currentFrame+1,spikeinfo.dim(4));
  case 30 %up arrow
    spikeinfo.currentSlice = min(spikeinfo.currentSlice+1,spikeinfo.dim(3));
  case 31 %down arrow
    spikeinfo.currentSlice = max(spikeinfo.currentSlice-1,1);
  otherwise
    return;
end
guidata(hFigure,spikeinfo);
   
spikePlotImage(hFigure,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotController   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotController(hFigure)

spikeinfo = guidata(hFigure);

% check all scans that have spikeinfo
spikeInfoScans = {};
for i = 1:viewGet(spikeinfo.v,'nScans',spikeinfo.groupNum);
  if i==spikeinfo.scanNum %for the current scan, get the current spikeinfo, not the one from the view
    thisSpikeinfo=spikeinfo;
  else
    thisSpikeinfo = viewGet(spikeinfo.v,'spikeinfo',i,spikeinfo.groupNum);
  end
  if ~isempty(thisSpikeinfo);
    if fieldIsNotDefined(thisSpikeinfo,'currentCriterion')
      criterion=thisSpikeinfo.criterion;
      spikeNumber = thisSpikeinfo.n;
    else
      criterion=thisSpikeinfo.currentCriterion;
      spikeNumber = thisSpikeinfo.currentSpikeNumber;
    end
    spikeInfoScans{end+1} = sprintf('%i: %s (%i spikes at z>%.2f)',...
      i,viewGet(spikeinfo.v,'description',i,thisSpikeinfo.groupNum),...
      spikeNumber,criterion);
    if i==spikeinfo.scanNum
      thisSpikeInfoScans = spikeInfoScans{end};
    end
  end
end
spikeInfoScans = putOnTopOfList(thisSpikeInfoScans,spikeInfoScans);

% now put up a control dialog
paramsInfo{1}  = {'scanNum',spikeInfoScans,'Scan number to view'};
paramsInfo{end+1} = {'recompute',[],'type=pushbutton','buttonString=Recompute spike detection','callback',@spikePlotRecomputeCallback,'callbackArg',hFigure,'Recomputer spike detection'};
if spikeinfo.n > 0
  paramsInfo{end+1} = {'spikeNum',1,sprintf('minmax=[1 %i]',spikeinfo.n),'incdec=[-1 1]','round=1','Which spike to display'};
  paramsInfo{end+1} = {'criterion',max(spikeinfo.criterion,spikeinfo.currentCriterion),sprintf('minmax=[%i inf]',spikeinfo.criterion),'incdec=[-.5 .5]','To change the value of the criterion used (must be > criterion used for computing the detection)'};
else
  paramsInfo{end+1} = {'noSpikes','No spikes found','editable=0','type=string','No spikes found in scan'};
end  
mrParamsDialog(paramsInfo,'mlrSpikeDetector',[],@spikePlotCallback,hFigure,{@spikePlotOKCallback,hFigure});


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotCallback(params,hFigure)

spikeinfo = guidata(hFigure);

% get the scanNum
scanNum = str2num(strtok(params.scanNum,':'));

if scanNum ~= spikeinfo.scanNum
  %keep view but remove it from the current spikeinfo
  v = spikeinfo.v;
  spikeinfo = rmfield(spikeinfo,'v');
  %save the current spikeinfo
  v = viewSet(v,'spikeinfo',spikeinfo,spikeinfo.scanNum,spikeinfo.groupNum);
  % load up thhe new spike info
  spikeinfo = viewGet(v,'spikeinfo',scanNum,spikeinfo.groupNum);
  spikeinfo.v=v;
  hFigure=spikePlot(spikeinfo);
  spikePlotController(hFigure);
else
  if isfield(params,'criterion') && params.criterion~=spikeinfo.currentCriterion
    setCriterion(hFigure,params.criterion);
    %change scan description
    %relaunch spikePlotcontroller (cannot use mrParamsSet to change the string of a popupmenu...)
    spikePlotController(hFigure);
  elseif isfield(params,'spikeNum') && params.spikeNum~=spikeinfo.currentSpike
    spikePlotImage(hFigure,params.spikeNum);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotOKCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotOKCallback(hFigure)

spikeinfo = guidata(hFigure);

%keep view but remove it from the current spikeinfo
v = spikeinfo.v;
spikeinfo = rmfield(spikeinfo,'v');
%save the current spikeinfo
v = viewSet(v,'spikeinfo',spikeinfo,spikeinfo.scanNum,spikeinfo.groupNum);
h = selectGraphWin(0,'replace');
closeGraphWin(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = spikePlotRecomputeCallback(hFigure)

val = [];
spikeinfo = guidata(hFigure);
%keep view but remove it from the current spikeinfo
v = spikeinfo.v;
spikeinfo = rmfield(spikeinfo,'v');
%save the current spikeinfo
v = viewSet(v,'spikeinfo',spikeinfo,spikeinfo.scanNum,spikeinfo.groupNum);
close;
h = selectGraphWin(0,'replace');
closeGraphWin(h);

mlrSpikeDetector(v,spikeinfo.scanNum,spikeinfo.groupNum,'recompute=1');
%%%%%%%%%%%%%%%%%%%%%%%%
%%   setCriterion   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function setCriterion(hFigure,thisCriterion)

spikeinfo = guidata(hFigure);

spikeinfo.currentCriterion = thisCriterion;

if isfield(spikeinfo,'hBox') && all(ishandle(spikeinfo.hBox(:)))
  delete(spikeinfo.hBox);
  spikeinfo = rmfield(spikeinfo,'hBox');
end

%draw a contour around each slice/frame that contains spike over the current criterion
boxXcoords = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
boxYcoords = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
cSpike=0;
if spikeinfo.n && isfield(spikeinfo.spikelocs{1},'zValue')
  for iSpike = 1:spikeinfo.n
    if any(spikeinfo.spikelocs{iSpike}.zValue>thisCriterion)
      cSpike = cSpike+1;
      spikeinfo.hBox(cSpike,:) = plot(spikeinfo.hSpikeMatrix,spikeinfo.time(iSpike)+boxXcoords,spikeinfo.slice(iSpike)+boxYcoords,'g');
    end
  end
else  %for an old spikeinfo structure, we don't have the zValue info
  mrWarnDlg('(mlrSpikeDetector:setCriterion) There is no zValue information in spikeinfo, please recompute spike detection for this scan.');
  for cSpike = 1:spikeinfo.n
      spikeinfo.hBox(cSpike,:) = plot(spikeinfo.hSpikeMatrix,spikeinfo.time(cSpike)+boxXcoords,spikeinfo.slice(cSpike)+boxYcoords,'g');
  end
end

spikeinfo.currentSpikeNumber = cSpike;
title(spikeinfo.hTseries,...
      {sprintf('%i spikes found at z>%0.1f)',cSpike,thisCriterion),...
      '(Note that spike detection is done on each FFT component not on the mean of the slice)'},...
      'interpreter','none');

guidata(hFigure,spikeinfo); 

% display images with possible artifacts for the current slice/frame
spikePlotImage(hFigure,spikeinfo.currentSpike);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotImage(hFigure,spikeNum)

spikeinfo = guidata(hFigure);

%if there is no spikenum passed (0), see if there are spikes for this slice and frame
if ~spikeNum
  spikeNum = find(spikeinfo.slice==spikeinfo.currentSlice & spikeinfo.time==spikeinfo.currentFrame);
  if isempty(spikeNum)
    spikeNum=0;
  end
else
%otherwise, set the the current slice and frame on this spike
  spikeinfo.currentFrame =spikeinfo.time(spikeNum);
  spikeinfo.currentSlice =spikeinfo.slice(spikeNum);
end
spikeinfo.currentSpike = spikeNum;
thisFrame = spikeinfo.currentFrame;
thisSlice = spikeinfo.currentSlice;

%draw the cursor
if isfield(spikeinfo,'hCursor') && all(ishandle(spikeinfo.hCursor))
  delete(spikeinfo.hCursor);
end
boxXcoords = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
boxYcoords = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
spikeinfo.hCursor = plot(spikeinfo.hSpikeMatrix,thisFrame+boxXcoords,thisSlice+boxYcoords,'m','linewidth',3);
drawnow

% display image
% read data
spikeinfo.v = viewSet(spikeinfo.v,'curGroup',spikeinfo.groupNum);
data = loadTSeries(spikeinfo.v,spikeinfo.scanNum,thisSlice,thisFrame);
imageg(squeeze(data),0.5,spikeinfo.imagePlot);
% and title
title(spikeinfo.imagePlot,{sprintf('Slice %i, Frame %i',thisSlice,thisFrame),'(Normalized intensity values)'});
fftimage = abs(fftshift(fft2(squeeze(data))));
% set gamma
imageg(fftimage,0.5,spikeinfo.fftPlot);
title(spikeinfo.fftPlot,{'FFT',sprintf('Components highlighted in green: z>%.2f',spikeinfo.currentCriterion)});

if spikeNum
  hold(spikeinfo.fftPlot,'on')
  % draw boxes around noise values
  if isfield(spikeinfo.spikelocs{1},'zValue')
    for iLoc = 1:length(spikeinfo.spikelocs{spikeNum}.x)
      if spikeinfo.spikelocs{spikeNum}.zValue(iLoc)>spikeinfo.currentCriterion
        plot(spikeinfo.fftPlot,spikeinfo.spikelocs{spikeNum}.y(iLoc)+boxXcoords,spikeinfo.spikelocs{spikeNum}.x(iLoc)+boxYcoords,'g','linewidth',1); %why should x and y be inverted (and does that matter ?)
      end
    end
  else  %for an old spikeinfo structure, we don't have the zValue info
    mrWarnDlg('(mlrSpikeDetector:setCriterion) There is no zValue information in spikeinfo, please recompute spike detection for this scan.');
    for iLoc = 1:length(spikeinfo.spikelocs{spikeNum}.x)
      plot(spikeinfo.fftPlot,spikeinfo.spikelocs{spikeNum}.y(iLoc)+boxXcoords,spikeinfo.spikelocs{spikeNum}.x(iLoc)+boxYcoords,'g','linewidth',1); %why should x and y be inverted (and does that matter ?)
    end
  end
  hold(spikeinfo.fftPlot,'off')
end

% % and fix up axis
% axis off
% axis equal

guidata(hFigure,spikeinfo);

% imageg.m
%
%      usage: imageg(data array,gamma)
%      usage: imageg(data structure,slicenum,volumenum,gamma)
%         by: justin gardner
%       date: 05/09/03
%    purpose: gamma correct image display for epi images
%
function d = imageg(d,gamma,handle)

if size(d,2) < size(d,1)
  d = d';
end

% scale image values to between 0 and 1
imagemax = max(max(d));
imagemin = min(min(d));
d = (d-imagemin)./(imagemax-imagemin);

% apply monitor gamma
d = d.^gamma;

% rescale to 0-255 uint8
d = floor(255*d);

image(d,'parent',handle);
colormap(gray(256));
axis(handle,'off');
axis(handle,'equal');
