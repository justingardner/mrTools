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

if exist('criterion')~=1,criterion = 10;,end
if exist('dispfigs')~=1,dispfigs = 1;,end

% try to load spikeinfo from view
spikeinfo = viewGet(v,'spikeinfo',scanNum,groupNum);

if isempty(spikeinfo)
  % ask user how to recalculate
  paramsInfo{1} = {'criterion',criterion,'incdec=[-1 1]','minmax=[0 inf]','Criterion for spike detection. Spike detection works by computing the mean and standard deviation of each fourier component of the images across time. If a fourier component on any single volume exceeds criterion standard deviations of the mean, it is considered to be a spike. Default value is 10. i.e. A fourier component has to be 10 standard deviations greater from the mean to be considered to be a spike.'};
  paramsInfo{2} = {'useMedian', 0, 'type=checkbox', 'Use the median and interquartile range to calculate the center and spread of the data.  This is useful if the data are very noisy.'};
  paramsInfo = mrParamsDialogSelectScans(v,groupNum,paramsInfo,scanNum);
  params = mrParamsDialog(paramsInfo,'mlrSpikeDetector params');
  if isempty(params),return,end
  % go through and recalculate
  scanNums = find(params.include);
  for s = scanNums
    spikeinfo = calcSpikeInfo(v,s,groupNum,params);
    v = viewSet(v,'spikeinfo',spikeinfo,s,groupNum);
  end
  saveSession;
  % now get back the one we want
  spikeinfo = viewGet(v,'spikeinfo',scanNum,groupNum);
end
if dispfigs
  spikePlotController(v,spikeinfo);
  spikePlot(v,spikeinfo);
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
spikeinfo.data = loadTSeries(v,scanNum);
% Dump junk frames
junkFrames = viewGet(v, 'junkframes', scanNum);
nFrames = viewGet(v, 'nFrames', scanNum);
spikeinfo.data = spikeinfo.data(:,:,:,junkFrames+1:junkFrames+nFrames);

spikeinfo.dim = size(spikeinfo.data);
disppercent(inf);

% skip some frames in the beginning to account
% for saturation
if junkFrames < 5
    startframe = min(5,spikeinfo.dim(4));
else
    startframe = 1;
end

% compute fourier transform of data
% calculating fourier transform of data
disppercent(-inf,'(mlrSpikeDetector) Calculating FFT');
data = zeros(spikeinfo.dim);
for i = startframe:spikeinfo.dim(4)
  disppercent(i/spikeinfo.dim(4));
  for j = 1:spikeinfo.dim(3)
    data(:,:,j,i) = abs(fftshift(fft2(squeeze(spikeinfo.data(:,:,j,i)))));
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
slice = [];time = [];numspikes = [];spikelocs = {};
disppercent(-inf,'(mlrSpikeDetector) Looking for spikes');
for i = startframe:spikeinfo.dim(4)
  disppercent(i/spikeinfo.dim(4));
  data(:,:,:,i) = squeeze(data(:,:,:,i))-meandata;
  % see if any voxels are larger then expected
  for slicenum = 1:spikeinfo.dim(3)
    [spikex spikey] = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
    if ~isempty(spikex)
      slice(end+1) = slicenum;
      time(end+1) = i;
      numspikes(end+1) = length(spikex);
      spikelocs{end+1}.x = spikex;
      spikelocs{end}.y = spikey;
      spikelocs{end}.linear = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
    end
  end
end
disppercent(inf);
if length(slice)
  disp(sprintf('======================================================'));
  disp(sprintf('(mlrSpikeDetector) Found %i spikes in scan %i, group %i',length(slice),scanNum,groupNum));
  disp(sprintf('======================================================'));
else
  disp(sprintf('(mlrSpikeDetector) No spikes in scan %i, group %i',scanNum,groupNum));
end

% compute timecourse means for each 
for slicenum = 1:spikeinfo.dim(3)
  % calculate mean timecourse for slice
  sliceMeans(slicenum,:) = mean(mean(spikeinfo.data(:,:,slicenum,:),1),2);
  % normalize to % signal change
  sliceMeans(slicenum,:) = 100*sliceMeans(slicenum,:)/mean(sliceMeans(slicenum,:));
  % compute standard deviation
  sliceStds(slicenum) = std(sliceMeans(slicenum,:));
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
spikeinfo = rmfield(spikeinfo,'data');

%%%%%%%%%%%%%%%%%%%
%%   spikePlot   %%
%%%%%%%%%%%%%%%%%%%
function spikePlot(v,spikeinfo);

if isempty(spikeinfo)
  return
end

% plot the slice means
selectGraphWin;
coloroffset = round(rand*10);
if spikeinfo.n == 0
  nrows = 1;ncols = 1;
  plot1 = 1;
else
  nrows = 2;ncols = 2;
  plot1 = 1:2;
end
for slicenum = 1:spikeinfo.dim(3)
  subplot(nrows,ncols,plot1)
  plot(spikeinfo.sliceMeans(slicenum,:),getcolor(slicenum+coloroffset,'.-'));
  hold on
end
title(sprintf('%s:%i %s (%s)\nSpikes found=%i (criterion=%0.1f)',viewGet(v,'groupName',spikeinfo.groupNum),spikeinfo.scanNum,viewGet(v,'description',spikeinfo.scanNum,spikeinfo.groupNum),spikeinfo.filename,spikeinfo.n,spikeinfo.criterion),'interpreter','none');

% plot lines where artificats may be occurring
ytop = 90;
ybot = 95;
if spikeinfo.n < 50
  for slicenum = 1:spikeinfo.dim(3)
    % plot lines where there is an artifact
    thisslice = find(spikeinfo.slice == slicenum);
    if ~isempty(thisslice)
      for i = 1:length(thisslice)
	plot([spikeinfo.time(thisslice(i)) spikeinfo.time(thisslice(i))],[ytop ybot],'Color',getcolor(slicenum+coloroffset));
      end
    end
    drawnow;
  end
else
  disp(sprintf('(mlrSpikeDetector) Too many spikes (%i) to display lines on graphs',spikeinfo.n));
end

% print out slice labels
xmax = min(get(gca,'XTick'));
ymax = max(get(gca,'yTick'));
for slicenum = 1:spikeinfo.dim(3)
  htext = text(xmax,ymax,sprintf('%i',slicenum));
  set(htext,'Color',getcolor(slicenum+coloroffset));
  position = get(htext,'Position');
  textextent = get(htext,'Extent');
  position(2) = position(2)-textextent(4);
  position(1) = position(1)+10;
  set(htext,'Position',position);
  xmax = xmax+textextent(3);
end

xlabel(sprintf('Volume number\nNote that spike detection is done on each FFT component not on the mean of the slice.'));
ylabel('Mean over each slice (% signal change)');
zoom on

fignum = 0;

% display images with possible artifacts
spikePlotImage(v,spikeinfo,1);
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotController   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotController(v,spikeinfo)

% check all scans that have spikeinfo
spikeInfoScans = {};
for i = 1:viewGet(v,'nScans',spikeinfo.groupNum);
  thisSpikeinfo = viewGet(v,'spikeinfo',i,spikeinfo.groupNum);
  if ~isempty(thisSpikeinfo);
    spikeInfoScans{end+1} = sprintf('%i: %s (%i spikes)',i,viewGet(v,'description',i,spikeinfo.groupNum),thisSpikeinfo.n);
  end
end
spikeInfoScans = putOnTopOfList(sprintf('%i: %s (%i spikes)',spikeinfo.scanNum,viewGet(v,'description',spikeinfo.scanNum,spikeinfo.groupNum),spikeinfo.n),spikeInfoScans);

% now put up a control dialog
spikeinfo.v = v;
paramsInfo{1}  = {'scanNum',spikeInfoScans,'Scan number to view'};
paramsInfo{end+1} = {'recompute',[],'type=pushbutton','buttonString=Recompute spike detection','callback',@spikePlotRecomputeCallback,'callbackArg',spikeinfo,'Recomputer spike detection'};
if spikeinfo.n > 0
  paramsInfo{end+1} = {'spikeNum',1,sprintf('minmax=[1 %i]',spikeinfo.n),'incdec=[-1 1]','round=1','Which spike to display'};
else
  paramsInfo{end+1} = {'noSpikes','No spikes found','editable=0','type=string','No spikes found in scan'};
end  
mrParamsDialog(paramsInfo,'mlrSpikeDetector',[],@spikePlotCallback,spikeinfo,@spikePlotOKCallback);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotCallback(params,spikeinfo)

% get the scanNum
scanNum = str2num(strtok(params.scanNum,':'));

if scanNum ~= spikeinfo.scanNum
  % load up that spike info
  v = spikeinfo.v;
  spikeinfo = viewGet(v,'spikeinfo',scanNum,spikeinfo.groupNum);
  spikePlotController(v,spikeinfo);
  spikePlot(v,spikeinfo);
elseif isfield(params,'spikeNum')
  selectGraphWin(1);hold on
  subplot(2,2,3);cla;subplot(2,2,4);cla;
  spikePlotImage(spikeinfo.v,spikeinfo,params.spikeNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotOKCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotOKCallback(spikeinfo)

selectGraphWin;
closeGraphWin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = spikePlotRecomputeCallback(spikeinfo)

val = [];
close;
selectGraphWin;
closeGraphWin;
spikeinfo.v = viewSet(spikeinfo.v,'spikeinfo',[],spikeinfo.scanNum,spikeinfo.groupNum);
mlrSpikeDetector(spikeinfo.v,spikeinfo.scanNum,spikeinfo.groupNum);
%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotImage(v,spikeinfo,imageNum)

if imageNum > spikeinfo.n
  return
end
% get the rows and cols
thisslice = spikeinfo.slice(imageNum);
thistime = spikeinfo.time(imageNum);
% make the figure
nrows = 2;ncols = 2;
% display image
subplot(nrows,ncols,3);
% read data
v = viewSet(v,'curGroup',spikeinfo.groupNum);
data = loadTSeries(v,spikeinfo.scanNum,thisslice,thistime);
imageg(squeeze(data),0.5);
% and title
title(sprintf('slice=%i TR=%i',thisslice,thistime));
subplot(nrows,ncols,4);
fftimage = abs(fftshift(fft2(squeeze(data))));
% set gamma
fftimage = imageg(fftimage,0.5);
% set noise values to higher values
fftimage(spikeinfo.spikelocs{imageNum}.linear) = 0.5*fftimage(spikeinfo.spikelocs{imageNum}.linear)+1.5*256;
image(fftimage);
% and title
title(sprintf('FFT\nslice=%i TR=%i',thisslice,thistime));
% and fix up axis
axis off
axis equal
zoom on
colormap([gray(256) ;gray(256)*[[0 0 0];[0 1 0];[0 0 1]]]);


% imageg.m
%
%      usage: imageg(data array,gamma)
%      usage: imageg(data structure,slicenum,volumenum,gamma)
%         by: justin gardner
%       date: 05/09/03
%    purpose: gamma correct image display for epi images
%
function d = imageg(d,arg1,arg2,arg3)

% make sure we have data passed in
if (nargin <1)
  help imageg;
  return
end

% default gamma value
gamma = 0.8;
    
% if a structure is passed in then arg1 and arg2 are
% the slicenum and volumenum
if (isstruct(d))
  if (nargin == 2) 
    d = d.data(:,:,arg1,1);
  elseif (nargin >= 3)
    d = d.data(:,:,arg1,arg2);
  end
  % get gamma setting
  if (nargin == 4)
    gamma = arg3;
  end
  % too many arguments
  if (nargin > 4)
    help imageg;
    return
  end
% not a structure then the only other argument is gamma
else
  if (nargin == 2)
    gamma = arg1;
  elseif (nargin > 2)
    help imageg;
    return
  end
end

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

% only display with no output arguments
if (nargout == 0)
  image(d);
  colormap(gray(256));
  axis off
  axis equal
end