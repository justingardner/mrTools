function corAnalPlot(view,overlayNum,scan,x,y,z)
%
% corAnal(view,view,overlayNum,scan,x,y,z)
% 
% djh 9/2005

% If corAnal is loaded, then use it. Otherwise, error.
co = viewGet(view,'co');
amp = viewGet(view,'amp');
ph = viewGet(view,'ph');
if isempty(co) | isempty(amp) | isempty(ph)
	mrErrorDlg('corAnal must be loaded.');
end
coval = co.data{scan}(x,y,z);
ampval = amp.data{scan}(x,y,z);
phval = ph.data{scan}(x,y,z);

% Analysis parameters
ncycles = viewGet(view,'ncycles',scan);
detrend = viewGet(view,'detrend',scan);
spatialnorm = viewGet(view,'spatialnorm',scan);
junkframes = viewGet(view,'junkFrames',scan);
nframes = viewGet(view,'nFrames',scan);
framePeriod = viewGet(view,'framePeriod',scan);
highpassPeriod = round(nframes/ncycles);


% get the number of ROIs, we first will 
% look to see if the user clicked on a 
% voxel that is contained in an ROI.
roi = [];
if ~strcmp(viewGet(view,'showROIs'),'hide');
  nROIs = viewGet(view,'numberOfROIs');
else
  nROIs = 0;
end
for roinum = 1:nROIs
  roicoords = getRoiCoordinates(view,roinum,scan);
  % see if this is a matching roi
  if ismember([x y z],roicoords','rows')
    % get the roi
    roi = viewGet(view,'roi',roinum);
  end
end

% now if we have an roi then load its time series
% and get the mean and plot that instead of the
% time series for the voxel
roiPlot = 0;
if ~isempty(roi)
  roi = loadROITSeries(view,roi,viewGet(view,'curScan'),viewGet(view,'curGroup'));
  tseries = mean(roi.tSeries)';
  ptseriesSte = std(100*roi.tSeries/mean(tseries))'/sqrt(roi.n);
  ptseriesSte = ptseriesSte(junkframes+1:junkframes+nframes);
  headerStr = sprintf('Times series from roi %s (n=%i)',roi.name,roi.n);
  roiPlot = 1;
else
  % Load tseries from file. Error if file doesn't exist.
  pathStr = viewGet(view,'tseriesPathStr',scan);
  if ~exist(pathStr,'file')
    mrErrorDlg(['File ',pathStr,' not found']);
  end
  [tseries,hdr] = cbiReadNifti(pathStr,{x,y,z,[]});
  headerStr = sprintf('Times series from voxel [%i %i %i] ',x,y,z);
  tseries = squeeze(tseries);
end

tseries = tseries(junkframes+1:junkframes+nframes);

% Remove dc, convert to percent, detrend, and spatial normalization
ptseries = percentTSeries(tseries,...
    'detrend', detrend,...
    'highpassPeriod', highpassPeriod,...
    'spatialNormalization', spatialnorm,...
    'subtractMean', 'Yes',...
    'temporalNormalization', 'No');

% calculate mean cycle
singleCycle = mean(reshape(ptseries,nframes/ncycles,ncycles)');
singleCycleSte = std(reshape(ptseries,nframes/ncycles,ncycles)')/sqrt(ncycles);

% calculate fft
ft = fft(ptseries);
absft = abs(ft(1:round((length(ft)/2))));
signalAmp = absft(ncycles+1);
noiseFreq = round(2*length(absft)/3):length(absft);
noiseAmp = mean(absft(noiseFreq));
snr = signalAmp/noiseAmp;

% if this is the roi we will need to recompute the corAnal
if roiPlot
  % code for this studiously copied from corAnal.m
  ft = ft(1:1+fix(size(ft, 1)/2), :);
  ampFT = 2*abs(ft)/nframes;
  sumAmp = sqrt(sum(ampFT.^2));
  if sumAmp ~= 0
    coval = ampFT(ncycles+1) ./ sumAmp;
  else
    coval = nan;
  end
  ampval = ampFT(ncycles+1);
  phval = -(pi/2) - angle(ft(ncycles+1,:));
  phval(phval<0) = phval(phval<0)+pi*2;
end

% Create model fit
model = ampval * sin(2*pi*ncycles/nframes * [0:nframes-1]' - phval);

% Change xlabel to be in seconds
% Change ylabel contingent upon detrend to include units

% Select window
selectGraphWin;

% Plot it
subplot(2,3,1:3)
time = linspace(framePeriod,nframes*framePeriod,nframes)';
if roiPlot
  errorbar(time,ptseries,ptseriesSte,'k.-','LineWidth',1);hold on
else
  plot(time,ptseries,'k.-','LineWidth',1);hold on
end
plot(time,model,'r-','LineWidth',1.5);

% Ticks and labels
fontSize = 14;
set(gca,'FontSize',fontSize);
headerStr = sprintf('%s co=%f amp=%f ph=%f (%i deg)\n%s',headerStr,coval,ampval,phval,round(180*phval/pi),viewGet(view,'description',viewGet(view,'curScan')));
title(headerStr);
xtickStep = nframes*framePeriod/ncycles;
xtick = ceil([0:xtickStep:nframes*framePeriod]);
set(gca,'xtick',xtick);
set(gca,'XLim',ceil([0,nframes*framePeriod]));
set(gca,'xgrid','on')
xlabel('Time (sec)');
switch spatialnorm
 case {'Divide by mean'}
  ylabel('fMRI response (% change in image intensity)');
 otherwise
  ylabel('fMRI response(arbitrary units)');
end


% Plot single cylce
subplot(2,3,4)
time = linspace(framePeriod,(nframes/ncycles)*framePeriod,(nframes/ncycles))';
errorbar(time,singleCycle,singleCycleSte,'k.-','LineWidth',1);hold on
plot(time,model(1:length(time)),'r-','LineWidth',1.5);

set(gca,'FontSize',fontSize);
title('Single cycle with STE across cycles');
xlabel('Time (sec)');
switch spatialnorm
 case {'Divide by mean'}
  ylabel('fMRI response (% change in image intensity)');
 otherwise
  ylabel('fMRI response(arbitrary units)');
end
axis tight
xaxis(0,max(time));

% plot fourier amplitude
subplot(2,3,5:6)
plot(0:(length(absft)-1),absft,'k.-');hold on
plot(ncycles,signalAmp,'ro');
plot(noiseFreq-1,absft(noiseFreq),'go');
set(gca,'FontSize',fontSize);
title(sprintf('Stimulus (red): %f Noise (Mean of green): %f CNR: %f',signalAmp,noiseAmp,snr));
ylabel('Magnitude');
xlabel('Fourier component number');

return


