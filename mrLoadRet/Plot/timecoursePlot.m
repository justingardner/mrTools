% timecoursePlot
%
%      usage: timecoursePlot()
%         by: justin gardner
%       date: 03/15/07
%    purpose: 
%
function timecoursePlot(view,overlayNum,scan,x,y,s,roi)


% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end

% select the window to plot into
selectGraphWin;
global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'MenuBar','none');
set(fignum,'NumberTitle','off');
set(fignum,'Name','mrDefaultInterrogator');

junkFrames = viewGet(view, 'junkFrames', s);
nFrames = viewGet(view,'nFrames',scan);

tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));
tSeries = tSeries(junkFrames+1:junkFrames+nFrames);


% get the mean and trend
model = [(1:nFrames);ones(1,nFrames)]';
wgts = model \ tSeries;
fit = model*wgts;

subplot(2,1,1);
plot(tSeries,'k.-');
hold on;
plot(fit, 'r-');
title(sprintf('Voxel: [%i, %i, %i], mean=%0.2f, trend=%0.2f (%% sig change)',x,y,s,wgts(2), 100*wgts(1)/wgts(2)));
xlabel('Volumes');
ylabel('fMRI Signal');
axis tight;
% draw borders between rund
concatInfo = viewGet(view,'concatInfo',scan);
if ~isempty(concatInfo)
  vline(concatInfo.runTransition(2:end,1)-1,'r-');
end

subplot(2,1,2);
fftTSeries = fft(tSeries);
% set mean to zero
fftTSeries(1) = 0;
% plot it 
plot(abs(fftTSeries(1:length(fftTSeries)/2)),'k.-');
title(sprintf('Voxel: [%i %i %i]',x,y,s));
xlabel('FFT components');
ylabel('FFT of fMRI Signal');
axis tight;