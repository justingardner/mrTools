% mrDefaultInterrogator.m
%
%      usage: mrDefaultInterrogator()
%         by: justin gardner
%       date: 03/15/07
%    purpose: 
%
function mrDefaultInterrogator(view,overlayNum,scan,x,y,s)


% check arguments
if ~any(nargin == [1:6])
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

tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));

subplot(2,1,1);
plot(tSeries,'k.-');
title(sprintf('Voxel: [%i %i %i]',x,y,s));
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