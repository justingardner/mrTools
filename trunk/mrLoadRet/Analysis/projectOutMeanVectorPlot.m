% projectOutMeanVectorPlot.m
%
%        $Id$
%      usage: projectOutMeanVectorPlot()
%         by: justin gardner
%       date: 05/02/08
%    purpose: 
%
function projectOutMeanVectorPlot(view,overlayNum,scan,x,y,s,roi)

% check arguments
if ~any(nargin == [1:7])
  help projectOutMeanVectorPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(projectOutMeanVectorPlot) projectionAnal not for scan %i',scan));
  return
end
frameNums = [];
% get overlay num and name
overlayNum = viewGet(view,'curOverlay');
overlayName = viewGet(view,'overlayName');

% figure out which d it corresponds to.
if strcmp(overlayName,'mag')
  dnum = 1;
elseif ~isempty(regexp(overlayName,'[0-9]+'))
  numloc = regexp(overlayName,'[0-9]+');
  dnum = str2num(overlayName(numloc:end));
  concatInfo = viewGet(view,'concatInfo',scan);
  frameNums = concatInfo.runTransition(dnum,:);
elseif strcmp(overlayName,'magmean');
  disp(sprintf('(projectOutMeanVectorPlot) No plot for mean'));
  return
else
  disp(sprintf('(projectOutMeanVectorPlot) Unknown overlay name %s',overlayName));
  return
end

% check to make sure overlay exists
if isempty(analysis.overlays(overlayNum).data{scan})
  disp(sprintf('(projectOutMeanVectorPlot) Overlay does not exist for scan %i',scan));
  return
end

% look for this voxel in projection info
linearCoord = sub2ind(viewGet(view,'scanDims'),x,y,s);
thisLinearCoord = find(d{dnum}.linearCoords == linearCoord);

% get tSeries
tSeries = squeeze(loadTSeries(view,scan,s,frameNums,x,y));

% de mean normalize the tseries
tSeries = tSeries-mean(tSeries);
tSeries = tSeries/norm(tSeries);

% get projection magnitude
normMag = analysis.overlays(overlayNum).data{scan}(x,y,s);
reconMag = d{dnum}.reconProjectionMagnitude(thisLinearCoord);

% this voxel didn't have anything projected out
if isempty(thisLinearCoord)
  reconMag = 0;
end

% compute original
original = tSeries+reconMag*d{dnum}.sourceMeanVector';
original = original/norm(original);

% compute what is left of tSeries along projection vector
leftoverMag = tSeries'*d{dnum}.sourceMeanVector';

% select the window to plot into
selectGraphWin;

% plot the tSeries
subplot(2,1,1);
plot(tSeries,'k.-');hold on
if reconMag == 0
  plot(d{dnum}.sourceMeanVector,'r.-');
  legend('tSeries','Projection vector');
else
  plot(d{dnum}.sourceMeanVector*reconMag,'r.-');
  plot(original,'g.-');
  legend('tSeries','Projection vector','tSeries before projection');
end

xlabel('Time (TR)');
ylabel('Normalized magnitude');
title(sprintf('Voxel: [%i %i %i] Projection magnitude: %0.3f Projection ROI: %s\nMagnitude of component left in projection direction: %f',x,y,s,normMag,d{dnum}.sourceName,leftoverMag),'Interpreter','none');

subplot(2,1,2);
plot(abs(fft(tSeries)),'k.-');hold on
plot(abs(fft(d{dnum}.sourceMeanVector)),'r.-');
if reconMag ~= 0
  plot(abs(fft(original)),'g.-');
end
xaxis(0,size(tSeries,1)/2);
xlabel('Frequency (Fourier component number)')
ylabel('Magnitude');


