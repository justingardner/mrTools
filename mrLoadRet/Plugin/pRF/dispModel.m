function display = dispModel(voxels)

%
voxels = [47 74 24; 37 67 25];
analysis = 'pRF_lV1_123456.mat';
scanNum = 7;

% Set curr group to Concat and load the Analysis
v = newView;
v = viewSet(v, 'currentGroup', 'Concatenation');
v = loadAnalysis(v, ['pRFAnal/' analysis]);
a = viewGet(v, 'Analysis');
d = viewGet(v, 'd', scanNum);

% Create new figure
f = figure;

%%% (1) Voxel 1 Receptive Field %%%
rf1 = subplot(2, 3, 1);

  %Get params
x=voxels(1,1); y=voxels(1,2); z=voxels(1,3);
whichVoxel = find(d.linearCoords == sub2ind(viewGet(v, 'scanDims', scanNum), x, y, z));
r = d.r(whichVoxel, :);
params = d.params(:, whichVoxel);

  %get model fit for this voxel
m = pRFFit(v, scanNum, x,y,z, 'stim', d.stim, 'getModelResponse=1', 'params', params, 'concatInfo', d.concatInfo, 'fitTypeParams', a.params.pRFFit, 'paramsInfo', d.paramsInfo);

  %plot voxel RF
imagesc(d.stimX(:,1),d.stimY(1,:),flipud(m.rfModel'));
set(rf1,'Box','off');
set(rf1,'Color',[0.8 0.8 0.8]);
set(rf1,'TickDir','out');
axis equal
axis tight
hold on
hline(0,'w:');vline(0,'w:');
title('Voxel 1 Receptive Field Position');

%% (2) Voxel 2 Receptive Field
rf2 = subplot(2, 3, 4);
x=voxels(2,1); y=voxels(2,2); z=voxels(2,3);
whichVoxel = find(d.linearCoords == sub2ind(viewGet(v, 'scanDims', scanNum), x, y, z));
params = d.params(:, whichVoxel);
m2 = pRFFit(v, scanNum, x,y,z, 'stim', d.stim, 'getModelResponse=1', 'params', params, 'concatInfo', d.concatInfo, 'fitTypeParams', a.params.pRFFit, 'paramsInfo', d.paramsInfo);

imagesc(d.stimX(:,1),d.stimY(1,:),flipud(m2.rfModel'));
set(rf2,'Box','off');
set(rf2,'Color',[0.8 0.8 0.8]);
set(rf2,'TickDir','out');
axis equal
axis tight
hold on
hline(0, 'w:');vline(0,'w:');
title('Voxel 2 Receptive Field Position');

%% (3) Stimulus position

[thisT, modelT] = makeStim(100, 50);

function [thisT, modelT] = makeStim(time, modelTime)
  stP = subplot(2, 3, 2);
  %clf(stP);
  cla(stP);
  im = [];
  thisT = time;
  thisScan = d.concatInfo.whichScan(thisT);
  thisVolume = d.concatInfo.whichVolume(thisT);

  modelT = modelTime;
  modelScan = d.concatInfo.whichScan(modelT);
  modelVol = d.concatInfo.whichVolume(modelT);

  im(:, :, 3) = flipud(0.7*d.stim{thisScan}.im(:, :, thisVolume)');
  im(:, :, 2) = flipud(0.7*d.stim{modelScan}.im(:, :, modelVol)');
  im(:, :, 1) = 0.30*m.rfModel' + 0.30*m2.rfModel';
  image(d.stimX(:, 1), d.stimY(1, :), im);
  axis image
  hold on
  hline(0, 'w:'); vline(0, 'w:');
  title(sprintf('timeT=%i, modelT=%i', thisT, modelT));
end

  % Get covariance matrix
residual(1, :) = m.tSeries - m.modelResponse;
residual(2, :) = m2.tSeries - m2.modelResponse;
covMat = residual*residual.';
radius = covMat(1, 2);

%%%% (3) and (6): Plot voxels' model response time course and time series
% Plot voxel 1 model timecourse and timeseries
subplot(2, 3, 3);
plot((1:477), m.modelResponse, 'r');
hold on; plot((1:477), m.tSeries, 'black');
hold on; vl1m = vline(1); tl1m = text(0,0,'');
hold on; vl1t = vline(1); tl1t = text(0,0,'');
title('Voxel 1 Model Response and Time Series');

% Plot voxel 2 model timecourse and timeseries
subplot(2, 3, 6);
plot((1:477), m2.modelResponse, 'r');
hold on; plot((1:477), m2.tSeries, 'black');
hold on; vl2m = vline(1); tl2m = text(0,0,'');
hold on; vl2t = vline(1); tl2t = text(0,0,'');
title('Voxel 2 model response and time series');


%% (4) 2D Gaussian model + tSeries
p4 = subplot(2, 3, 5);

timepoint = 100;
modelTime = 50;

x1 = m.tSeries(timepoint);
y1 = m2.tSeries(timepoint);

x2 = m.modelResponse(modelTime);
y2 = m2.modelResponse(modelTime);

timePlot = plot(x1, y1, '*c');
hold on
cPlot = circle(x2, y2);
xlim([min(m.tSeries), max(m.tSeries)]);
ylim([min(m2.tSeries), max(m2.tSeries)]);
xlabel(sprintf('Voxel 1: (%d, %d, %d)', voxels(1,1), voxels(1,2), voxels(1,3)));
ylabel(sprintf('Voxel 2: (%d, %d, %d)', voxels(2,1), voxels(2,2), voxels(2,3)));
title('Voxel 1 vs Voxel 2: Percent Signal Change');

timeSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(m.tSeries),...
                       'SliderStep', [1/477, 10/477],...
                       'Value', timepoint, 'Callback', @slider1Callback,...
                       'Parent', f, 'Position', [.05, .05, .5, .1], 'Units', 'normalized');
function slider1Callback(hObject, eventdata)
  newVal = round(get(hObject, 'Value'));
  % Move the * around
  subplot(2, 3, 5);
  delete(timePlot);
  x1 = m.tSeries(newVal);
  y1 = m2.tSeries(newVal);
  timePlot = plot(x1, y1, '*b');

  % Move stim bar
  [thisT, modelT] = makeStim(newVal, modelT);

  % Move cyan line around time series course
  subplot(2, 3, 3);
  delete([vl1t tl1t]);
  vl1t = vline(newVal, '-b');
  tl1t = text(newVal, 0.92, sprintf('time t: %d', newVal), 'color', 'b');
  subplot(2, 3, 6);
  delete([vl2t tl2t]);
  vl2t = vline(newVal, '-b');
  tl2t = text(newVal, 0.96, sprintf('time t: %d', newVal), 'color', 'b');
end

modelSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(m.modelResponse),...
                        'SliderStep', [1/477 10/477], 'Value', modelTime, 'Callback', @slider2Callback,...
                        'Parent', f, 'Units', 'normalized', 'Position', [.05, .15, .5, .1]);
function slider2Callback(hObject, eventdata)
  newVal = round(get(hObject, 'Value'));
  % Move the circle around
  subplot(2, 3, 5);
  delete(cPlot);
  x2 = m.modelResponse(newVal);
  y2 = m2.modelResponse(newVal);
  cPlot = circle(x2, y2);

  % Move the stimulus bar
  [thisT, modelT] = makeStim(thisT, newVal);

  % Move green line around model time course
  subplot(2, 3, 3);
  delete([vl1m tl1m]);
  vl1m = vline(newVal, '-g');
  tl1m = text(newVal, 0.95, sprintf('Model t: %d', newVal), 'color', 'g');
  subplot(2, 3, 6);
  delete([vl2m tl2m]);
  vl2m = vline(newVal, '-g');
  tl2m = text(newVal, 0.97, sprintf('Model t: %d', newVal), 'color', 'g');
end 

align([timeSlider, modelSlider, f], 'center', 'distribute');


function cPlot = circle(x, y)
  r = radius;
  ang = 0:0.01:2*pi;
  xp = r*cos(ang);
  yp = r*sin(ang);
  cPlot = plot(x+xp, y+yp, 'g');
end

end
