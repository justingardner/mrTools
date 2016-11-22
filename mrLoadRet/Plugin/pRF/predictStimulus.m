%%  predictStimulus.m
%%
%%
%%      usage: predictStimulus(fits, d) --> plot figs for avg (last element of fits)
%%             predictStimulus(fits, d, index) --> plot figs for i'th fold in fits
%%         by: akshay jagadeesh
%%       date: 11/20/2016
%%    purpose: Given model fits, this function predicts the stimulus most likely to have elicited
%%             the given model response and compares this to the actual observed stimulus. Outputs the
%%             percentage of timepoints that are correctly predicted (within 5) using maximum likelihood. 
%%


function predictStimulus(fits, d, index)

numFolds = length(fits)-1;
if(ieNotDefined('index'))
  fit = fits(numFolds+1);
else
  disp(sprintf('getting %d th fit', index)); 
  fit = fits(index);
end

timeLen = size(fit.modelResponse, 2);

%%% 1. Plot prediction probability chart (batik)
f=figure; subplot(2, 2, 1); imagesc(log(fit.probTable)); axis ij; colorbar; 
xlabel('Model Time Series'); ylabel('Observed Time series');
hold on; hl = hline(200, '-b');
hold on; pt = plot(0,0,'*g');
    % Plot best predicted points
count = 0;
for i = 1:timeLen
  tS = fit.probTable(i,:);
  [~, mI] = max(tS);
  maxInd(i) = mI;
  if( abs(mI - i) <= 5)
    count = count + 1;
  end
end
hold on; plot(maxInd, (1:timeLen), '*g', 'MarkerSize', 1);
title(sprintf('Avg Model Prediction Probability:\n %d / %d timepoints (%0.1f %%) correctly predicted within 5 using Max Likelihood', count, timeLen, 100*count/timeLen));

%%% 3. Plot observed and model response time course
subplot(2, 2, 3); plot((1:timeLen), fit.leftOut, 'k');
hold on; plot((1:timeLen), fit.modelResponse, 'r'); xlim([1 timeLen]); title('Left out time series + Model Response');
hold on; vl1 = vline(0, 'k'); hold on; vl2 = vline(0, 'k');

%%% 4. Plot voxel receptive fields
%rf = subplot(2,2,4); plot(fits(1).fitParams(:,4), fits(1).fitParams(:,5), '*');
%title('Voxel Receptive Field positions'); 
%set(rf, 'Box', 'off'); set(rf, 'Color', [.8 .8 .8]); set(rf,'TickDir','out');
%axis equal; axis tight; xlim([-30 30]); ylim([-20 20]); hold on; hline(0, 'w:'); vline(0,'w:');

timeSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(fit.modelResponse),...
                       'SliderStep', [1/timeLen 10/timeLen], 'Value', 200, 'Callback', @timeSliderCallback,...
                       'Parent',f, 'Units', 'normalized', 'Position', [.25, -.05,.5,.1]);
function timeSliderCallback(hObject, evendata)
  newVal = round(get(hObject, 'Value'));

  %plot the stimulus predicted by the model time point with highest posterior probability
  timeSeries = fit.probTable(newVal, :);
  [maxVal, maxIndex] = max(timeSeries);
  stimPl = subplot(2, 2, [2 4]);
  cla(stimPl);
  im = [];
  predScan = d.concatInfo.whichScan(maxIndex);
  predVol = d.concatInfo.whichVolume(maxIndex);
  trueScan = d.concatInfo.whichScan(newVal);
  trueVol = d.concatInfo.whichVolume(newVal);
  im(:, :, 3) = flipud(0.7*d.stim{trueScan}.im(:,:,trueVol)');
  im(:, :, 2) = flipud(0.7*d.stim{predScan}.im(:,:,predVol)');
  im(:, :, 1) = 0; 
  image(d.stimX(:,1), d.stimY(1,:), im);
  axis image; hold on;
  hline(0, 'w:'); vline(0, 'w:');
  hold on; plot(fits(1).fitParams(:, 4), fits(1).fitParams(:, 5), '*w', 'MarkerSize', 1);
  title(sprintf('Stim predictions: Timepoint %i is best predicted by Modelpoint %i', newVal, maxIndex));
  
  % Plot line and point of highest prediction
  subplot(2,2,1);
  delete([hl pt]);
  hl = hline(newVal, '-b');
  set(hl, 'LineWidth', 5);
  pt = plot(maxIndex, newVal, '*g');
  
  % Plot vertical line on time series plots
  subplot(2, 2, 3);
  delete([vl1 vl2]);
  vl1 = vline(newVal, '-b');
  set(vl1, 'LineWidth', 5);
  hold on;
  vl2 = vline(maxIndex, '-g');
  set(vl2, 'LineWidth', 5);
  title(sprintf('Time Series: %d; Model Time: %d', newVal, maxIndex));
end

bl3 = uicontrol('Parent',f,'Style','text','Units', 'normalized', 'Position', [.45, .05 ,.1,.015],...
                'String', sprintf('Time Slider (blue)'));

end
