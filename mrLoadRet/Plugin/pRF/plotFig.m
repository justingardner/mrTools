% plotFig
%
%      does analysis given struct of fits for two models
%      - fitsGauss, fitsNorm: struct 
%
%
function plotFig(fitsGauss, fitsNorm, sim, num, stim, concatInfo)

if any(num==1) % draw true rf's
  dispTrueRFs = 1;
end
if any(num==1.5)
  linkRFs = 1;
end
if any(num==2) % draw model predicted rf's
  dispModelRFs = 1;
end
if any(num==3)
  dispModelR2 = 1;
end
if any(num==4)
  drawR2 = 1;
end
if any(num==5)
  dispBestModelResponse = 1;
end
if any(num==6)
  dispSuppFieldSize = 1; 
end
if any(num==7)
  dispWorstModelResponse = 1;
end

paramsG = getFieldArrFromStruct(fitsGauss, 'params');
paramsN = getFieldArrFromStruct(fitsNorm, 'params');
r2G = getFieldArrFromStruct(fitsGauss, 'r2');
r2N = getFieldArrFromStruct(fitsNorm, 'r2');

if ~ieNotDefined('dispTrueRFs') 
  figure;
  for i = 1:196
    simRFWidth = sqrt(sim.x(i)^2 + sim.y(i)^2)*0.1;
    circle(sim.x(i), sim.y(i), simRFWidth, 'b'); hold on;
    text(sim.x(i), sim.y(i), num2str(i));
  end
  title('True Voxel Receptive Fields (Simulation)');
  hline(0, ':'); vline(0,':');
  
end
if ~ieNotDefined('linkRFs')
  figure;
  plot(sim.x, sim.y, '*k'); hold on;
  plot(paramsG(:,1), paramsG(:,2), '*b'); hold on
  for i = 1:196
    plot([paramsG(i,1); sim.x(i)], [paramsG(i,2); sim.y(i)], 'k'); hold on
  end
  legend('True RFs', 'Gaussian RFs'); 
  title('Mapping of Gaussian Model RFs to true RFs'); hline(0,':'); vline(0,':');

  figure;
  plot(sim.x, sim.y, '*k'); hold on;
  plot(paramsN(:,1), paramsN(:,2), '*g'); hold on
  for i = 1:196
    plot([paramsN(i,1); sim.x(i)], [paramsN(i,2); sim.y(i)], 'k'); hold on
  end
  legend('True RFs', 'Normalization RFs');
  title('Mapping of Normalization Model RFs to true RFs'); hline(0,':'); vline(0,':');
end
if ~ieNotDefined('dispModelRFs')
  figure;
  for i = 1:196
    rfWidth = sqrt(paramsG(i,1)^2 + paramsG(i,2)^2)*0.1;
    circle(paramsG(i,1), paramsG(i,2), paramsG(i,3), 'b'); hold on;
    text(paramsG(i,1), paramsG(i,2), num2str(i));
  end
  title('Gaussian Model Predicted Receptive Fields');
  hline(0, ':'); vline(0,':');

  figure;
  for i = 1:196
    rfWidth = sqrt(paramsN(i,1)^2 + paramsN(i,2)^2)*0.1;
    circle(paramsN(i,1), paramsN(i,2), paramsN(i,3), 'b');
    circle(paramsN(i,1), paramsN(i,2), paramsN(i,3)*paramsN(i,4), 'g');
    text(paramsN(i,1), paramsN(i,2), num2str(i));
  end
  title('Normalization Model Predicted Receptive & Suppressive Fields');
  hline(0,':'), vline(0,':');
end
if ~ieNotDefined('dispModelR2')
  % Plot model RF predictions colored by r2
  figure;
  scatter(paramsG(:,1), paramsG(:,2), 20, r2G);
  title('Gaussian model RF predictions, colored by R2');
  hline(0,':'); vline(0,':'); colorbar;

  figure;
  scatter(paramsN(:,1), paramsN(:,2), 20, r2N);
  title('Normalization model RF predictions, colored by R2');
  hline(0,':'); vline(0,':'); colorbar;
end
if ~ieNotDefined('drawR2')
  % Draw R2 by voxel for each model
  figure;
  plot(r2G, 'b'); hold on; plot(r2N, 'g');
  title('Model r2 plotted by voxel');
  xlabel('Voxel Number'); ylabel('R-squared');
end
if ~ieNotDefined('dispBestModelResponse') || ~ieNotDefined('dispWorstModelResponse')
  if ieNotDefined('stim') || ieNotDefined('concatInfo')
    disp('Either stimulus or concatInfo not provided. Quitting...');
    return
  end
  
  [bestR2N bestVoxN]  = max(r2N); bestN = fitsNorm{bestVoxN};
  [bestR2G bestVoxG] = max(r2G); bestG = fitsGauss{bestVoxG};
  
  G = @(x,y,sigma) exp(-(((stim.x - x).^2) + ((stim.y - y).^2))/(2*sigma^2));
  respNorm = @(c, x, y, sigma, stdRatio, b1, b2) (b1*c*G(x,y,sigma))./(1+ b2*c*G(x,y,stdRatio*sigma));
  Rx = 0:.05:1;
  Ry = [];

  modelResponseG = [];
  modelResponseN = [];
  tSeries = [];
  hrf = getCanonicalHRF();
  disppercent(-inf, 'Computing model response based on model params for 196 voxels');
  for i = 1:196
    gParams(i,:) = fitsGauss{i}.params;
    rfModel = exp(-(((stim.x - gParams(i,1)).^2) + ((stim.y - gParams(i,2)).^2))/(2*(gParams(i,3)^2)));
    rfModel = convolveModelWithStimulus(rfModel, stim.im);
    modelGwithHRF = convolveModelResponseWithHRF(rfModel, hrf);
    modelResponseG(i,:) = applyConcatFiltering(modelGwithHRF/std(modelGwithHRF), concatInfo);

    nParams(i,:) = fitsNorm{i}.params;
    rfModel = exp(-(((stim.x - nParams(i,1)).^2) + ((stim.y - nParams(i,2)).^2))/(2*(nParams(i,3)^2)));
    rfModel = convolveModelWithStimulus(rfModel, stim.im);
    rfModel2 = exp(-(((stim.x - nParams(i,1)).^2) + ((stim.y-nParams(i,2)).^2))/(2*((nParams(i,4)*nParams(i,3))^2)));
    rfModel2 = convolveModelWithStimulus(rfModel2, stim.im);
    modelNwithHRF = convolveModelResponseWithHRF((nParams(i,5)*rfModel)./(nParams(i,6)*rfModel2), hrf);
    modelResponseN(i,:) = applyConcatFiltering(modelNwithHRF, concatInfo);

    tSeries(i,:) = applyConcatFiltering(sim.tSeries(i,:), concatInfo);
    disppercent(i/196);
  end
  disppercent(inf);

  if ~ieNotDefined('dispBestModelResponse')
  figure;
  plot(sim.filteredTSeries(bestVoxN,:), 'k'); hold on;
  plot(modelResponseN(bestVoxN,:), 'g'); hold on;
  plot(modelResponseG(bestVoxN,:), 'b');
  title('Model response and true time series for well-fit voxel');
  xlabel('Time (volumes)'); ylabel('Response');
  legend('timeseries', 'normalization', 'gaussian');
  end
  if ~ieNotDefined('dispWorstModelResponse')
    figure;
    plot(sim.filteredTSeries(194,:), 'k'); hold on;
    plot(modelResponseN(194,:), 'g'); hold on;
    title(sprintf('Voxel 194 model response vs time series (R2: %d)', r2N(194)));
    xlabel('Time (volumes)'); ylabel('Response');
    legend('timeseries', 'normalization');

    figure; 
    plot(sim.filteredTSeries(189,:), 'k'); hold on;
    plot(modelResponseN(189,:), 'g'); hold on;
    title(sprintf('Voxel 189 model response vs time series (R2:%d)', r2N(189)));
    xlabel('Time (volumes)'); ylabel('Response');
    legend('timeseries', 'normalization'); 

    figure;
    plot(sim.filteredTSeries(1,:), 'k'); hold on;
    plot(modelResponseN(1,:), 'g'); hold on;
    title(sprintf('Voxel 1 model response vs time series (R2:%d)', r2N(1)));
    xlabel('Time (volumes)'); ylabel('Response');
    legend('timeseries', 'normalization');
  end
end  
if ~ieNotDefined('dispSuppFieldSize')

  figure;
  rfWidthN = paramsN(:,3);
  suppWidthN = paramsN(:,4).*rfWidthN;
  eccN = sqrt(paramsN(:,1).^2 + paramsN(:,2).^2);
  plot(eccN, suppWidthN, '*g'); hold on;
  fit = polyfit(eccN, suppWidthN, 1);
  line = @(x) fit(1)*x + fit(2);
  plot([min(eccN), max(eccN)], [line(min(eccN)) line(max(eccN))], 'k'); text(5, 7, sprintf('Slope: %d', fit(1)));
  title('Relationship between suppressive field size and eccentricity');
  xlabel('Eccentricity'); ylabel('Suppressive field width');

  figure; 
  plot(eccN, rfWidthN, '*g'); hold on;
  fit = polyfit(eccN, rfWidthN, 1); line = @(x) fit(1)*x + fit(2);
  plot([min(eccN), max(eccN)], [line(min(eccN)) line(max(eccN))], 'k'); text(5, 7, sprintf('Slope: %d', fit(1)));
  title('Relationship between RF field size and eccentricity (normalization model)');
  xlabel('Eccentricity'); ylabel('RF Field Width');

  figure;
  eccG = sqrt(paramsG(:,1).^2 + paramsG(:,2).^2);
  rfWidthG = paramsG(:,3);
  plot(eccG, rfWidthG, '*b'); hold on;
  fit = polyfit(eccG, rfWidthG, 1); line = @(x) fit(1)*x + fit(2);
  plot([min(eccN), max(eccG)], [line(min(eccG)) line(max(eccG))], 'k'); text(5, 7, sprintf('Slope: %d', fit(1)));
  title('Relationship between RF field size and eccentricity (gaussian model)');
  xlabel('Eccentricity'); ylabel('RF Field Width');

  keyboard
end
if num == -1
  %%%%
  % Plot model RFs and true RFs using circle tool
  f = figure;
  for i = 1:196
    %% Plot true RF's
    plot(sim.x(i), sim.y(i), '*r', 'MarkerSize', 10); hold on;
    %simRFWidth = sqrt(sim.x(i)^2 + sim.y(i)^2)*0.1;
    %circle(sim.x(i), sim.y(i), simRFWidth, 'r'); hold on;

    %% Plot gaussian model RF's
    fit = fitsGauss{i};
    %circle(fit.x, fit.y, fit.std, 'b');

    %% Plot normalization model rf's
    fit2 = fitsNorm{i};
    circle(fit2.x, fit2.y, fit2.std, 'g'); 
    circle(fit2.x, fit2.y, fit2.params(4)*fit2.std, 'c')

    %% Plot Lines connecting model rf center to true rf center
    %plot([fit.x; sim.x(i)], [fit.y; sim.y(i)],'-b','MarkerSize',3); hold on
    plot([fit2.x; sim.x(i)], [fit2.y; sim.y(i)], '-g', 'MarkerSize', 3); hold on
  end
  xlim([-100 100]); ylim([-100 100]);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%% < end of main section > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   helper methods below  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------
% getFieldArrFromStruct(structure, fieldNameStr, 
%--------------------------------------
function fieldArray = getFieldArrFromStruct(structure, fieldNameStr)

fieldArray = [];
for i = 1:length(structure)
  fieldArray(i,:) = structure{i}.(fieldNameStr);
end

%---------------------------------------
% circle
%     draw a circle given center (x,y), radius and color
%----------------------------------------
function circle(x,y,r, color)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp, color, 'MarkerSize', 4);
hold on


%---------------------------------------
% getCanonicalHRF()
%     get the canonical hemodnamic response function 
%----------------------------------------
function hrf = getCanonicalHRF()
offset = 0;
timelag = 1;
tau = 0.6;
exponent = 6;
sampleRate = 0.5;
amplitude = 1;
hrf.time = 0:sampleRate:25;

exponent = round(exponent);
gammafun = (((hrf.time - timelag)/tau).^(exponent-1).*exp(-(hrf.time-timelag)/tau))./(tau*factorial(exponent-1));
gammafun(find((hrf.time-timelag) <0)) = 0;

if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);

hrf.hrf = gammafun;
hrf.hrf = hrf.hrf / max(hrf.hrf);


%---------------------------------------
% convolveModelWithStimulus 
%      convolve the model response with the stimulus
%----------------------------------------
function modelResponse = convolveModelWithStimulus(rfModel, stim)
nStimFrames = size(stim, 3);
modelResponse = zeros(1,nStimFrames);
for frameNum = 1:nStimFrames
  modelResponse(frameNum) = sum(sum(rfModel.*stim(:,:,frameNum)));
end

%---------------------------------------
% convolveModelResponseWithHRF 
%      convolve the model response with the HRF
%----------------------------------------
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)
n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

%---------------------------------------
% applyConcatFiltering
%      applies detrending, hipassfilter, projection, removes mean.
%----------------------------------------
function tSeries = applyConcatFiltering(tSeries, concatInfo)
runnum=1;
tSeries = tSeries(:);

% apply detrending

if ~isfield(concatInfo,'filterType') || ~isempty(findstr('detrend',lower(concatInfo.filterType)))
  tSeries = eventRelatedDetrend(tSeries);
end

% apply hipass filter
if isfield(concatInfo,'hipassfilter') && ~isempty(concatInfo.hipassfilter{runnum})
  if ~isequal(length(tSeries),length(concatInfo.hipassfilter{runnum}))
          disp(sprintf('(applyConcatFiltering) Mismatch dimensions of tSeries (length: %i) and concat filter (length: %i)',length(tSeries),length(concatInfo.hipassfilter{runnum})));
  else
    tSeries = real(ifft(fft(tSeries) .* repmat(concatInfo.hipassfilter{runnum}', 1, size(tSeries,2)) ));
  end
end

tSeries = tSeries-repmat(mean(tSeries,1),size(tSeries,1),1);
tSeries = tSeries(:)';
