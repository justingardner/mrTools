% pRFSimulate
%
%    simulates N v1 voxels given 2 rf types and outputs expected time series given stimImage
%
%    - stim: struct containing fields x (192x108), y (192x108), and im (192x108x time)
%    - numVox: number of voxels we want to simulate
%
function [simulation, fits_gauss,fits_norm, r2_gauss, r2_norm] = pRFSimulate(stim)

% stimImage = stim.im;
% numVols = size(stimImage, 2);
% [myscreen, stimImage] = 
modelType = 'gaussian-divs';
plotFigs = 0;

%%% Preset Model Parameters %%%
modelParams.rfRatio = 0.1;
modelParams.normRatio = 1.5; % set norm pool width to be 1.5x the rf width
modelParams.cssExp = 0.75;
modelParams.noiseAmplitude = 0.1;

% Get canonical hemodynamic response function
hrf = getCanonicalHRF();

%%% Create voxels in simulation
% Num voxels: 196 --> same as in prefit
xVals = -32.5:5:32.5; xVals = repmat(xVals, 1, 14);
yVals = -32.5:5:32.5; yVals = repmat(yVals, 14, 1); yVals = yVals(:);

if(length(xVals) == length(yVals))
  numVox = length(xVals);
  disp(sprintf('Number of voxels in simulation: %d', numVox));
end
if(plotFigs == 1); plot(xVals, yVals, '*'); end

simulation = [];
for i=1:numVox

  x_sim = xVals(i); y_sim = yVals(i);
  eccentricity = sqrt(x_sim^2 + y_sim^2);
  rfWidth = modelParams.rfRatio*eccentricity;
  
  % Calculate the Gaussian-derived RF from Dumoulin & Wandell 2008
  rfModel = exp(-(((stim.x - x_sim).^2) + ((stim.y - y_sim).^2))/(2*(rfWidth^2)));
  if(mod(i,10)==0 && plotFigs == 1)
    disp(i); figure; 
    subplot(1,2,1); plot(x_sim, y_sim, '*'); xlim([-50,50]); ylim([-50,50]);
    subplot(1,2,2); imagesc(rfModel'); axis xy; xlim([0,192]); ylim([0,108]); 
  end

  % Calculate simulated time series using the Ratio-of-Gaussians model
  resp = convolveModelWithStimulus(rfModel, stim.im);
  suppFieldWidth = modelParams.normRatio*rfWidth;
  beta1 = 1; beta2 = 1;
  % Calculate suppressive field model
  rfModel2 = exp(-(((stim.x - x_sim).^2) + ((stim.y - y_sim).^2))/(2*(suppFieldWidth^2)));
  resp2 = convolveModelWithStimulus(rfModel2, stim.im);
  % Estimate model response as the ratio of the 2 gaussian RF
  modelResponse = beta1*resp ./ (1 + beta2*resp2);

  % Convolve with HRF 
  modelHRF = convolveModelResponseWithHRF(modelResponse, hrf);

  % Then add white noise and save the simulated time series data.
  simulation(i,:) = addWhiteNoise(modelHRF, modelParams.noiseAmplitude);
  %noiseless(i,:) = modelHRF;
  %gaussian(i,:) = convolveModelResponseWithHRF(resp, hrf);
  %simulation(i,:) = addWhiteNoise(convolveModelResponseWithHRF(resp./(max(resp) - min(resp)), hrf));
  %simulation(i,:) = convolveModelResponseWithHRF(resp./(max(resp) - min(resp)),hrf);
end

%%% To Do:
% 1. Compute stimulus image with 464 volumes (equal to s0315 retinotopy)
%     a. DONE -- figure out why mglDoubleBars framegrab is acting up.
% 2. Figure out how to map the huge number of prefit computed voxels to my simulated voxels
%     -- Prefit computes 7623 (33x33x7) voxels using [prefitx prefity prefitrfHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
%     so just map my simulated voxels to values from -0.4 to +0.4 in intervals of 0.025

%%% Model Fit: now fit the pRF model to the simulation using 2 different models (rftypes)
v = newView;
v = viewSet(v, 'currentGroup', 'Concatenation');
v = loadAnalysis(v, ['pRFAnal/' 'pRF_v1.mat']);
analParams = v.analyses{1}.params;
analParams.pRFFit.quickPrefit=1; % Use quick prefit for 196 voxels instead of 7623
analParams.pRFFit.verbose = 0; %set verbose to 0 to silence all that annoying printing
analParams.pRFFit.prefitOnly=1;
analParams2 = analParams;
analParams2.pRFFit.rfType='gaussian-divs';
d = viewGet(v, 'd', analParams.scanNum);
concatInfo = viewGet(v, 'concatInfo', 1);
v = newView;

stimA{1} = stim;
numVoxels = numVox;

% Compute prefit for all voxels
disppercent(-inf, sprintf('(pRFSimulate) Computing fits for %d simulated voxels', numVoxels));
for i = 1:numVoxels
  tSeries = applyConcatFiltering(simulation(i,:), concatInfo);
  disp(sprintf('Voxel %d', i));

  % Run for gaussian model
  fit = pRFFit(v, [], 1, [], [], 'stim', stimA, 'tSeries', tSeries, 'concatInfo', concatInfo, 'fitTypeParams', analParams.pRFFit, 'quickPrefit', true, 'verbose=0');
  r2_gauss(i) = fit.r2;
  fits_gauss{i} = fit;

  % Run for normalization model
  fit2 = pRFFit(v, [], 1, [], [], 'stim', stimA, 'tSeries', tSeries, 'concatInfo', concatInfo, 'fitTypeParams', analParams2.pRFFit, 'quickPrefit', true, 'verbose=0');
  r2_norm(i) = fit2.r2;
  fits_norm{i} = fit2;
  disppercent(i/numVoxels);
end
disppercent(inf);
disp(sprintf('(pRFSimulate) Successfully fit model to %d simulated voxels', numVoxels));

%keyboard

return

%% Do full pRFFits
%for i = 1:numVoxels
%  fit = pRFFit(v, [], [], [], [], 'fitTypeParams', analParams.pRFFit, 'returnPrefit', true);
%  keyboard
%  modelFit1(i,:) = pRFFit(v, [], [],[],[], 'stim', stim, 'tSeries', simulation(i,:), 'getModelResponse=1', 'rfType=gaussian',...
%                          'fitTypeParams', analParams.pRFFit, 'paramsInfo', d.paramsInfo);
%  modelFit2(i,:) = pRFFit(v, [], [],[],[], 'stim', stim, 'tSeries', simulation(i,:), 'getModelResponse=1', 'rfType=gaussian_divs',...
%                          'fitTypeParams', analParams.pRFFit, 'paramsInfo', d.paramsInfo);
%end

%--------------------------------------------
%      addWhiteNoise
% adds gaussian white noise to the signal
function response = addWhiteNoise(signal, noiseAmplitude)
response=signal+sqrt(noiseAmplitude).*randn(1,size(signal,2));

%----------------------------------------------
%           getCanonicalHRF
%   Gets the canonical hrf given params and sample rate
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    ~~~ convolveModelWithStimulus ~~~    %%

function modelResponse = convolveModelWithStimulus(rfModel, stim)

nStimFrames = size(stim, 3);
modelResponse = zeros(1,nStimFrames);
for frameNum = 1:nStimFrames
  modelResponse(frameNum) = sum(sum(rfModel.*stim(:,:,frameNum)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% project out the mean vector
