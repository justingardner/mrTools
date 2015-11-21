% numFit.m
%
%        $Id:$ 
%      usage: numFit()
%         by: Guillaume Riesen / Justin Gardner
%       date: 11/20/15
%    purpose: Fit numerosity model to a voxel
%
function retval = numFit(v,overlayNum,scan,x,y,s,roi)

% get time series
tSeries = squeeze(loadTSeries(v,scan,s,[],x,y));tSeries = tSeries(:);
concatInfo = viewGet(v,'concatInfo');

% get stimvols
stimvol = getStimvol(v,'_all_','taskNum=1','phaseNum=2');

% now average over trials
trialVols = stimvol{1};
trialLen = min(diff(trialVols));
nTrials = length(trialVols)-1;
for iTrials = 1:nTrials
  meanTSeries(iTrials,1:trialLen) = tSeries(trialVols(iTrials):trialVols(iTrials)+trialLen-1);
end
tSeries = mean(meanTSeries,1);

% percent signal change
tSeries = (tSeries-mean(tSeries))/mean(tSeries);

% get basic paramters (like for canonical HDR)
[v fitParams] = pRF(v,[],'justGetParams=1','defaultParams=1','scanList',scan);
fitParams = fitParams.pRFFit;

% set the inital parameters to start off the model with
fitParams.initParams = [1 3];

% set algorithm for nonlinear optimization and some parameters
fitParams.algorithm = 'nelder-mead';
fitParams.optimParams = optimset('MaxIter',inf);
fitParams.framePeriod = viewGet(v,'framePeriod');

% FIX, need to set this to the stimulus sequence
fitParams.stim =  []

% Just for testing
r = getModelResidual(fitParams.initParams,tSeries,fitParams);
keyboard

% find best parameters
[params fval exitflag] = fminsearch(@getModelResidual,fitParams.initParams,fitParams.optimParams,tSeries,fitParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getModelResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual modelResponse rfModel r] = getModelResidual(params,tSeries,fitParams,justGetModel)

residual = [];
if nargin < 4, justGetModel = 0;end

% get the model response
% convert parameter array into a parameter strucutre
p = getFitParams(params,fitParams);

% FIX, compute the RF
%rfModel = getRFModel(p,fitParams);

% FIX, compute the RF response to the stimulus (using fitParams.stim which should have the stimulus sequence)
modelResponse = zeros(1,length(tSeries));
%modelResponse = computeModelResponse(rfModel, fitParams.stim,...)

% get a model hrf
hrf = getCanonicalHRF(p.canonical,fitParams.framePeriod);

% and convolve in time.
modelResponse = convolveModelResponseWithHRF(modelResponse,hrf);

% with no filtering, just remove mean
modelResponse = modelResponse - mean(modelResponse);

% compute correlation of model with time series
r = corr(tSeries(:),modelResponse(:));

dispFit = 1;
if dispFit
  mlrSmartfig('numFit','reuse');clf;
  plot(100*tSeries,'k.-');hold on;
  plot(100*modelResponse,'r-');
  xlabel('Time (vols)');
  ylabel('BOLD (% signal change)');
  title(sprintf('Center: %f Width: %f',p.center,p.width));
  drawnow
end
  
% for nelder-mead just compute correlation and return negative
% since nelder-mean is optimizing for the least value
if strcmp(lower(fitParams.algorithm),'nelder-mead')
  residual = -r;
end


%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)

p.center = params(1);
p.width = params(2);

% use a fixed single gaussian
p.canonical.type = 'gamma';
p.canonical.lengthInSeconds = 25;
p.canonical.timelag = fitParams.timelag;
p.canonical.tau = fitParams.tau;
p.canonical.exponent = fitParams.exponent;
p.canonical.offset = 0;
p.canonical.diffOfGamma = fitParams.diffOfGamma;
p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
p.canonical.timelag2 = fitParams.timelag2;
p.canonical.tau2 = fitParams.tau2;
p.canonical.exponent2 = fitParams.exponent2;
p.canonical.offset2 = 0;

%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

