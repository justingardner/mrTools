% pRFPreProcess.m
%
%        $Id:$ 
%      usage: pRFPreProcess()
%         by: justin gardner
%       date: 11/30/15
%    purpose: Preprocess model time series - convolve with
%             hemodynamic response, apply concat filtering and apply beta
%             weights
%
function [response p hrf] = pRFPreProcess(params,fitParams,response,tSeries)

% FIX, FIX, FIX - concatInfo stuff and beta each scan stuff
% also handle short time series

% check here for response shorter then expected
if length(response) ~= fitParams.concatInfo.runTransition(end,end)
  % make up a concat info for this length run 
  concatInfo.n = 1;
  concatInfo.runTransition = [1 length(response)];
  concatInfo.totalJunkedFrames = 0;
  fitParams.applyFiltering = false;
else
  % normal processing use real concatInfo
  concatInfo = fitParams.concatInfo;
end

% get the parameters for the hemodynamic function
p = getHemodynamicFitParams(params,fitParams);

% for each scan in the concat
for iScan = 1:concatInfo.n
  % get this scan of the response
  thisResponse = response(concatInfo.runTransition(iScan,1):concatInfo.runTransition(iScan,2));

  % get a model hrf
  hrf = getCanonicalHRF(p,fitParams.framePeriod);

  % and convolve in time.
  thisModelResponse = convolveModelResponseWithHRF(thisResponse,hrf);

  % drop junk frames here
  thisResponse = thisResponse(concatInfo.totalJunkedFrames(iScan)+1:end);

  % apply concat filtering
  if fitParams.applyFiltering
    thisResponse = applyConcatFiltering(thisResponse,concatInfo,iScan);
  else
    % with no filtering, just remove mean
    thisResponse = thisResponse - mean(thisResponse);
  end

  % If we have been passed in tSeries and need to scale
  % each scan, then do that here.
  if (nargin>=4) && fitParams.betaEachScan
    % scale and offset the model to best match the tSeries
    [thisResponse thisResidual] = scaleAndOffset(thisResponse,thisTSeries);
  else
    thisResidual = [];
  end
  
  % put back into response
  response(concatInfo.runTransition(iScan,1):concatInfo.runTransition(iScan,2)) = thisResponse;
end

% normalize
response = response(:)'-mean(response);
response = response/sqrt(response*response');

% scale the whole time series if we have been passed in tSeries
if (nargin>=4) && ~fitParams.betaEachScan
  [response residual] = scaleAndOffset(response,tSeries);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getHemodynamicFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = getHemodynamicFitParams(params,fitParams)

if fitParams.fitHemodynamic
  % get the hemodynamic parameters
  hemoParams = params(fitParams.nParams+1:end);
  % canonical is being fit
  p.lengthInSeconds = 25;
  p.timelag = hemoParams(1);
  p.tau = hemoParams(2);
  p.exponent = fitParams.exponent;
  p.offset = 0;
  p.diffOfGamma = fitParams.diffOfGamma;
  if fitParams.diffOfGamma
    p.amplitudeRatio = hemoParams(3);
    p.timelag2 = hemoParams(4);
    p.tau2 = hemoParams(5);
    p.exponent2 = fitParams.exponent2;
    p.offset2 = 0;
  end
else
  % use a fixed single gaussian
  p.lengthInSeconds = 25;
  p.timelag = fitParams.timelag;
  p.tau = fitParams.tau;
  p.exponent = fitParams.exponent;
  p.offset = 0;
  p.diffOfGamma = fitParams.diffOfGamma;
  p.amplitudeRatio = fitParams.amplitudeRatio;
  p.timelag2 = fitParams.timelag2;
  p.tau2 = fitParams.tau2;
  p.exponent2 = fitParams.exponent2;
  p.offset2 = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applyConcatFiltering    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tSeries = applyConcatFiltering(tSeries,concatInfo,runnum)

% apply the same filter as original data
% check for what filtering was done
tSeries = tSeries(:);

% apply detrending (either if concatInfo does not say what it did or if
% the filterType field has detrend in it)
if ~isfield(concatInfo,'filterType') || ~isempty(findstr('detrend',lower(concatInfo.filterType)))
  tSeries = eventRelatedDetrend(tSeries);
end

% apply hipass filter
if isfield(concatInfo,'hipassfilter') && ~isempty(concatInfo.hipassfilter{runnum})
  % check for length match
  if ~isequal(length(tSeries),length(concatInfo.hipassfilter{runnum}))
    disp(sprintf('(pRFFit:applyConcatFiltering) Mismatch dimensions of tSeries (length: %i) and concat filter (length: %i)',length(tSeries),length(concatInfo.hipassfilter{runnum})));
  else
    tSeries = real(ifft(fft(tSeries) .* repmat(concatInfo.hipassfilter{runnum}', 1, size(tSeries,2)) ));
  end
end

% project out the mean vector
if isfield(concatInfo,'projection') && ~isempty(concatInfo.projection{runnum})
  projectionWeight = concatInfo.projection{runnum}.sourceMeanVector * tSeries;
  tSeries = tSeries - concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
end

% now remove mean
tSeries = tSeries-repmat(mean(tSeries,1),size(tSeries,1),1);

% make back into the right dimensions
tSeries = tSeries(:)';

%%%%%%%%%%%%%%%%%%%%%%%%
%    scaleAndOffset    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [modelResponse residual] = scaleAndOffset(modelResponse,tSeries)

modelResponse = modelResponse(:);
tSeries = tSeries(:);

designMatrix = modelResponse;
designMatrix(:,2) = 1;

% get beta weight for the modelResponse
if ~any(isnan(modelResponse))
  beta = pinv(designMatrix)*tSeries;
  beta(1) = max(beta(1),0);
  modelResponse = designMatrix*beta;
  residual = tSeries-modelResponse;
else
  residual = tSeries;
end

