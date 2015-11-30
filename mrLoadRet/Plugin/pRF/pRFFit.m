% pRFFit
%
%      usage: pRFFit(v,scanNum,x,y,s,<dispFit=0>)
%         by: justin gardner
%       date: 11/14/2011
%    purpose: interrogator that fits pRF model to selected voxel
%
function fit = pRFFit(varargin)

fit = [];
% parse input arguments - note that this is set
% up so that it can also be called as an interrogator
[v scanNum x y z fitParams tSeries] = parseArgs(varargin);
if isempty(v),return,end
 
% get the tSeries
if ~isempty(x)
  % if tSeries was not passed in then load it
  if isempty(tSeries)
    % load using loadTSeries
    tSeries = squeeze(loadTSeries(v,scanNum,z,[],x,y));
  end

  % convert to percent tSeries. Note that we  detrend here which is not necessary for concats,
  % but useful for raw/motionCorrected time series. Also, it is very important that
  % the tSeries is properly mean subtracted
  if ~isfield(fitParams.concatInfo,'hipassfilter')
    tSeries = percentTSeries(tSeries,'detrend','Linear','spatialNormalization','Divide by mean','subtractMean', 'Yes', 'temporalNormalization', 'No');
  end
  
  % if there are any nans in the tSeries then don't fit
  if any(isnan(tSeries))
    if fitParams.verbose
      disp(sprintf('(pRFFit) Nan found in tSeries for voxel [%i %i %i] in scan %s:%i. Abandoning fit',x,y,z,viewGet(v,'groupName'),scanNum));
    end
    fit=[];return
  end
else
  tSeries = [];
end

% handle junk frames (i.e. ones that have not already been junked)
if ~isempty(fitParams.junkFrames) && ~isequal(fitParams.junkFrames,0)
  % drop junk frames
  disp(sprintf('(pRFFit) Dropping %i junk frames',fitParams.junkFrames));
  tSeries = tSeries(fitParams.junkFrames+1:end);
  if ~isfield(fitParams.concatInfo,'totalJunkedFramesIncludesJunked');
    fitParams.concatInfo.totalJunkedFrames = fitParams.concatInfo.totalJunkedFrames+fitParams.junkFrames;
    fitParams.concatInfo.totalJunkedFramesIncludesJunked = 1;
  end
end

% test to see if scan lengths and stim lengths match
tf = true;
for iScan = 1:fit.concatInfo.n
  sLength = fit.concatInfo.runTransition(iScan,2) - fit.concatInfo.runTransition(iScan,1) + 1;
  if sLength ~= size(fitParams.stim{iScan}.im,3)
    mrWarnDlg(sprintf('(pRFFit) Data length of %i for scan %i (concatNum:%i) does not match stimfile length %i',fit.concatInfo.runTransition(iScan,2),scanNum,iScan,size(fitParams.stim{iScan}.im,3)));
    %tf = false;
  end
end

if ~tf,fit = [];return,end

keyboard
global gPRFModels;
[params.pRFFit tSeries] = feval(gPRFModels(params.modelNum).initVoxel,params.pRFFit,tSeries);

% do prefit.
  % normalize tSeries to 0 mean unit length
  % get best r2 for all the models
%  [maxr bestModel] = max(r);
%  fitParams.initParams(1) = fitParams.prefit.x(bestModel);
%  fitParams.initParams(2) = fitParams.prefit.y(bestModel);
%  fitParams.initParams(3) = fitParams.prefit.rfHalfWidth(bestModel);
%  if fitParams.prefitOnly
    % return if we are just doing a prefit
%    fit = getFitParams(fitParams.initParams,fitParams);
%n    fit.rfType = fitParams.rfType;
%    fit.params = fitParams.initParams;
%    fit.r2 = maxr^2;
%    fit.r = maxr;
%    [fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);
%    fit.overlayParams = [fit.r2 fit.polarAngle fit.eccentricity fit.std];
    % display
%    if fitParams.verbose
%      disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f',fitParams.dispstr,x,y,z,fit.r2,r2d(fit.polarAngle),fit.eccentricity,fit.std));
%    end
%    return
%  end
%%%%%end


% set output arguments
fit = getFitParams(params,fitParams);
fit.rfType = fitParams.rfType;
fit.params = params;

% compute r^2
[residual modelResponse rfModel fit.r] = getModelResidual(params,tSeries,fitParams);
if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
  fit.r2 = 1-sum((residual-mean(residual)).^2)/sum((tSeries-mean(tSeries)).^2);
elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
  fit.r2 = residual^2;
end

% compute polar coordinates
[fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);
fit.overlayParams = [fit.r2 fit.polarAngle fit.eccentricity fit.std];

% display
if fitParams.verbose
  disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f',fitParams.dispstr,x,y,z,fit.r2,r2d(fit.polarAngle),fit.eccentricity,fit.std));
end

%%%%%%%%%%%%%%%%%%%%%%
%    setFitParams    %
%%%%%%%%%%%%%%%%%%%%%%
function fitParams = setFitParams(fitParams);

% set rfType
if ~isfield(fitParams,'rfType') || isempty(fitParams.rfType)
  fitParams.rfType = 'gaussian';
end

% get stimulus x,y and t
fitParams.stimX = fitParams.stim{1}.x;
fitParams.stimY = fitParams.stim{1}.y;
fitParams.stimT = fitParams.stim{1}.t;

% set stimulus extents
fitParams.stimExtents(1) = min(fitParams.stimX(:));
fitParams.stimExtents(3) = max(fitParams.stimX(:));
fitParams.stimExtents(2) = min(fitParams.stimY(:));
fitParams.stimExtents(4) = max(fitParams.stimY(:));
fitParams.stimWidth = fitParams.stimExtents(3)-fitParams.stimExtents(1);
fitParams.stimHeight = fitParams.stimExtents(4)-fitParams.stimExtents(2);

if ~isfield(fitParams,'initParams')
  % check the rfType to get the correct min/max arrays
  switch (fitParams.rfType)
   case 'gaussian'
    % parameter names/descriptions and other information for allowing user to set them
    fitParams.paramNames = {'x','y','rfWidth'};
    fitParams.paramDescriptions = {'RF x position','RF y position','RF width (std of gaussian)'};
    fitParams.paramIncDec = [1 1 1];
    fitParams.paramMin = [-inf -inf 0];
    fitParams.paramMax = [inf inf inf];
    % set min/max and init
    fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0];
    fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf];
    fitParams.initParams = [0 0 4];
   case 'gaussian-hdr'
    % parameter names/descriptions and other information for allowing user to set them
    fitParams.paramNames = {'x','y','rfWidth','timelag','tau'};
    fitParams.paramDescriptions = {'RF x position','RF y position','RF width (std of gaussian)','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
    fitParams.paramIncDec = [1 1 1 0.1 0.5];
    fitParams.paramMin = [-inf -inf 0 0 0];
    fitParams.paramMax = [inf inf inf inf inf];
    % set min/max and init
    fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0];
    fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf 3 inf];
    fitParams.initParams = [0 0 4 fitParams.timelag fitParams.tau];
    % add on parameters for difference of gamma
    if fitParams.diffOfGamma
      % parameter names/descriptions and other information for allowing user to set them
      fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
      fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
      fitParams.paramIncDec = [fitParams.paramsIncDec(:) 0.1 0.1 0.5];
      fitParams.paramMin = [fitParams.paramMin(:) 0 0 0];
      fitParams.paramMax = [fitParams.paramMax(:) inf inf inf];
      % set min/max and init
      fitParams.minParams = [fitParams.minParams 0 0 0];
      fitParams.maxParams = [fitParams.maxParams inf 6 inf];
      fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
    end
   otherwise
    disp(sprintf('(pRFFit:setFitParams) Unknown rfType %s',rfType));
    return
  end
  
  % round constraints
  fitParams.minParams = round(fitParams.minParams*10)/10;
  fitParams.maxParams = round(fitParams.maxParams*10)/10;

  % handle constraints here
  % Check if fit algorithm is one that allows constraints
  algorithmsWithConstraints = {'levenberg-marquardt'};
  if any(strcmp(fitParams.algorithm,algorithmsWithConstraints))
    % if constraints allowed then allow user to adjust them here (if they set defaultConstraints)
    if isfield(fitParams,'defaultConstraints') && ~fitParams.defaultConstraints
      % create a dialog to allow user to set constraints
      paramsInfo = {};
      for iParam = 1:length(fitParams.paramNames)
	paramsInfo{end+1} = {sprintf('min%s',fitParams.paramNames{iParam}) fitParams.minParams(iParam) sprintf('Minimum for parameter %s (%s)',fitParams.paramNames{iParam},fitParams.paramDescriptions{iParam}) sprintf('incdec=[%f %f]',-fitParams.paramIncDec(iParam),fitParams.paramIncDec(iParam)) sprintf('minmax=[%f %f]',fitParams.paramMin(iParam),fitParams.paramMax(iParam))};
	paramsInfo{end+1} = {sprintf('max%s',fitParams.paramNames{iParam}) fitParams.maxParams(iParam) sprintf('Maximum for parameter %s (%s)',fitParams.paramNames{iParam},fitParams.paramDescriptions{iParam})  sprintf('incdec=[%f %f]',-fitParams.paramIncDec(iParam),fitParams.paramIncDec(iParam)) sprintf('minmax=[%f %f]',fitParams.paramMin(iParam),fitParams.paramMax(iParam))};
      end
      params = mrParamsDialog(paramsInfo,'Set parameter constraints');
      % if params is not empty then set them
      if isempty(params)
	disp(sprintf('(pRFFit) Using default constraints'));
      else
	% get the parameter constraints back from the dialog entries
	for iParam = 1:length(fitParams.paramNames)
	  fitParams.minParams(iParam) = params.(sprintf('min%s',fitParams.paramNames{iParam}));
	  fitParams.maxParams(iParam) = params.(sprintf('max%s',fitParams.paramNames{iParam}));
	end
      end
`    end
    % Now display parameter constraints
    for iParam = 1:length(fitParams.paramNames)
      disp(sprintf('(pRFFit) Parameter %s [min:%f max:%f] (%i:%s)',fitParams.paramNames{iParam},fitParams.minParams(iParam),fitParams.maxParams(iParam),iParam,fitParams.paramDescriptions{iParam}));
    end
  else
    % no constraints allowed
    disp(sprintf('(pRFFit) !!! Fit constraints ignored for algorithm: %s (if you want to constrain the fits, then use: %s) !!!',fitParams.algorithm,cell2mat(algorithmsWithConstraints)));
  end
end

fitParams.nParams = length(fitParams.initParams);

% optimization parameters
if ~isfield(fitParams,'algorithm') || isempty(fitParams.algorithm)
  fitParams.algorithm = 'nelder-mead';
end
fitParams.optimParams = optimset('MaxIter',inf,'Display',fitParams.optimDisplay);

% compute number of frames
fitParams.nFrames = size(fitParams.stim{1}.im,3);

% parameters for converting the stimulus
params = {'xFlip','yFlip','timeShiftStimulus'};
for i = 1:length(params)
  if ~isfield(fitParams,params{i}) || isempty(fitParams.(params{i}))
    fitParams.(params{i}) = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getModelResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual modelResponse rfModel r] = getModelResidual(params,tSeries,fitParams,justGetModel)

residual = [];
if nargin < 4, justGetModel = 0;end



%%%%%%%%%%%%%%%%%%%%%%
%    dispModelFit    %
%%%%%%%%%%%%%%%%%%%%%%
function dispModelFit(params,fitParams,modelResponse,tSeries,rfModel)

mlrSmartfig('pRFFit_getModelResidual','reuse');
clf
subplot(4,4,[1:3 5:7 9:11 13:15]);
%plot(fitParams.stimT(fitParams.junkFrames+1:end),tSeries,'k-');
plot(tSeries,'k-');
hold on
%plot(fitParams.stimT(fitParams.junkFrames+1:end),modelResponse,'r-');
plot(modelResponse,'r-');
xlabel('Time (sec)');
ylabel('BOLD (% sig change)');
p = getFitParams(params,fitParams);
titleStr = sprintf('x: %s y: %s rfHalfWidth: %s',mlrnum2str(p.x),mlrnum2str(p.y),mlrnum2str(p.std));
titleStr = sprintf('%s\n(timelag: %s tau: %s exponent: %s)',titleStr,mlrnum2str(p.canonical.timelag),mlrnum2str(p.canonical.tau),mlrnum2str(p.canonical.exponent));
if p.canonical.diffOfGamma
  titleStr = sprintf('%s - %s x (timelag2: %s tau2: %s exponent2: %s)',titleStr,mlrnum2str(p.canonical.amplitudeRatio),mlrnum2str(p.canonical.timelag2),mlrnum2str(p.canonical.tau2),mlrnum2str(p.canonical.exponent2));
end
title(titleStr);
axis tight

subplot(4,4,[8 12 16]);
imagesc(fitParams.stimX(:,1),fitParams.stimY(1,:),flipud(rfModel'));
axis image;
hold on
hline(0);vline(0);

subplot(4,4,4);cla
p = getFitParams(params,fitParams);
canonical = getCanonicalHRF(p.canonical,fitParams.framePeriod);
plot(canonical.time,canonical.hrf,'k-')
if exist('myaxis') == 2,myaxis;end

%%%%%%%%%%%%%%%%%%%%%%%%
%    scaleAndOffset    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [modelResponse residual] = scaleAndOffset(modelResponse,tSeries)

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

%%%%%%%%%%%%%%%%%%%
%    parseArgs    %
%%%%%%%%%%%%%%%%%%%
function [v scanNum x y s fitParams tSeries] = parseArgs(args);

v = [];scanNum=[];x=[];y=[];s=[];fitParams=[];tSeries = [];

% check for calling convention from interrogator
if (length(args) >= 7) && isnumeric(args{6})
  v = args{1};
  %overlayNum = args{2};
  scanNum = args{3};
  x = args{4};
  y = args{5};
  s = args{6};
  %roi = args{7};
  fitParams.dispFit = true;
  fitParams.optimDisplay = 'final';
  fitParams.algorithm = 'nelder-mead';
  fitParams.getModelResponse = false;
  fitParams.prefit = [];
  fitParams.xFlipStimulus = 0;
  fitParams.yFlipStimulus = 0;
  fitParams.timeShiftStimulus = 0;
  fitParams.betaEachScan = false;
  fitParams.justGetStimImage = false;
  fitParams.returnPrefit = false;
  fitParams.verbose = 1;
  fitParams.timelag = 1;
  fitParams.tau = 0.6;
  fitParams.exponent = 6;
  clearConstraints = false;
  getArgs({args{8:end}},{'fitTypeParams=[]'});
  if isempty(fitTypeParams)
    % no fit type params, check if we have them set in
    % the global (this is useful so that when called as an 
    % interrogator we don't have to keep setting them
    global gpRFFitTypeParams
    % if user is holding shift, then reget parameters
    if ~isempty(gcf) && any(strcmp(get(gcf,'CurrentModifier'),'shift'))
      gpRFFitTypeParams = [];
    end
    % get the parameters from the user interface if not already set
    if isempty(gpRFFitTypeParams)
      fitTypeParams = pRFGUI('pRFFitParamsOnly=1','v',v);
      if isempty(fitTypeParams) 
	v = [];
	return
      end
      gpRFFitTypeParams = fitTypeParams;
      % flag to clear the constraints
      clearConstraints = true;
    else
      % otherwise grab old ones
      disp(sprintf('(pRFFit) Using already set parameters to compute pRFFit. If you want to use different parameters, hold shift down as you click the next voxel'));
      fitTypeParams = gpRFFitTypeParams;
    end
  end
  if ~isempty(fitTypeParams)
    % if fitTypeParams is passed in (usually from pRF / pRFGUI) then
    % grab parameters off that structure
    fitTypeParamsFields = fieldnames(fitTypeParams);
    for i = 1:length(fitTypeParamsFields)
      fitParams.(fitTypeParamsFields{i}) = fitTypeParams.(fitTypeParamsFields{i});
    end
  end
  
% normal calling convention
elseif length(args) >= 5
  v = args{1};
  scanNum = args{2};
  x = args{3};
  y = args{4};
  s = args{5};
  % parse anymore argumnets
  dispFit=[];stim = [];getModelResponse = [];params = [];concatInfo = [];prefit = [];
  xFlip=[];yFlip=[];timeShiftStimulus=[];rfType=[];betaEachScan=[];fitTypeParams = [];
  dispIndex = [];dispN = [];returnPrefit = [];tSeries=[];quickPrefit=[];junkFrames=[];
  verbose = [];justGetStimImage = [];framePeriod = [];
  getArgs({args{6:end}},{'dispFit=0','stim=[]','getModelResponse=0','params=[]','concatInfo=[]','prefit=[]','xFlipStimulus=0','yFlipStimulus=0','timeShiftStimulus=0','rfType=gaussian','betaEachScan=0','fitTypeParams=[]','justGetStimImage=[]','verbose=1','dispIndex=[]','dispN=[]','returnPrefit=0','quickPrefit=0','tSeries=[]','junkFrames=[]','framePeriod=[]','paramsInfo=[]'});
  % default to display fit
  fitParams.dispFit = dispFit;
  fitParams.stim = stim;
  fitParams.optimDisplay = 'off';
  fitParams.getModelResponse = getModelResponse;
  fitParams.params = params;
  fitParams.concatInfo = concatInfo;
  fitParams.prefit = prefit;
  fitParams.xFlipStimulus = xFlipStimulus;
  fitParams.yFlipStimulus = yFlipStimulus;
  fitParams.timeShiftStimulus = timeShiftStimulus;
  fitParams.rfType = rfType;
  fitParams.betaEachScan = betaEachScan;
  fitParams.justGetStimImage = justGetStimImage;
  fitParams.verbose = verbose;
  fitParams.returnPrefit = returnPrefit;
  fitParams.junkFrames = junkFrames;
  fitParams.framePeriod = framePeriod;
  % now read in all the fields in the paramsInfo
  if ~isempty(paramsInfo)
    paramsInfoFields = fieldnames(paramsInfo);
    for iField = 1:length(paramsInfoFields)
      fitParams.(paramsInfoFields{iField}) = paramsInfo.(paramsInfoFields{iField});
    end
  end
  if ~isempty(fitTypeParams)
    % if fitTypeParams is passed in (usually from pRF / pRFGUI) then
    % grab parameters off that structure
    fitTypeParamsFields = fieldnames(fitTypeParams);
    for i = 1:length(fitTypeParamsFields)
      fitParams.(fitTypeParamsFields{i}) = fitTypeParams.(fitTypeParamsFields{i});
    end
  end
  if ~isempty(dispIndex) && ~isempty(dispN)
    % create a display string. Note that we use sprintf twice here so that
    % we can create a string with the proper amount of space padding the index
    % so that each row always displays as the same length string
    prefitOnlyStr = '';
    if isfield(fitParams,'prefitOnly') && fitParams.prefitOnly
      prefitOnlyStr = ' (prefit only)';
    end
    fitParams.dispstr = sprintf(sprintf('Voxel %%%i.f/%%i%%s: ',length(sprintf('%i',dispN))),dispIndex,dispN,prefitOnlyStr);
    end
  if getModelResponse && isempty(params)
    disp(sprintf('(pRFFit) Must pass in params when using getModelResponse'));
    fitParams.getModelResponse = false;
  end
else
  help pRFFit;
end

% some default parameters
if ~isfield(fitParams,'prefitOnly') || isempty(fitParams.prefitOnly)
  fitParams.prefitOnly = false;
end
if ~isfield(fitParams,'dispstr')
  fitParams.dispstr = '';
end
if ~isfield(fitParams,'quickPrefit') || isempty(fitParams.quickPrefit)
  fitParams.quickPrefit = false;
end
if ~isfield(fitParams,'verbose') || isempty(fitParams.verbose)
  fitParams.verbose = true;
end

% get some info about the scanNum
if ~isfield(fitParams,'framePeriod') || isempty(fitParams.framePeriod)
  fitParams.framePeriod = viewGet(v,'framePeriod');
end
if ~isfield(fitParams,'junkFrames') || isempty(fitParams.junkFrames)
  fitParams.junkFrames = viewGet(v,'junkFrames',scanNum);
end



%%%%%%%%%%%%%
%%   r2d   %%
%%%%%%%%%%%%%
function degrees = r2d(angle)

degrees = (angle/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = convolveModelWithStimulus(rfModel,stim,nFrames)

% get number of frames
nStimFrames = size(stim.im,3);

% preallocate memory
modelResponse = zeros(1,nFrames);

for frameNum = 1:nStimFrames
  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
  modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
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

