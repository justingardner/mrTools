% pR_somatoFFit
%
%      usage: pRF_somatoFit(v,scanNum,x,y,s,<dispFit=0>)
%         by: ds / completely base on code by justin gardner
%       date: 201602 [orig 11/14/2011]
%    purpose: interrogator that fits pRF model to selected voxel
%
function fit = pRF_somatoFit(varargin)

fit = [];
% parse input arguments - note that this is set
% up so that it can also be called as an interrogator

[v, scanNum, x, y, z, fitParams, tSeries, hrfprf] = parseArgs(varargin);
%[v, scanNum, x, y, z, fitParams, tSeries] = parseArgs(varargin);
if isempty(v),return,end

% get concat info
if ~isfield(fitParams,'concatInfo') || isempty(fitParams.concatInfo)
    fitParams.concatInfo = viewGet(v,'concatInfo',scanNum);
end

if ~ieNotDefined('hrfprf')
    hrfprfcheck = 1;
else
    hrfprfcheck = 0;
end

% if there is no concatInfo, then make one that will
% treat the scan as a single scan
if isempty(fitParams.concatInfo)
    nFrames = viewGet(v,'nFrames',scanNum);
    fitParams.concatInfo.isConcat = false;
    fitParams.concatInfo.n = 1;
    fitParams.concatInfo.whichScan = ones(1,nFrames);
    fitParams.concatInfo.whichVolume = 1:nFrames;
    fitParams.concatInfo.runTransition = [1 nFrames];
    fitParams.concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
    if length(fitParams.concatInfo.totalJunkedFrames > 1)
        % first check for consistency in totalJunkedFrames
        if length(unique(fitParams.concatInfo.totalJunkedFrames)) > 1
            disp(sprintf('(pRFFit) totalJunkedFrames are different for different members of component scans - could be an average in which different scans with different number of junked frames were removed. This could cause a problem in computing what the stimulus was for the average. The total junked frames count was: %s, but we will use %i as the actual value for computing the stimulus',num2str(fitParams.concatInfo.totalJunkedFrames),floor(median(fitParams.concatInfo.totalJunkedFrames))));
        end
        fitParams.concatInfo.totalJunkedFrames = floor(median(fitParams.concatInfo.totalJunkedFrames));
    end
else
    fitParams.concatInfo.isConcat = true;
    if ~isfield(fitParams.concatInfo,'totalJunkedFrames')
        fitParams.concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
    end
end

% get the stimulus movie if it wasn't passed in
if ~isfield(fitParams,'stim') || isempty(fitParams.stim)
    fitParams.stim = getStim(v,scanNum,fitParams);
end
if isempty(fitParams.stim),return,end

% if we are being called to just return the stim image
% then return it here
if fitParams.justGetStimImage
    fit = fitParams.stim;
    return
end

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
            disp(sprintf('(pRF_somatoFit) Nan found in tSeries for voxel [%i %i %i] in scan %s:%i. Abandoning fit',x,y,z,viewGet(v,'groupName'),scanNum));
        end
        fit=[];return
    end
else
    tSeries = [];
end

% handle junk frames (i.e. ones that have not already been junked)
if ~isempty(fitParams.junkFrames) && ~isequal(fitParams.junkFrames,0)
    % drop junk frames
    disp(sprintf('(pRF_somatoFit) Dropping %i junk frames',fitParams.junkFrames));
    tSeries = tSeries(fitParams.junkFrames+1:end);
    if ~isfield(fitParams.concatInfo,'totalJunkedFramesIncludesJunked');
        fitParams.concatInfo.totalJunkedFrames = fitParams.concatInfo.totalJunkedFrames+fitParams.junkFrames;
        fitParams.concatInfo.totalJunkedFramesIncludesJunked = 1;
    end
end


% set up the fit routine params
fitParams = setFitParams(fitParams);

% just return model response for already calculated params
if fitParams.getModelResponse
    
    
    if ~ieNotDefined('hrfprf')
        fitParams.hrfprf = hrfprf;
        % get model fit
        [residual, fit.modelResponse, fit.rfModel, ~, realhrf] = getModelResidual(fitParams.params,tSeries,fitParams, [], hrfprfcheck);
        % get the real hrf from the inputted ones
        fit.p = getFitParams(fitParams.params,fitParams);
        fit.canonical = realhrf;
        %fit.canonical = getCanonicalHRF(fit.p.canonical,fitParams.framePeriod);
        % return tSeries
        fit.tSeries = tSeries;
        return;
    else
        
        % get model fit
        [residual fit.modelResponse fit.rfModel] = getModelResidual(fitParams.params,tSeries,fitParams, [], hrfprfcheck);
        % get the canonical
        fit.p = getFitParams(fitParams.params,fitParams);
        fit.canonical = getCanonicalHRF(fit.p.canonical,fitParams.framePeriod);
        % return tSeries
        fit.tSeries = tSeries;
        return;
        
    end
end

% return some fields
fit.stim = fitParams.stim;
fit.stimX = fitParams.stimX;
fit.stimY = fitParams.stimY;
fit.stimT = fitParams.stimT;
fit.concatInfo = fitParams.concatInfo;
fit.nParams = fitParams.nParams;
paramsInfoFields = {'minParams','maxParams','initParams','paramNames','paramDescriptions'};
for iField = 1:length(paramsInfoFields)
    fit.paramsInfo.(paramsInfoFields{iField}) = fitParams.(paramsInfoFields{iField});
end

% test to see if scan lengths and stim lengths match
tf = true;
for iScan = 1:fit.concatInfo.n
    sLength = fit.concatInfo.runTransition(iScan,2) - fit.concatInfo.runTransition(iScan,1) + 1;
    if sLength ~= size(fitParams.stim{iScan}.im,3)
        mrWarnDlg(sprintf('(pRF_somatoFit) Data length of %i for scan %i (concatNum:%i) does not match stimfile length %i',fit.concatInfo.runTransition(iScan,2),scanNum,iScan,size(fitParams.stim{iScan}.im,3)));
        tf = false;
    end
end

if ~tf,fit = [];return,end

% do prefit. This computes (or is passed in precomputed) model responses
% for a variety of parameters and calculates the correlation between
% the models and the time series. The one that has the best correlation
% is then used as the initial parameters for the nonlinear fit. This
% helps prevent getting stuck in local minima
if isfield(fitParams,'prefit') && ~isempty(fitParams.prefit)
    params = fitParams.initParams;
    % calculate model if not already calculated
    if ~isfield(fitParams.prefit,'modelResponse')
        % get number of workers
        nProcessors = mlrNumWorkers;
        mlrDispPercent(-inf,sprintf('(pRF_somatoFit) Computing %i prefit model responses using %i processors',fitParams.prefit.n,nProcessors));
        % first convert the x/y and width parameters into sizes
        % on the actual screen
        %fitParams.prefit.x = fitParams.prefit.x;% *fitParams.stimWidth;
        %fitParams.prefit.y = fitParams.prefit.y; %*fitParams.stimHeight;
        %fitParams.prefit.rfHalfWidth = fitParams.prefit.rfHalfWidth; % *max(fitParams.stimWidth,fitParams.stimHeight);
        %fitParams.prefit.x = fitParams.prefit.x *fitParams.stimWidth;
        %fitParams.prefit.y = fitParams.prefit.y *fitParams.stimHeight;
        %fitParams.prefit.rfHalfWidth = fitParams.prefit.rfHalfWidth *max(fitParams.stimWidth,fitParams.stimHeight);
        %fitParams.prefit.hrfDelay = fitParams.prefit.hrfDelay; %*max(fitParams.stimWidth,fitParams.stimHeight);
        % init modelResponse
        allModelResponse = nan(fitParams.prefit.n,fitParams.concatInfo.runTransition(end,end));
        % compute all the model response, using parfor loop
        % parfor i = 1:fitParams.prefit.n
        
        parfor i = 1:fitParams.prefit.n
            % fit the model with these parameters
            %[residual modelResponse rfModel] = getModelResidual([fitParams.prefit.x(i) fitParams.prefit.y(i) fitParams.prefit.rfHalfWidth(i) fitParams.prefit.hrfDelay(i) params(4:end)],tSeries,fitParams,1);
            %[residual modelResponse rfModel] = getModelResidual([fitParams.prefit.x(i) fitParams.prefit.y(i) fitParams.prefit.rfHalfWidth(i) params(4:end)],tSeries,fitParams,1, crossValcheck);
            [residual modelResponse rfModel] = getModelResidual([fitParams.prefit.x(i) fitParams.prefit.y(i) fitParams.prefit.rfHalfWidth(i) params(4:end)],tSeries,fitParams,1);
            % normalize to 0 mean unit length
            allModelResponse(i,:) = (modelResponse-mean(modelResponse))./sqrt(sum(modelResponse.^2))';
            if fitParams.verbose
                disp(sprintf('(pRF_somatoFit) Computing prefit model response %i/%i: Center [%6.2f,%6.2f] rfHalfWidth=%5.2f ',i,fitParams.prefit.n,fitParams.prefit.x(i),fitParams.prefit.y(i),fitParams.prefit.rfHalfWidth(i) ));
            end
        end
        mlrDispPercent(inf);
        fitParams.prefit.modelResponse = allModelResponse;
        clear allModelResponse;
    end
    % save in global, so that when called as an interrogator
    % we don't have to keep computing fitParams
    global gpRFFitTypeParams
    gpRFFitTypeParams.prefit = fitParams.prefit;
    % return some computed fields
    fit.prefit = fitParams.prefit;
    if fitParams.returnPrefit,return,end
    % normalize tSeries to 0 mean unit length
    tSeriesNorm = (tSeries-mean(tSeries))/sqrt(sum(tSeries.^2));
    % calculate r for all modelResponse by taking inner product
    r = fitParams.prefit.modelResponse*tSeriesNorm;
    % get best r2 for all the models
    [maxr bestModel] = max(r);
    fitParams.initParams(1) = fitParams.prefit.x(bestModel);
    fitParams.initParams(2) = fitParams.prefit.y(bestModel);
    fitParams.initParams(3) = fitParams.prefit.rfHalfWidth(bestModel);
    %fitParams.initParams(4) = fitParams.prefit.hrfDelay(bestModel);
    if fitParams.prefitOnly
        % return if we are just doing a prefit
        fit = getFitParams(fitParams.initParams,fitParams);
        fit.rfType = fitParams.rfType;
        fit.params = fitParams.initParams;
        fit.r2 = maxr^2;
        fit.r = maxr;
        % [fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);
        fit.prefDigit = fit.y; % Is this flipped?
        fit.prefPD = fit.x;
        fit.rfHalfWidth = fit.std;
        %fit.hrfDelay = fit.params(4);
        
        
        %% anon function see pRFFit.m
        %         mod = 'somato'; % this variable should be set in the GUI - the user can choose the stimulus / modality
        %         overlayNames = getMetaData(v,params,mod,'overlayNames');
        %         % r2
        %         eval(sprintf('fit.%s = fit.r2',overlayNames{1}));
        %         % x
        %         eval(sprintf('fit.%s = fit.x',overlayNames{2}));
        %
        %         if numel(overlayNames) == 4
        %             % y
        %             eval(sprintf('fit.%s = fit.y',overlayNames{3}));
        %             % hw
        %             eval(sprintf('fit.%s = fit.std',overlayNames{4}));
        %         else
        %             % hw
        %             eval(sprintf('fit.%s = fit.std',overlayNames{3}));
        %         end
        %%%%%%
        
        % display
        if fitParams.verbose
            % disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f',fitParams.dispstr,x,y,z,fit.r2,r2d(fit.polarAngle),fit.eccentricity,fit.std));
            disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f prefDigit=%6.1f prefPD=%6.1f rfHalfWidth=%6.1f ',fitParams.dispstr,x,y,z,fit.r2,fit.prefDigit,fit.prefPD,fit.std ));
            
        end
        return
    end
end

% this works if we want to run the pRF on already computed values, e.g.
% cross validation. But in the case where we want to fit a few parameters,
% e.g. fit for the first 3 gaussian parameters, but give it some
% precomputed HRF params, then we need to be selective in our params
%
% An easy way would be to let the fit happen as normal, but then overwrite
% the params(end-5:end) with those precomputed HRF params - but this is
% very dirty and slow...
%
%
if ~ieNotDefined('hrfprf')
    fitParams.hrfprf = hrfprf;
    if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
        [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@getModelResidual,fitParams.initParams,fitParams.minParams,fitParams.maxParams,fitParams.optimParams,tSeries,fitParams);
    elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
        [params fval exitflag] = fminsearch(@getModelResidual,fitParams.initParams,fitParams.optimParams,(tSeries-mean(tSeries))/var(tSeries.^2),fitParams);
        %[params fval exitflag] = fmincon(@getModelResidual,fitParams.initParams,[],[],[],[],[-5 -5 -5],[5 5 5],[],fitParams.optimParams,(tSeries-mean(tSeries))/var(tSeries.^2),fitParams);
    else
        disp(sprintf('(pRF_somatoFit) Unknown optimization algorithm: %s',fitParams.algorithm));
        return
    end
    
else
    
    % now do nonlinear fit
    if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
        [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@getModelResidual,fitParams.initParams,fitParams.minParams,fitParams.maxParams,fitParams.optimParams,tSeries,fitParams);
    elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
        [params fval exitflag] = fminsearch(@getModelResidual,fitParams.initParams,fitParams.optimParams,(tSeries-mean(tSeries))/var(tSeries.^2),fitParams);
    else
        disp(sprintf('(pRF_somatoFit) Unknown optimization algorithm: %s',fitParams.algorithm));
        return
    end
    
end

% set output arguments
fit = getFitParams(params,fitParams);
fit.rfType = fitParams.rfType;
fit.params = params;

% compute r^2
[residual modelResponse rfModel fit.r] = getModelResidual(params,tSeries,fitParams, [], hrfprfcheck);
%fit.r = r;
if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
    fit.r2 = 1-sum((residual-mean(residual)).^2)/sum((tSeries-mean(tSeries)).^2);
elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
    fit.r2 = residual^2;
end

fit.modelResponse = modelResponse;
fit.tSeries = tSeries;
if strcmp(lower(fitParams.algorithm),'levenberg-marquardt') % this may crash if running Nelder-Mead
    fit.residual = residual;
end

% compute polar coordinates
% [fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);
switch fit.rfType
    case {'gaussian-1D'} % WITHIN DIGIT MODEL


        % hang on, this is stupid
        % this is ignoring the fact that we made the rf model Gaussian in
        % the first place!!!

        %         if size(fitParams.stimX,1) == 4
        %             [~, index] = max([fit.amp1 fit.amp2 fit.amp3 fit.amp4]);
        %             fit.prefDigit = index;
        %             allMeans = [fit.meanOne fit.meanTwo fit.meanThr fit.meanFour];
        %             fit.prefPD = allMeans(index);% Take the prefDigit and specify the
        %             allStd = [fit.stdOne fit.stdTwo fit.stdThr fit.stdFour];
        %             fit.rfHalfWidth = allStd(index);
        %         else
        %             [~, index] = max([fit.amp1 fit.amp2 fit.amp3]);
        %             fit.prefDigit = index;
        %             allMeans = [fit.meanOne fit.meanTwo fit.meanThr];
        %             fit.prefPD = allMeans(index);% Take the prefDigit and specify the
        %             allStd = [fit.stdOne fit.stdTwo fit.stdThr ];
        %             fit.rfHalfWidth = allStd(index);
        %         end
        %keyboard
        %[ma] 2019

        if size(fitParams.stimX,2) == 5 % for TW Touchmap
            [~,index] = max(rfModel);
        else
            [~,index] = max(max(rfModel)); % max digit amp of Gaussian
        end
        % ds thinks.. this is the correct way to read off the "other" direction
        % reason the above is not quite correct is that at this point the RF model is flipped.
        % so other way of changing this would be to apply rules to transpose(rfModel) ?!
        %[~,index] = max(max(rfModel,[],2),[], 1); % max digit amp of Gaussian

        fit.prefDigit = index;
        %figure, plot(rfModel)
        %disp(fit.prefDigit)
        %fit.prefPD = index;

        %          %keyboard
        %          % this finds prefPD from X, not y axis!
        %          % problem: r2 for 1D Gaussian is wrong, because it's using a
        %          % non-interp Gaussian model for the fit, here we do it too late,
        %          % just for outputs.
        %          % You'd have to change the stimfile to have 100 x points.
        thisX = linspace(1,size(rfModel,2),100);

        if size(fitParams.stimX,2) == 5
            Pone = [fit.amp1 fit.meanOne fit.stdOne 0];
            thisStd = fit.stdOne;
        else


            if fit.prefDigit == 1
                Pone = [fit.amp1 fit.meanOne fit.stdOne 0];
                thisStd = fit.stdOne;
            elseif fit.prefDigit == 2
                Pone = [fit.amp2 fit.meanTwo fit.stdTwo 0];
                thisStd = fit.stdTwo;
            elseif fit.prefDigit == 3
                Pone = [fit.amp3 fit.meanThr fit.stdThr 0];
                thisStd = fit.stdThr;
            elseif fit.prefDigit == 4
                Pone = [fit.amp4 fit.meanFour fit.stdFour 0];
                thisStd = fit.stdFour;
            elseif fit.prefDigit == 5
                Pone = [fit.amp5 fit.meanFive fit.stdFive 0];
                thisStd = fit.stdFive;
            end

        end

        %          if fit.prefPD == 1
        %              Pone = [fit.amp1 fit.meanOne fit.stdOne 0];
        %              thisStd = fit.stdOne;
        %          elseif fit.prefPD == 2
        %              Pone = [fit.amp2 fit.meanTwo fit.stdTwo 0];
        %              thisStd = fit.stdTwo;
        %          elseif fit.prefPD == 3
        %              Pone = [fit.amp3 fit.meanThr fit.stdThr 0];
        %              thisStd = fit.stdThr;
        %          elseif fit.prefPD == 4
        %              Pone = [fit.amp4 fit.meanFour fit.stdFour 0];
        %              thisStd = fit.stdFour;
        %          end

        thisZ = gauss(Pone,thisX)';
        %thisZ = gauss(Pone,thisX);
        [~,index2] = max(thisZ);
        fit.prefPD = thisX(index2); % this is trickier, because we need to rewrap the PD values onto a 1-4 grid
        %fit.prefDigit = thisX(index2);
        %%%%%fit.prefPD = mean(rfModel(:,index)); % mean of that digit = PD
        fit.prefPD = abs(fit.prefPD - (size(rfModel,2)+1) ); % I'm not sure why this is necessary, but it needs to be unflipped outside of mrTools
        %fit.prefDigit = abs(fit.prefDigit - (size(rfModel,2)+1));



        %fit.rfHalfWidth = std(rfModel(:,index)); % std of digit Gaussian
        fit.rfHalfWidth = thisStd;

    case {'gaussian-1D-transpose'} % BETWEEN DIGIT MODEL


        if size(fitParams.stimX,2) == 5 % for TW Touchmap
            [~,index] = max(rfModel);
        else
            [~,index] = max(max(rfModel,[],2),[], 1); % max digit amp of Gaussian
        end


        %[~,index] = max(max(rfModel,[],2),[], 1); % max digit amp of Gaussian
        fit.prefPD = index;

        thisX = linspace(1,size(rfModel,2),100);

        if size(fitParams.stimX,2) == 5
            Pone = [fit.amp1 fit.meanOne fit.stdOne 0];
            thisStd = fit.stdOne;
        else


            if fit.prefPD == 1
                Pone = [fit.amp1 fit.meanOne fit.stdOne 0];
                thisStd = fit.stdOne;
            elseif fit.prefPD == 2
                Pone = [fit.amp2 fit.meanTwo fit.stdTwo 0];
                thisStd = fit.stdTwo;
            elseif fit.prefPD == 3
                Pone = [fit.amp3 fit.meanThr fit.stdThr 0];
                thisStd = fit.stdThr;
            elseif fit.prefPD == 4
                Pone = [fit.amp4 fit.meanFour fit.stdFour 0];
                thisStd = fit.stdFour;
            elseif fit.prefPD == 5
                Pone = [fit.amp5 fit.meanFive fit.stdFive 0];
                thisStd = fit.stdFive;

            end

        end
        thisZ = gauss(Pone,thisX)';
        [~,index2] = max(thisZ);

        % weird buglet here wrt colormaps for base/tips pRF - ma jan2020
        fit.prefDigit = thisX(index2); % original
        %prefDigit_tmp = thisX(index2);
        %fit.prefDigit = abs(prefDigit_tmp - (size(rfModel,2)+1) );


        % ignore this for base/tips pRF fitting

        % this is trickier, because we need to rewrap the PD values onto a 1-4 grid
        if size(fitParams.stimX,2) == 5
            fit.prefPD = fit.prefPD;
        else
            fit.prefPD = abs(fit.prefPD - (size(rfModel,2)+1) );
        end
        fit.rfHalfWidth = thisStd;



    case {'gaussian','gaussian-hdr'}
        fit.prefDigit = fit.y;
        fit.prefPD = fit.x;

        fit.rfHalfWidth = fit.std;
    case {'gaussian-hdr-double'}

        % test
        %         stimsY = meshgrid(1:0.03:4,1:0.03:4);
        %         stimsX = transpose(stimsY);
        %         bloop = exp(-(((stimsX-fit.x).^2)/(2*(fit.stdx^2))+((stimsY-fit.y).^2)/(2*(fit.stdy^2))));
        %         figure, imagesc(bloop)
        if size(fitParams.stimX,2) == 5
            fit.prefDigit = fit.x;
            fit.prefPD = fit.y;
        else


            fit.prefDigit = fit.y;
            fit.prefPD = fit.x; % for some reason this is flipped top to bottom.
        end

        fit.rfHalfWidthX = fit.stdx;
        fit.rfHalfWidthY = fit.stdy;

    case {'gaussian-1D-orthotips'}


        fit.prefDigit = mean(rfModel);
        fit.prefPD = fit.y;
        fit.rfHalfWidth = std(rfModel);

        %     case {'gaussian-tips', 'gaussian-tips-hdr'}
        %         fit.prefDigit = fit.meanOne;
        %         fit.prefPD = fit.y;
        %         fit.rfHalfWidth = fit.std;
        %
    case {'sixteen-hdr','nine-param','nine-param-hdr', 'five-hdr'}
        %figure, imagesc(rfModel)
        [~,fit.prefDigit] = max(mean(rfModel,1)); % in 1st dim
        [~,fit.prefPD]    = max(mean(flipud(rfModel),2)); %orthog dim
        %[~,fit.prefPD]    = max(mean(rfModel,2)); %orthog dim
        fit.rfHalfWidth = std(rfModel(:));
    otherwise
        fit.prefDigit = fit.x; %Flipped this, otherwise prefPD and prefDigit are incorrectly labelled in GUI, i.e. fit.prefDigit = fit.y
        fit.prefPD = fit.y;
        fit.rfHalfWidth = fit.std;

end


% display
if fitParams.verbose
    
    if strcmpi(fit.rfType,'gaussian-hdr-double')
        disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f prefDigit=%6.1f prefPD=%6.1f STDX=%6.1f STDY=%6.1f ',fitParams.dispstr,x,y,z,fit.r2,fit.prefDigit,fit.prefPD,fit.rfHalfWidthX, fit.rfHalfWidthY ));

    else
        
        % disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f polarAngle=%6.1f eccentricity=%6.1f rfHalfWidth=%6.1f',fitParams.dispstr,x,y,z,fit.r2,r2d(fit.polarAngle),fit.eccentricity,fit.std));
        disp(sprintf('%s[%2.f %2.f %2.f] r2=%0.2f prefDigit=%6.1f prefPD=%6.1f rfHalfWidth=%6.1f ',fitParams.dispstr,x,y,z,fit.r2,fit.prefDigit,fit.prefPD,fit.rfHalfWidth ));
    end
    
    
end

end


%%%%%%%%%%%%%%%%%%%%%%
%%    setFitParams    %
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
% need to change this to 4 for more params??
%fitParams.stimWidth = fitParams.stimExtents(3)-fitParams.stimExtents(1);
%fitParams.stimHeight = fitParams.stimExtents(4)-fitParams.stimExtents(2);

if ~isfield(fitParams,'initParams')
    % check the rfType to get the correct min/max arrays
    switch (fitParams.rfType)
        case 'gaussian'
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'x','y','rfWidth'};
            fitParams.paramDescriptions = {'RF x position (digit)','RF y position (PD)','RF width (std of gaussian)'};
            fitParams.paramIncDec = [1 1 1 ];
            fitParams.paramMin = [0 0 0  ];
            fitParams.paramMax = [inf inf 10 ];
            % set min/max and init
            fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 ];
            fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf ];
            fitParams.initParams = [0 0 4 ];
        case 'gaussian-hdr'
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'x','y','rfWidth','timelag','tau'};
            fitParams.paramDescriptions = {'RF x position (digit)','RF y position (PD)','RF width (std of gaussian)','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1 0.1 0.5];
            fitParams.paramMin = [-inf -inf 0 0 0];
            fitParams.paramMax = [inf inf inf 10 10];
            % set min/max and init
            fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0];
            fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf 10 10];
            fitParams.initParams = [0 0 4 fitParams.timelag fitParams.tau];
            % add on parameters for difference of gamma
            if fitParams.diffOfGamma
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' 20 20 20 ];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams 20 20 20 ];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
        case 'gaussian-hdr-double'
            fitParams.paramNames = {'x','y','stdx','stdy','timelag','tau'};
            fitParams.paramDescriptions = {'RF x position (digit)','RF y position (PD)','stdx','stdy','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1 1 0.1 0.5];
            fitParams.paramMin = [-inf -inf 0 0 0 0];
            fitParams.paramMax = [inf inf inf inf 10 10];
            % set min/max and init
            fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0 0];
            fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf inf 10 10];
            fitParams.initParams = [0 0 1 1 fitParams.timelag fitParams.tau];
            % add on parameters for difference of gamma
            if fitParams.diffOfGamma
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' 20 20 20 ];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams 20 20 20 ];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
            
        case 'gaussian-surround'
            fitParams.paramNames = {'x','y','rfWidth','timelag','tau','surrAmp', 'surrWidth'};
            fitParams.paramDescriptions = {'RF x position (digit)','RF y position (PD)','RF width (std of gaussian)','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)','amplitude of the surround', 'width of the surround'};
            fitParams.paramIncDec = [1 1 1 0.1 0.5 0.5 1];
            fitParams.paramMin = [-inf -inf 0 0 0 -inf -inf];
            fitParams.paramMax = [inf inf inf inf inf inf inf];
            % set min/max and init
            fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0 -inf -inf];
            fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf 3 inf inf inf];
            fitParams.initParams = [0 0 4 fitParams.timelag fitParams.tau 0 0];
            % add on parameters for difference of gamma
            if fitParams.diffOfGamma
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            % nine params
        case 'nine-param'
            fitParams.rfType = 'nine-param';
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'a1','a2','a3','b1','b2','b3','c1','c2','c3'};
            fitParams.paramDescriptions = {'Weight (a1)','Weight (a2)','Weight (a3)', ...
                'Weight (b1)','Weight (b2)','Weight (b3)', ...
                'Weight (c1)','Weight (c2)','Weight (c3)'};
            fitParams.paramIncDec = [1 1 1   1 1 1  1 1 1];
            %fitParams.paramIncDec = [-0.5+rand(1,9)];
            fitParams.paramMin = [-inf -inf -inf   -inf -inf -inf   -inf -inf -inf];
            fitParams.paramMax = [inf inf inf   inf inf inf   inf inf inf];
            % set min/max and init
            fitParams.minParams = fitParams.stimExtents([ones(1,9), 2*ones(1,9)]) ;
            fitParams.maxParams = fitParams.stimExtents([3*ones(1,9), 4*ones(1,9)]);
            fitParams.initParams = [ones(1,9)];
            %fitParams.initParams = [-0.5+rand(1,9)];
        case 'nine-param-hdr'
            fitParams.rfType = 'nine-param';
            
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'a1','a2','a3','b1','b2','b3','c1','c2','c3','timelag','tau'};
            fitParams.paramDescriptions = {'Weight (a1)','Weight (a2)','Weight (a3)', ...
                'Weight (b1)','Weight (b2)','Weight (b3)', ...
                'Weight (c1)','Weight (c2)','Weight (c3)','Time before start of rise of hemodynamic function', 'Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1  1 1 1  1 1 1 0.1 0.5];
            %fitParams.paramIncDec = [-0.5+rand(1,9)];
            %fitParams.paramMin = [-inf -inf -inf  -inf -inf -inf  -inf -inf -inf 0 0];
            %fitParams.paramMax = [inf inf inf  inf inf inf  inf inf inf  inf inf];
            fitParams.paramMin = [zeros(1,9), 0 0];
            fitParams.paramMax = [4.*ones(1,9) 10 10];
            
            % set min/max and init
            %fitParams.minParams = [fitParams.stimExtents([ones(1,9), 2*ones(1,9)]) 0 0] ;
            %fitParams.maxParams = [fitParams.stimExtents([3*ones(1,9), 4*ones(1,9)]) 10  10 ];
            fitParams.minParams = [ones(1,9), 0 0];
            fitParams.maxParams = [3*ones(1,9), 10 10];
            
            % random seed
            %fitParams.initParams = [[0.5+3 .* rand(1,9)] fitParams.timelag fitParams.tau];
            fitParams.initParams = [ones(1,9) fitParams.timelag fitParams.tau];
            %fitParams.initParams = [-0.5+rand(1,9) fitParams.timelag fitParams.tau];
            % add on parameters for difference of gamma
            
            if fitParams.diffOfGamma
                
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' 20 20 20];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams 20 20 20];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
            % Want to implement a 1D Gaussian ALONG a digit, across maybe 5
            % discrete sites. That would be 15 stim sites (5x3). We want 6
            % params, a mean and SD for the gaussian at each site.
            
        case {'gaussian-1D','gaussian-1D-transpose'}
            
            if  size(fitParams.stimX,1) == 4
                % Not PROPERLY implemented yet
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {'mean1', 'mean2', 'mean3', 'mean4'...
                    'sd1', 'sd2', 'sd3','sd4',...
                    'amp2', 'amp3','amp4',...
                    'timelag','tau'};
                fitParams.paramDescriptions = {'mean of digit1', 'mean of digit2', 'mean of digit3', 'mean of digit4', ...
                    'SD1', 'SD2', 'SD3','SD4', ...
                    'amplitude2', 'amplitude3', 'amplitude4', ...
                    'Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
                fitParams.paramIncDec = [1 1 1 1  1 1 1 1  1 1 1  0.1 0.5];
                fitParams.paramMin = [0.5 0.5 0.5 0.5 0 0 0 0    -inf -inf -inf    0 0]; %Centers only go between 0 and 3 (in units of 'digit')
                fitParams.paramMax = [4.5 4.5 4.5 4.5 inf inf inf inf  inf inf inf   inf inf];
                % set min/max and init
                fitParams.minParams = [0.5 0.5 0.5 0.5 0 0 0 0  -inf -inf -inf   0 0];
                fitParams.maxParams = [4.5 4.5 4.5 4.5 inf inf inf inf  inf inf inf  3 inf];
                fitParams.initParams = [0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 fitParams.timelag fitParams.tau];
                
                % adding a random initialisation
                %fitParams.initParams = [ [0.5+4 .* rand(1,4)] [0+10 .* rand(1,7)] fitParams.timelag fitParams.tau];
                % adding a fixed initilaistaon
                %fitParams.initParams = [4.5 4.5 4.5 4.5 1 1 1 1 1 1 1 fitParams.timelag fitParams.tau];
                
                % add on parameters for difference of gamma
                if fitParams.diffOfGamma
                    %disp(sprintf('(diffOfGamma) !!! Not implemented yet!!!'));
                    % parameter names/descriptions and other information for allowing user to set them
                    fitParams.paramNames = {fitParams.paramNames{:} 'amplitude2' 'timelag2','tau2'};
                    fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                    fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                    fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                    fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                    % set min/max and init
                    fitParams.minParams = [fitParams.minParams 0 0 0];
                    fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                    fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
                end
                
            elseif  size(fitParams.stimX,1) == 3 %Gaussian 1D 3x3
                fitParams.paramNames = {'mean1', 'mean2', 'mean3'...
                    'sd1', 'sd2', 'sd3',...
                    'amp2', 'amp3',...
                    'timelag','tau'};
                fitParams.paramDescriptions = {'mean of digit1', 'mean of digit2', 'mean of digit3',  ...
                    'SD1', 'SD2', 'SD3', ...
                    'amplitude2', 'amplitude3', ...
                    'Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
                fitParams.paramIncDec = [1 1 1   1 1 1   1 1   0.1 0.5];
                fitParams.paramMin = [0.5 0.5 0.5 0 0 0   -inf -inf    0 0]; %Centers only go between 0 and 3 (in units of 'digit')
                fitParams.paramMax = [4.5 4.5 4.5  inf inf inf   inf inf    inf inf];
                % set min/max and init
                fitParams.minParams = [0.5 0.5 0.5 0 0 0  -inf -inf   0 0];
                fitParams.maxParams = [4.5 4.5 4.5  inf inf inf  inf inf   3 inf];
                fitParams.initParams = [0.5 0.5 0.5  1 1 1 1  1 fitParams.timelag fitParams.tau];
                % random initilaisation
                %fitParams.initParams = [ [0.5+3 .* rand(1,3)] [0+10 .* rand(1,5)] fitParams.timelag fitParams.tau];
                % add on parameters for difference of gamma
                if fitParams.diffOfGamma
                    %disp(sprintf('(diffOfGamma) !!! Not implemented yet!!!'));
                    % parameter names/descriptions and other information for allowing user to set them
                    fitParams.paramNames = {fitParams.paramNames{:} 'amplitude2' 'timelag2','tau2'};
                    fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                    fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                    fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                    fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                    % set min/max and init
                    fitParams.minParams = [fitParams.minParams 0 0 0];
                    fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                    fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
                end
                
            elseif size(fitParams.stimX,2) == 5 % pRF on patients TW
                
                %                 fitParams.paramNames = {'mean1', 'mean2', 'mean3', 'mean4', 'mean5',...
                %                     'sd1', 'sd2', 'sd3','sd4', 'sd5',...
                %                     'amp2', 'amp3','amp4','amp5',...
                %                     'timelag','tau'};
                fitParams.paramNames = {'mean1',...
                    'sd1', ...
                    'amp1',...
                    'timelag','tau'};
                %                 fitParams.paramDescriptions = {'mean of digit1', 'mean of digit2', 'mean of digit3', 'mean of digit4', 'mean of digit5',  ...
                %                     'SD1', 'SD2', 'SD3','SD4','SD5', ...
                %                     'amplitude2', 'amplitude3', 'amplitude4', 'amplitude5',...
                %                     'Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
                fitParams.paramDescriptions = {'mean of digit1',  ...
                    'SD1', ...
                    'amplitude1',...
                    'Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
                %                 fitParams.paramIncDec = [1 1 1 1 1   1 1 1 1 1   1 1 1 1   0.1 0.5];
                %                 fitParams.paramMin = [0.5 0.5 0.5 0.5 0.5  0 0 0 0 0   -inf -inf -inf -inf    0 0]; %Centers only go between 0 and 3 (in units of 'digit')
                %                 fitParams.paramMax = [5.5 5.5 5.5 5.5 5.5 inf inf inf inf inf   inf inf inf inf inf    inf inf];
                %                 % set min/max and init
                %                 fitParams.minParams = [0.5 0.5 0.5 0.5 0.5  0 0 0 0 0  -inf -inf -inf -inf   0 0];
                %                 fitParams.maxParams = [5.5 5.5 5.5 5.5 5.5  inf inf inf inf inf   inf inf inf inf  3 inf];
                %                 fitParams.initParams = [0.5 0.5 0.5 0.5 0.5  1 1 1 1 1   1 1 1 1 fitParams.timelag fitParams.tau];
                
                fitParams.paramIncDec = [1    1    1    0.1 0.5];
                fitParams.paramMin = [0.5   0    -inf     0 0];
                fitParams.paramMax = [5.5  inf   inf     inf inf];
                % set min/max and init
                fitParams.minParams = [0.5   0   -inf    0 0];
                fitParams.maxParams = [5.5   inf    inf   3 inf];
                fitParams.initParams = [0.5   1    1  fitParams.timelag fitParams.tau];
                
                % random initilaisation
                %fitParams.initParams = [ [0.5+3 .* rand(1,3)] [0+10 .* rand(1,5)] fitParams.timelag fitParams.tau];
                % add on parameters for difference of gamma
                if fitParams.diffOfGamma
                    %disp(sprintf('(diffOfGamma) !!! Not implemented yet!!!'));
                    % parameter names/descriptions and other information for allowing user to set them
                    fitParams.paramNames = {fitParams.paramNames{:} 'amplitude2' 'timelag2','tau2'};
                    fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                    fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                    fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                    fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                    % set min/max and init
                    fitParams.minParams = [fitParams.minParams 0 0 0];
                    fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                    fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
                end
                
                
                
                
            end
            
            
        case 'gaussian-1D-orthotips'
            
            %if  size(fitParams.stimX,1) == 4
            fitParams.paramNames = {'mean1', 'sd','amp','timelag','tau'};
            fitParams.paramDescriptions = {'mean of G1','SD1','amplitude1','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1  1  1   0.1 0.5];
            fitParams.paramMin = [0.5 -inf -inf  0 0]; %Centers only go between 0 and 3 (in units of 'digit')
            fitParams.paramMax = [5.5 inf inf  inf inf];
            % set min/max and init
            fitParams.minParams = [0.5 0 -inf  0 0];
            fitParams.maxParams = [5.5 inf inf 3 inf];
            fitParams.initParams = [0.5 1 1 fitParams.timelag fitParams.tau];
            
            if fitParams.diffOfGamma
                
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amplitude2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
            %else
            %                 fitParams.paramNames = {'mean1', 'sd1','amp1','timelag','tau'};
            %                 fitParams.paramDescriptions = {'mean of G1','SD1','amplitude1','Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
            %                 fitParams.paramIncDec = [1  1  1   0.1 0.5];
            %                 fitParams.paramMin = [0.5 -inf -inf  0 0]; %Centers only go between 0 and 3 (in units of 'digit')
            %                 fitParams.paramMax = [3.5 inf inf  inf inf];
            %                 % set min/max and init
            %                 fitParams.minParams = [0.5 -inf -inf  0 0];
            %                 fitParams.maxParams = [3.5 inf inf 3 inf];
            %                 fitParams.initParams = [0.5 1 1 fitParams.timelag fitParams.tau];
            %
            %                 if fitParams.diffOfGamma
            %                     %disp(sprintf('(diffOfGamma) !!! Not implemented yet!!!'));
            %                     % parameter names/descriptions and other information for allowing user to set them
            %                     fitParams.paramNames = {fitParams.paramNames{:} 'amplitude2' 'timelag2','tau2'};
            %                     fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
            %                     fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
            %                     fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
            %                     fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
            %                     % set min/max and init
            %                     fitParams.minParams = [fitParams.minParams 0 0 0];
            %                     fitParams.maxParams = [fitParams.maxParams inf 6 inf];
            %                     fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            %                 end
            
            %end
            
        case 'sixteen-hdr'
            fitParams.rfType = 'sixteen-hdr';
            
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'a1','a2','a3', 'a4','b1','b2','b3', 'b4','c1','c2','c3', ...
                'c4', 'd1', 'd2', 'd3', 'd4','timelag','tau'};
            fitParams.paramDescriptions = {'Weight (a1)','Weight (a2)','Weight (a3)', ...
                'Weight (a4)', 'Weight (b1)','Weight (b2)','Weight (b3)', ...
                'Weight (b4)', 'Weight (c1)','Weight (c2)','Weight (c3)',...
                'Weight (c4)', ...
                'Weight (d1', 'Weight (d2)', 'Weight (d3)', 'Weight (d4)'...
                'Time before start of rise of hemodynamic function', 'Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1 0.1 0.5];
            %fitParams.paramIncDec = [-0.5+rand(1,9)];
            %fitParams.paramMin = [-inf -inf -inf  -inf -inf -inf  -inf -inf -inf 0 0];
            %fitParams.paramMax = [inf inf inf  inf inf inf  inf inf inf  inf inf];
            fitParams.paramMin = [zeros(1,16) 0 0];
            fitParams.paramMax = [4.*ones(1,16) 10 10];
            
            % set min/max and init
            fitParams.minParams = [zeros(1,16) 0 0] ;
            fitParams.maxParams = [4.*ones(1,16) 10 10 ];
            fitParams.initParams = [ones(1,16) fitParams.timelag fitParams.tau];
            
            % try random seed
            %fitParams.initParams = [[1+3 .* rand(1,16)] fitParams.timelag fitParams.tau];
            % try max seed
            %fitParams.initParams = [4.*ones(1,16) fitParams.timelag fitParams.tau];
            
            % add on parameters for difference of gamma
            
            if fitParams.diffOfGamma
                
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' 20 20 20];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams 20 20 20];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
            %         case 'gaussian-tips'
            %             fitParams.paramNames = {'mean','amp','rfwidth'};
            %             fitParams.paramDescriptions = {'mean','amp','rfwidth'};
            %             fitParams.paramIncDec = [1 1 1 ];
            %             fitParams.paramMin = [0 0 0]; %Centers only go between 0 and 3 (in units of 'digit')
            %             fitParams.paramMax = [inf inf inf];
            %             % set min/max and init
            %             fitParams.minParams = [ 0 0 0 ];
            %             fitParams.maxParams = [inf inf inf];
            %             fitParams.initParams = [1 1 1];
            %         case 'gaussian-tips-hdr'
            %             fitParams.paramNames = {'mean','amp','rfwidth', 'timelag', 'tau'};
            %             fitParams.paramDescriptions = {'mean','amp','rfwidth', 'timlag', 'tau'};
            %             fitParams.paramIncDec = [1 1 1 0.1 0.5 ];
            %             fitParams.paramMin = [0 0 0 0 0]; %Centers only go between 0 and 3 (in units of 'digit')
            %             fitParams.paramMax = [inf inf inf inf inf];
            %             % set min/max and init
            %             fitParams.minParams = [ 0 0 0 0 0];
            %             fitParams.maxParams = [inf inf inf inf inf];
            %             fitParams.initParams = [1 1 1 fitParams.timelag fitParams.tau];
            %             % add on parameters for difference of gamma
            %             if fitParams.diffOfGamma
            %                 % parameter names/descriptions and other information for allowing user to set them
            %                 fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
            %                 fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
            %                 fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
            %                 fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
            %                 fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
            %                 % set min/max and init
            %                 fitParams.minParams = [fitParams.minParams 0 0 0];
            %                 fitParams.maxParams = [fitParams.maxParams inf 6 inf];
            %                 fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            %             end
            
        case 'five-hdr'
            fitParams.rfType = 'five-hdr';
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'a1','a2','a3','a4', 'a5', 'timelag','tau'};
            fitParams.paramDescriptions = {'Weight (a1)','Weight (a2)','Weight (a3)', ...
                'Weight (a4)','Weight (a5)', ...
                'Time before start of rise of hemodynamic function', ...
                'Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1 1 1  0.1 0.5];
            fitParams.paramMin = [zeros(1,5), 0 0];
            fitParams.paramMax = [5.*ones(1,5) inf inf];
            % set min/max and init
            fitParams.minParams = [1 1 1 1 1 0 0] ;
            fitParams.maxParams = [5 5 5 5 5 3 inf ];
            fitParams.initParams = [ones(1,5) fitParams.timelag fitParams.tau];
            %randInit = 1 + (5-1).*rand(5,1);
            %randInit = randInit';
            %fitParams.initParams = [randInit fitParams.timelag fitParams.tau];
            if fitParams.diffOfGamma
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
        case 'four-hdr'
            fitParams.rfType = 'four-hdr';
            % parameter names/descriptions and other information for allowing user to set them
            fitParams.paramNames = {'a1','a2','a3','a4', 'timelag','tau'};
            fitParams.paramDescriptions = {'Weight (a1)','Weight (a2)','Weight (a3)', ...
                'Weight (a4)', ...
                'Time before start of rise of hemodynamic function', ...
                'Width of the hemodynamic function (tau parameter of gamma)'};
            fitParams.paramIncDec = [1 1 1 1   0.1 0.5];
            fitParams.paramMin = [zeros(1,4), 0 0];
            fitParams.paramMax = [4.*ones(1,4) inf inf];
            % set min/max and init
            fitParams.minParams = [1 1 1 1  0 0] ;
            fitParams.maxParams = [4 4 4 4 3 inf ];
            fitParams.initParams = [ones(1,4) fitParams.timelag fitParams.tau];
            if fitParams.diffOfGamma
                % parameter names/descriptions and other information for allowing user to set them
                fitParams.paramNames = {fitParams.paramNames{:} 'amp2' 'timelag2','tau2'};
                fitParams.paramDescriptions = {fitParams.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
                fitParams.paramIncDec = [fitParams.paramIncDec(:)' 0.1 0.1 0.5];
                fitParams.paramMin = [fitParams.paramMin(:)' 0 0 0];
                fitParams.paramMax = [fitParams.paramMax(:)' inf inf inf];
                % set min/max and init
                fitParams.minParams = [fitParams.minParams 0 0 0];
                fitParams.maxParams = [fitParams.maxParams inf 6 inf];
                fitParams.initParams = [fitParams.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
            end
            
            
            
            
        otherwise
            disp(sprintf('(pRF_somatoFit:setFitParams) Unknown rfType %s',rfType));
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
                disp(sprintf('(pRF_somatoFit) Using default constraints'));
            else
                % get the parameter constraints back from the dialog entries
                for iParam = 1:length(fitParams.paramNames)
                    fitParams.minParams(iParam) = params.(sprintf('min%s',fitParams.paramNames{iParam}));
                    fitParams.maxParams(iParam) = params.(sprintf('max%s',fitParams.paramNames{iParam}));
                end
            end
        end
        % Now display parameter constraints
        for iParam = 1:length(fitParams.paramNames)
            disp(sprintf('(pRF_somatoFit) Parameter %s [min:%f max:%f] (%i:%s)',fitParams.paramNames{iParam},fitParams.minParams(iParam),fitParams.maxParams(iParam),iParam,fitParams.paramDescriptions{iParam}));
        end
    else
        % no constraints allowed
        disp(sprintf('(pRF_somatoFit) !!! Fit constraints ignored for algorithm: %s (if you want to constrain the fits, then use: %s) !!!',fitParams.algorithm,cell2mat(algorithmsWithConstraints)));
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getModelResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual, modelResponse, rfModel, r, hrf] = getModelResidual(params,tSeries,fitParams,justGetModel, hrfprfcheck)

%residual = [];
if nargin < 4, justGetModel = 0;end


% get the model response
% convert parameter array into a parameter strucutre
p = getFitParams(params,fitParams);

% compute an RF
rfModel = getRFModel(p,fitParams);
% % % include somato model here:
% rfModel = getSomatoRFModel(p, fitParams);

%tempfix for crossVal
if ieNotDefined('hrfprfcheck')
    hrfprfcheck = 0;
end

% if crossValcheck == 1
% rfModel = rfModel';
% else
% end
% init model response
modelResponse = [];residual = [];

% create the model for each concat
for i = 1:fitParams.concatInfo.n
    % get model response
    thisModelResponse = convolveModelWithStimulus(rfModel,fitParams.stim{i});
    
    % get a model hrf
    
    if isfield(fitParams, 'hrfprf')
        hrf.hrf = fitParams.hrfprf;
        hrf.time = 0:fitParams.framePeriod:p.canonical.lengthInSeconds; % this if 24 seconds, tr 2s fyi...
        % normalize to amplitude of 1
        hrf.hrf = hrf.hrf / max(hrf.hrf);
    else
        hrf = getCanonicalHRF(p.canonical,fitParams.framePeriod);
    end
    
    
    %hrf = getCanonicalHRF(p.canonical,fitParams.framePeriod);
    
    % and convolve in time.
    thisModelResponse = convolveModelResponseWithHRF(thisModelResponse,hrf);
    
    % drop junk frames here
    thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);
    
    % apply concat filtering
    if isfield(fitParams,'applyFiltering') && fitParams.applyFiltering
        thisModelResponse = applyConcatFiltering(thisModelResponse,fitParams.concatInfo,i);
    else
        % with no filtering, just remove mean
        thisModelResponse = thisModelResponse - mean(thisModelResponse);
    end
    
    %if ~justGetModel
    %if isempty(justGetModel)
    if isempty(justGetModel)
        justGetModel = 0;
    end
    
    if justGetModel == 0
        % compute correlation of this portion of the model response with time series
        thisTSeries = tSeries(fitParams.concatInfo.runTransition(i,1):fitParams.concatInfo.runTransition(i,2));
        thisTSeries = thisTSeries - mean(thisTSeries);
        
        % check here for length
        if length(thisTSeries) ~= length(thisModelResponse)
            disp(sprintf('(pRFFit:getModelResidual) Voxel tSeries length of %i does not match model length of %i. This can happen, for instance, if the tSense factor was not set correctly or junk frames was not set correctly.',length(thisTSeries),length(thisModelResponse)));
            keyboard
        end
        
        r(i) = corr(thisTSeries(:),thisModelResponse(:));
        
        if fitParams.betaEachScan
            % scale and offset the model to best match the tSeries
            [thisModelResponse thisResidual] = scaleAndOffset(thisModelResponse',thisTSeries(:));
        else
            thisResidual = [];
        end
    else
        thisResidual = [];
    end
    
    % make into a column array
    modelResponse = [modelResponse;thisModelResponse(:)];
    residual = [residual;thisResidual(:)];
end

% return model only
if justGetModel,return,end

% scale the whole time series
if ~isfield(fitParams, 'hrfprf')
    if ~fitParams.betaEachScan
        [modelResponse residual] = scaleAndOffset(modelResponse,tSeries(:));
    end
end


% display the fit
if fitParams.dispFit
    dispModelFit(params,fitParams,modelResponse,tSeries,rfModel);
end


% scale and offset (manual)
if fitParams.getModelResponse == 1
    
    if ~any(isnan(modelResponse))
        mref = mean(tSeries);
        stdRef = std(tSeries);
        mSig = mean(modelResponse);
        stdSig = std(modelResponse);
        modelResponse = ((modelResponse - mSig)/stdSig) * stdRef + mref;
        
        residual = tSeries-modelResponse;
    else
        residual = tSeries;
    end
    
elseif fitParams.getModelResponse ~= 1
    if hrfprfcheck == 1
        if ~any(isnan(modelResponse))
            %     warning('off', 'MATLAB:rankDeficientMatrix');
            X = modelResponse(:);
            X(:,2) = 1;
            
            b = X \ tSeries; % backslash linear regression
            %b = pinv(X) * tSeries;
            modelResponse = X * b;
            residual = tSeries-modelResponse;
        else
            residual = tSeries;
        end
        %modelResponse = newSig;
    end
end



% for nelder-mead just compute correlation and return 1-4
if strcmp(lower(fitParams.algorithm),'nelder-mead')
    residual = -corr(modelResponse,tSeries);
    %  disp(sprintf('(pRFFit:getModelResidual) r: %f',residual));
end

end


%%%%%%%%%%%%%%%%%%%%%%
%%    dispModelFit    %
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

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%    scaleAndOffset    %
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
end

%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)

p.rfType = fitParams.rfType;

switch (fitParams.rfType)
    case 'gaussian'
        p.x = params(1);
        p.y = params(2);
        p.std = params(3);
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
    case 'gaussian-hdr'
        p.x = params(1);
        p.y = params(2);
        p.std = params(3);
        % use a fixed single gaussian
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(4);
        p.canonical.tau = params(5);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(6);
            p.canonical.timelag2 = params(7);
            p.canonical.tau2 = params(8);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
    case 'gaussian-hdr-double'
        p.x = params(1);
        p.y = params(2);
        p.stdx = params(3);
        p.stdy = params(4);
        % use a fixed single gaussian
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(5);
        p.canonical.tau = params(6);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(7);
            p.canonical.timelag2 = params(8);
            p.canonical.tau2 = params(9);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
    case 'gaussian-surround'
        p.x = params(1);
        p.y = params(2);
        p.std = params(3);
        p.surrAmp = params(6);
        p.surrWidth = params(7);
        % use a fixed single gaussian
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(4);
        p.canonical.tau = params(5);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(8);
            p.canonical.timelag2 = params(9);
            p.canonical.tau2 = params(10);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
    case 'nine-param'
        p.weights = params(1:9);
        %         [~, p.x] = max(mean(params([1 2 3; 4 5 6; 7 8 9]),1));
        %         [~, p.y] = max(mean(params([1 2 3; 4 5 6; 7 8 9]),2));
        %
        %         [~, p.x] = max(mean(params([1 2 3; 4 5 6; 7 8 9])));
        %         [~, p.y] = max(mean(params([1 2 3; 4 5 6; 7 8 9]')));
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
    case 'nine-param-hdr'
        p.weights = params(1:9);
        %         [~, p.x] = max(mean(params([1 2 3; 4 5 6; 7 8 9]),1));
        %         [~, p.y] = max(mean(params([1 2 3; 4 5 6; 7 8 9]),2));
        %         [~, p.x] = max(mean(params([1 2 3; 4 5 6; 7 8 9])));
        %         [~, p.y] = max(mean(params([1 2 3; 4 5 6; 7 8 9]')));
        %         p.std = var(params(:));
        %
        
        
        %pRFweights = p.weights;
        
        %[com, momentOfInertia ]  = centerOfMass(pRFweights);
        
        %p.x = com(1);
        %p.y = com(2);
        %p.std = momentOfInertia(1);
        
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(10);
        p.canonical.tau = params(11);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(12);
            p.canonical.timelag2 = params(13);
            p.canonical.tau2 = params(14);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
    case {'gaussian-1D','gaussian-1D-transpose'}
        
        if size(fitParams.stimX,1) == 4
            p.meanOne = params(1);
            p.meanTwo = params(2);
            p.meanThr = params(3);
            p.meanFour = params(4);
            p.amp1 = 1;
            p.amp2 = params(9);
            p.amp3 = params(10);
            p.amp4 = params(11);
            p.stdOne = params(5);
            p.stdTwo = params(6);
            p.stdThr = params(7);
            p.stdFour = params(8);
            
            % use a fixed single gaussian
            p.canonical.type = 'gamma';
            p.canonical.lengthInSeconds = 25;
            p.canonical.timelag = params(12);
            p.canonical.tau = params(13);
            p.canonical.exponent = fitParams.exponent;
            p.canonical.offset = 0;
            p.canonical.diffOfGamma = fitParams.diffOfGamma;
            if fitParams.diffOfGamma
                p.canonical.amplitudeRatio = params(14);
                p.canonical.timelag2 = params(15);
                p.canonical.tau2 = params(16);
                p.canonical.exponent2 = fitParams.exponent2;
                p.canonical.offset2 = 0;
            end
            
        elseif size(fitParams.stimX,1) == 3
            p.meanOne = params(1);
            p.meanTwo = params(2);
            p.meanThr = params(3);
            
            p.amp1 = 1;
            p.amp2 = params(7);
            p.amp3 = params(8);
            
            p.stdOne = params(4);
            p.stdTwo = params(5);
            p.stdThr = params(6);
            
            
            % use a fixed single gaussian
            p.canonical.type = 'gamma';
            p.canonical.lengthInSeconds = 25;
            p.canonical.timelag = params(9);
            p.canonical.tau = params(10);
            p.canonical.exponent = fitParams.exponent;
            p.canonical.offset = 0;
            p.canonical.diffOfGamma = fitParams.diffOfGamma;
            if fitParams.diffOfGamma
                p.canonical.amplitudeRatio = params(11);
                p.canonical.timelag2 = params(12);
                p.canonical.tau2 = params(13);
                p.canonical.exponent2 = fitParams.exponent2;
                p.canonical.offset2 = 0;
            end
            
        elseif size(fitParams.stimX,2) == 5
            %             p.meanOne = params(1);
            %             p.meanTwo = params(2);
            %             p.meanThr = params(3);
            %             p.meanFour = params(4);
            %             p.meanFive = params(5);
            %             p.amp1 = 1;
            %             p.amp2 = params(11);
            %             p.amp3 = params(12);
            %             p.amp4 = params(13);
            %             p.amp5 = params(14);
            %             p.stdOne = params(6);
            %             p.stdTwo = params(7);
            %             p.stdThr = params(8);
            %             p.stdFour = params(9);
            %             p.stdFive = params(10);
            %
            %             % use a fixed single gaussian
            %             p.canonical.type = 'gamma';
            %             p.canonical.lengthInSeconds = 25;
            %             p.canonical.timelag = params(15);
            %             p.canonical.tau = params(16);
            %             p.canonical.exponent = fitParams.exponent;
            %             p.canonical.offset = 0;
            %             p.canonical.diffOfGamma = fitParams.diffOfGamma;
            %             if fitParams.diffOfGamma
            %                 p.canonical.amplitudeRatio = params(17);
            %                 p.canonical.timelag2 = params(18);
            %                 p.canonical.tau2 = params(19);
            %                 p.canonical.exponent2 = fitParams.exponent2;
            %                 p.canonical.offset2 = 0;
            %             end
            p.meanOne = params(1);
            p.amp1 = params(3);
            p.stdOne = params(2);
            
            % use a fixed single gaussian
            p.canonical.type = 'gamma';
            p.canonical.lengthInSeconds = 25;
            p.canonical.timelag = params(4);
            p.canonical.tau = params(5);
            p.canonical.exponent = fitParams.exponent;
            p.canonical.offset = 0;
            p.canonical.diffOfGamma = fitParams.diffOfGamma;
            if fitParams.diffOfGamma
                p.canonical.amplitudeRatio = params(6);
                p.canonical.timelag2 = params(7);
                p.canonical.tau2 = params(8);
                p.canonical.exponent2 = fitParams.exponent2;
                p.canonical.offset2 = 0;
            end
            
            
        end
        
    case 'sixteen-hdr'
        p.weights = params(1:16);
        %[~, p.x] = max(mean(params([1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16])));
        
        %[~, p.y] = max(mean(params([1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]')));
        %p.std = var(params(1:16));
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(17);
        p.canonical.tau = params(18);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(19);
            p.canonical.timelag2 = params(20);
            p.canonical.tau2 = params(21);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
        
    case 'gaussian-1D-orthotips'
        
        p.meanOne = params(1);
        p.std = params(2);
        p.amp = params(3);
        p.y = 1;
        
        % use a fixed single gaussian
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(4);
        p.canonical.tau = params(5);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(6);
            p.canonical.timelag2 = params(7);
            p.canonical.tau2 = params(8);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
        %
        %     case 'gaussian-tips'
        %         p.meanOne = params(1);
        %         p.amp = params(2);
        %         p.std = params(3);
        %         p.y = 1;
        %         % use a fixed single gaussian
        %         p.canonical.type = 'gamma';
        %         p.canonical.lengthInSeconds = 25;
        %         p.canonical.timelag = fitParams.timelag;
        %         p.canonical.tau = fitParams.tau;
        %         p.canonical.exponent = fitParams.exponent;
        %         p.canonical.offset = 0;
        %         p.canonical.diffOfGamma = fitParams.diffOfGamma;
        %         p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
        %         p.canonical.timelag2 = fitParams.timelag2;
        %         p.canonical.tau2 = fitParams.tau2;
        %         p.canonical.exponent2 = fitParams.exponent2;
        %         p.canonical.offset2 = 0;
        %
        %     case 'gaussian-tips-hdr'
        %         p.meanOne = params(1);
        %         p.amp = params(2);
        %         p.std = params(3);
        %         p.y = 1;
        %         % use a fixed single gaussian
        %         p.canonical.type = 'gamma';
        %         p.canonical.lengthInSeconds = 25;
        %         p.canonical.timelag = params(4);
        %         p.canonical.tau = params(5);
        %         p.canonical.exponent = fitParams.exponent;
        %         p.canonical.offset = 0;
        %         p.canonical.diffOfGamma = fitParams.diffOfGamma;
        %         if fitParams.diffOfGamma
        %             p.canonical.amplitudeRatio = params(6);
        %             p.canonical.timelag2 = params(7);
        %             p.canonical.tau2 = params(8);
        %             p.canonical.exponent2 = fitParams.exponent2;
        %             p.canonical.offset2 = 0;
        %         end
        %
    case 'five-hdr'
        p.weights = params(1:5);
        %[~, p.x] = max(params([1 2 3 4 5]));
        %[~, p.y] = max(mean(params([1 2 3 4 5]')));
        %p.x = max(mean(params([1 2 3 4 5])));
        %p.y = 1;
        %p.std = var(params(1:5));
        
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(6);
        p.canonical.tau = params(7);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(8);
            p.canonical.timelag2 = params(9);
            p.canonical.tau2 = params(10);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
        
    case 'four-hdr'
        p.weights = params(1:4);
        [~, p.x] = max(params([1 2 3 4 ]));
        %[~, p.y] = max(mean(params([1 2 3 4 5]')));
        %p.x = max(mean(params([1 2 3 4 5])));
        p.y = 1;
        p.std = var(params(1:4));
        
        p.canonical.type = 'gamma';
        p.canonical.lengthInSeconds = 25;
        p.canonical.timelag = params(5);
        p.canonical.tau = params(6);
        p.canonical.exponent = fitParams.exponent;
        p.canonical.offset = 0;
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.amplitudeRatio = params(7);
            p.canonical.timelag2 = params(8);
            p.canonical.tau2 = params(9);
            p.canonical.exponent2 = fitParams.exponent2;
            p.canonical.offset2 = 0;
        end
        
        
    otherwise
        disp(sprintf('(pRFFit) Unknown rfType %s',rfType));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = convolveModelWithStimulus(rfModel,stim)

% get number of frames
nFrames = size(stim.im,3);

% preallocate memory
modelResponse = zeros(1,nFrames);

% check matrix dims
if isequal(size(rfModel), size(stim.im(:,:,1)) )
    
    for frameNum = 1:nFrames
        % multipy the stimulus frame by frame with the rfModel
        % and take the sum
        modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
    end
    
elseif ~isequal(size(rfModel), size(stim.im(:,:,1)) )
    
    rfModel = rfModel';
    
    for frameNum = 1:nFrames
        % multipy the stimulus frame by frame with the rfModel
        % and take the sum
        modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

end

%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
    fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

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

end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);

end

%%%%%%%%%%%%%%%%%%%%
%%   getRFModel   %%
%%%%%%%%%%%%%%%%%%%%
function rfModel = getRFModel(params,fitParams)

rfModel = [];

% now gernerate the rfModel
switch fitParams.rfType
    case {'gaussian','gaussian-hdr','gaussian-hdr-double','gaussian-1D','gaussian-1D-transpose','gaussian-surround', 'gaussian-1D-orthotips'}
        rfModel = makeRFGaussian(params,fitParams);
    case {'nine-param','nine-param-hdr'}
        rfModel = makeRFNineParam(params,fitParams);
    case {'sixteen-hdr'}
        rfModel = makeRFSixteenParam(params, fitParams);
    case {'five-hdr', 'four-hdr'}
        rfModel = makeRFFiveParam(params, fitParams);
        %rfModel = params.weights';
    otherwise
        disp(sprintf('(pRFFit:getRFModel) Unknown rfType: %s',fitParams.rfType));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeRFGaussian   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFGaussian(params,fitParams)

% compute rf
switch fitParams.rfType
    case {'gaussian-1D','gaussian-1D-transpose'}
        %params.x = max(params.x);
        %rfModel = exp(-(((fitParams.stimX-params.x).^2)/(2*(params.std^2))+((fitParams.stimY-params.y).^2)/(2*(params.std^2))));
        %oneD = fitParams.stimX(:,1);
        
        % Want to implement a one dimensional gaussian along a digit
        % This will give a 1x3 matrix, need a 3x3 for the stimulus
        % convolution. Not sure if this makes sense...
        %rfModel = normpdf(oneD, mean(oneD), std(oneD));
        %rfModel = [rfModel rfModel rfModel];
        %x = zeros(3);
        %rfModel = [x(:,1) rfModel x(:,1)]; %digit 2
        %rfModel = [rfModel x(:,1:2)]; % digit1
        %rfModel = [x(:,1:2), rfModel]; % digit3
        %rfModel = normpdf(fitParams.stimX, params.mean, params.std);
        
        % Usage: p = [height center SD offset]
        
        if size(fitParams.stimX,1) == 4
            X = linspace(1,4,4);
            params.amp1 = 1;
            pone = [params.amp1 params.meanOne params.stdOne 0];
            ptwo = [params.amp2 params.meanTwo params.stdTwo 0];
            pthr = [params.amp3 params.meanThr params.stdThr 0];
            pFour = [params.amp4 params.meanFour params.stdFour 0];
            %[R,Y] = meshgrid(1:4,1:4);
            %
            % this WAS the version % Z = [gauss(pone,X); gauss(ptwo,X); gauss(pthr,X); gauss(pFour,X)];
            % @DS suggests this .. this means each digit inhabits a COLUMN
            % in the model and should be consistent with pRFStimImage
            % convention.
            Z = [gauss(pone,X)', gauss(ptwo,X)', gauss(pthr,X)', gauss(pFour,X)'];
            
            if strcmpi(fitParams.rfType, 'gaussian-1D-transpose')
                %rfModel = (Z); %reality check (flip base to tip) using flipud
                rfModel = transpose(Z); % try for orthoGaussian
            else
                rfModel = Z;
            end
            
            %surf(R,Y,Z)
        elseif size(fitParams.stimX,1) == 3
            X = linspace(1,3,3);
            params.amp1 = 1;
            pone = [params.amp1 params.meanOne params.stdOne 0];
            ptwo = [params.amp2 params.meanTwo params.stdTwo 0];
            pthr = [params.amp3 params.meanThr params.stdThr 0];
            
            %[R,Y] = meshgrid(1:3,1:3);
            Z = [gauss(pone,X)', gauss(ptwo,X)', gauss(pthr,X)'];
            %rfModel = (Z);
            if strcmpi(fitParams.rfType, 'gaussian-1D-transpose')
                %rfModel = (Z); %reality check (flip base to tip) using flipud
                rfModel = transpose(Z); % try for orthoGaussian
            else
                rfModel = Z;
            end
            %rfModel = transpose(Z);
            % add a transpose here, if we want to do 1D in the orthogonal
            % direction!
            
            
        elseif size(fitParams.stimX,2) == 5
%             X = linspace(1,5,5);
%             params.amp1 = 1;
%             pone = [params.amp1 params.meanOne params.stdOne 0];
%             ptwo = [params.amp2 params.meanTwo params.stdTwo 0];
%             pthr = [params.amp3 params.meanThr params.stdThr 0];
%             pFour = [params.amp4 params.meanFour params.stdFour 0];
%             pFive = [params.amp5 params.meanFive params.stdFive 0];
%             Z = [gauss(pone,X)', gauss(ptwo,X)', gauss(pthr,X)', gauss(pFour,X)', gauss(pFive,X)'];
%             
%             if strcmpi(fitParams.rfType, 'gaussian-1D-transpose')
%                 %rfModel = (Z); %reality check (flip base to tip) using flipud
%                 rfModel = transpose(Z); % try for orthoGaussian
%             else
%                 rfModel = Z;
%             end
            X = linspace(1,5,5);
            %params.amp1 = 1;
            pone = [params.amp1 params.meanOne params.stdOne 0];
            %ptwo = [params.amp2 params.meanTwo params.stdTwo 0];
            %pthr = [params.amp3 params.meanThr params.stdThr 0];
            %pFour = [params.amp4 params.meanFour params.stdFour 0];
            %pFive = [params.amp5 params.meanFive params.stdFive 0];
            %Z = [gauss(pone,X)', gauss(ptwo,X)', gauss(pthr,X)', gauss(pFour,X)', gauss(pFive,X)'];
            Z = gauss(pone,X)';
            
            if strcmpi(fitParams.rfType, 'gaussian-1D-transpose')
                %rfModel = (Z); %reality check (flip base to tip) using flipud
                rfModel = transpose(Z); % try for orthoGaussian
            else
                rfModel = Z;
            end
            
            
            
            
            
            
        end
        
        
    case {'gaussian-surround'}
        % adds inhibitory surround to the 2D gaussian
        p = [1.5 params.x params.y sigma1 sigma2];
        [ZDOG, X, Y] = dogFit(1:3,1:3,p);
        rfModel = ZDOG;
        %surf(X,Y,ZDOG)
        
    case {'gaussian-1D-orthotips'} % CHECK THIS!
        
        mylen = length(fitParams.stimX);
        
        % if size(fitParams.stimX,1) == 4
        X = 1:mylen;
        pone = [params.amp params.meanOne params.std 0];
        Z = gauss(pone,X);
        rfModel = Z; %reality check (flip base to tip) using flipud
        %surf(R,Y,Z)
        %         else
        %             X = linspace(1,1,mylen);
        %             params.amp1 = 1;
        %             pone = [params.amp1 params.meanOne params.stdOne 0];
        %             %[R,Y] = meshgrid(1:3,1:3);
        %             Z = gauss(pone,X)';
        %             rfModel = Z;
        %             % add a transpose here, if we want to do 1D in the orthogonal
        %             % direction!
        %         end
        
        %
        %     case {'gaussian-tips'}
        %         X = 1:5;
        %         params.y = 1;
        %         pone = [params.amp params.meanOne params.std 0];
        %         [R,Y] = meshgrid(1:5,1:1);
        %         Z = gauss(pone,X);
        %         rfModel = Z;
        %     case {'gaussian-tips-hdr'}
        %         X = 1:5;
        %         params.y = 1;
        %         pone = [params.amp params.meanOne params.std 0];
        %         [R,Y] = meshgrid(1:5,1:1);
        %         Z = gauss(pone,X);
        %         rfModel = Z;
        
        %figure; plot(R,Z)
        
    case {'gaussian-hdr-double'}
        rfModel = exp(-(((fitParams.stimX-params.x).^2)/(2*(params.stdx^2))+((fitParams.stimY-params.y).^2)/(2*(params.stdy^2))));
        % do we need to flip this top to bottom?
        
    otherwise
        rfModel = exp(-(((fitParams.stimX-params.x).^2)/(2*(params.std^2))+((fitParams.stimY-params.y).^2)/(2*(params.std^2))));
        % try adding std
        
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeSomatoPRF    %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFNineParam(params,fitParams)
% makeRFNineParam - turn parameters into a pRF (here, 3x3)
%
%  9 parameters -- all independent, needs fixing...
%
% this function makes an appropriately shaped pRF from parameters

% turn the list into a grid
rfModel = reshape(params.weights, [3 3]);

% other versions of this might take another list of params, pIn (e.g. x0,
% y0, sigma0) into a 3x3 pOut. It depends what shape we want to impose.
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeSomatoPRF16    %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFSixteenParam(params,fitParams)
% makeRFNineParam - turn parameters into a pRF (here,4x4)
%
%  16 parameters -- all independent, needs fixing...
%
% this function makes an appropriately shaped pRF from parameters

% turn the list into a grid
rfModel = reshape(params.weights, [4 4]);

% other versions of this might take another list of params, pIn (e.g. x0,
% y0, sigma0) into a 3x3 pOut. It depends what shape we want to impose.

end

function rfModel = makeRFFiveParam(params,fitParams)
rfModel = reshape(params.weights, [1 5]);
end

%%%%%%%%%%%%%%%%%%%
%%    parseArgs    %
%%%%%%%%%%%%%%%%%%%
function [v ,scanNum, x, y, s, fitParams, tSeries, hrfprf] = parseArgs(args)

v = [];scanNum=[];x=[];y=[];s=[];fitParams=[];tSeries = [];hrfprf = [];

% check for calling convention from interrogator
if (length(args) >= 7) && isnumeric(args{6})
    v = args{1};
    %overlayNum = args{2};
    scanNum = args{3};
    x = args{4};
    y = args{5};
    s = args{6};
    %roi = args{7};
    hrfprf = args{7};
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
    getArgs({args{6:end}},{'dispFit=0','stim=[]','getModelResponse=0','params=[]','concatInfo=[]','prefit=[]','xFlipStimulus=0','yFlipStimulus=0','timeShiftStimulus=0','rfType=gaussian','betaEachScan=0','fitTypeParams=[]','justGetStimImage=[]','verbose=1','dispIndex=[]','dispN=[]','returnPrefit=0','quickPrefit=0','tSeries=[]','junkFrames=[]','framePeriod=[]','paramsInfo=[]', 'hrfprf=[]'});
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
        disp(sprintf('(pRF_somatoFit) Must pass in params when using getModelResponse'));
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


if isempty(fitParams.prefit) || (fitParams.prefit.quickPrefit ~= fitParams.quickPrefit)
    % set the values over which to first prefit
    % the best of these parameters will then be used
    % to init the non-linear optimization. Note that the
    % values here are expressed as a factor of the screen
    % dimensions (1 being the width/height of the screen)
    % Later when the prefit is calculated, they will be multiplied
    % by the screenWidth and screenHeight
    
    % make sure here that x and y points go through 0 symmetrically
    %[prefitx prefity prefitrfHalfWidth prefithrfDelay] = ndgrid(1:1:3,1:1:3,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75], 6);
    
    %[prefitx prefity prefitrfHalfWidth] = [1 2 3 4 5, 1 1 1 1, 1];
    
    %[prefitx prefity prefitrfHalfWidth prefithrfDelay] = ndgrid(1:1:3,1:1:3,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75], 6);
    % [prefitx prefity prefitrfHalfWidth prefithrfDelay] = ...
    % ndgrid(1:1:4,1:1:4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75], 6);
    
    switch fitParams.rfType
%         case 'five-hdr'
%             if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
%             if fitParams.quickPrefit
%                 prefitx = [1 2 3 4 5];
%                 prefity = [1 1 1 1 1];
%                 prefitrfHalfWidth = [1 1 1 1 1];
%             else
%                 prefitx = [1 2 3 4 5];
%                 prefity = [1 1 1 1 1];
%                 prefitrfHalfWidth = [1 1 1 1 1];
%             end
            
        case 'four-hdr'
            if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
            if fitParams.quickPrefit
                prefitx = [1 2 3 4 ];
                prefity = [1 1 1 1 ];
                prefitrfHalfWidth = [1 1 1 1 ];
            else
                prefitx = [1 2 3 4 ];
                prefity = [1 1 1 1 ];
                prefitrfHalfWidth = [1 1 1 1 ];
            end
            
        case {'gaussian','gaussian-hdr','gaussian-hdr-double','gaussian-1D','gaussian-surround','gaussian-1D-transpose'}
            if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
            % set the values over which to first prefit
            % the best of these parameters will then be used
            % to init the non-linear optimization. Note that the
            % values here are expressed as a factor of the screen
            % dimensions (1 being the width/height of the screen)
            % Later when the prefit is calculated, they will be multiplied
            % by the screenWidth and screenHeight
            if fitParams.quickPrefit
                
                % make sure here that x and y points go through 0 symmetrically
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.375:0.125:0.375,-0.375:0.125:0.375,[0.025 0.05 0.15 0.4]);
            else
                
                % trying something here...
                % 16
                [prefitx prefity prefitrfHalfWidth] = ndgrid(0:0.1:5.5,0:0.1:5.5,[0 1 2 3 4 5 6]);
                % 9
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(0:0.1:3.5,0:0.1:3.5,[0 1 2 3 4 5 6]);
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
            end
            
        case{'gaussian-1D-orthotips'}
            
            [prefitx prefity prefitrfHalfWidth] = ndgrid(0:0.1:5.5, [1 1 1 1 1], [0 1 2 3 4 5 6]);
            
            
            
            %             if fitParams.quickPrefit
            %                 prefitx = [1 2 3 4];
            %                 prefity = [1 1 1 1];
            %                 prefitrfHalfWidth = [1 1 1 1];
            %             else
            %                 prefitx = [1 2 3 4 ];
            %                 prefity = [1 1 1 1 ];
            %                 prefitrfHalfWidth = [1 1 1 1 ];
            %             end
            
            % IGNORE THIS
            %             if fitParams.quickPrefit
            %                 prefitx = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];
            %                 prefity = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]'
            %                 prefitrfHalfWidth = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
            %             else
            %                 prefitx = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];
            %                 prefity = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]';
            %                 prefitrfHalfWidth = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
            %             end
            
        case {'sixteen-hdr','five-hdr'}
            if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
            % set the values over which to first prefit
            % the best of these parameters will then be used
            % to init the non-linear optimization. Note that the
            % values here are expressed as a factor of the screen
            % dimensions (1 being the width/height of the screen)
            % Later when the prefit is calculated, they will be multiplied
            % by the screenWidth and screenHeight
            if fitParams.quickPrefit
                
                % make sure here that x and y points go through 0 symmetrically
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.375:0.125:0.375,-0.375:0.125:0.375,[0.025 0.05 0.15 0.4]);
            else
                [prefitx prefity prefitrfHalfWidth] = ndgrid(0:0.1:5.5,0:0.1:5.5,[0 1 2 3 4 5 6]);
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
            end
            
            
            %             if fitParams.quickPrefit
            %                 prefitx = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];
            %                 prefity = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]';
            %                 prefitrfHalfWidth = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
            %             else
            %                 prefitx = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];
            %                 prefity = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]';
            %                 prefitrfHalfWidth = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
            %             end
        case 'nine-param-hdr'
            
            if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
            % set the values over which to first prefit
            % the best of these parameters will then be used
            % to init the non-linear optimization. Note that the
            % values here are expressed as a factor of the screen
            % dimensions (1 being the width/height of the screen)
            % Later when the prefit is calculated, they will be multiplied
            % by the screenWidth and screenHeight
            if fitParams.quickPrefit
                
                % make sure here that x and y points go through 0 symmetrically
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.375:0.125:0.375,-0.375:0.125:0.375,[0.025 0.05 0.15 0.4]);
            else
                prefitx = [1 2 3 ; 1 2 3 ; 1 2 3 ];
                prefity = [1 2 3 ; 1 2 3 ; 1 2 3 ]';
                prefitrfHalfWidth = [1 1 1 ; 1 1 1;  1 1 1;  1 1 1];
                %[prefitx, prefity, prefitrfHalfWidth] = ndgrid(0:0.1:4.5,0:0.1:4.5,[0 1 2 3 4 5 6]);
                %[prefitx prefity prefitrfHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
            end
            %             if fitParams.verbose,fprintf('\n(pRF_somatoFit) Doing quick prefit');end
            %             if fitParams.quickPrefit
            %                 prefitx = [1 2 3 ; 1 2 3 ; 1 2 3 ];
            %                 prefity = [1 2 3 ; 1 2 3 ; 1 2 3 ]';
            %                 prefitrfHalfWidth = [1 1 1 ; 1 1 1;  1 1 1;  1 1 1];
            %             else
            %                 prefitx = [1 2 3 ; 1 2 3 ; 1 2 3 ];
            %                 prefity = [1 2 3 ; 1 2 3 ; 1 2 3 ]';
            %                 prefitrfHalfWidth = [1 1 1 ; 1 1 1;  1 1 1;  1 1 1];
            %             end
            
    end
    
    fitParams.prefit.quickPrefit = fitParams.quickPrefit;
    fitParams.prefit.n = length(prefitx(:));
    fitParams.prefit.x = prefitx(:);
    fitParams.prefit.y = prefity(:);
    fitParams.prefit.rfHalfWidth = prefitrfHalfWidth(:);
    %fitParams.prefit.hrfDelay = prefithrfDelay(:);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    checkStimForAverages    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim ignoreMismatchStimfiles] = checkStimForAverages(v,scanNum,groupNum,stim,concatInfo,stimImageDiffTolerance)

ignoreMismatchStimfiles = false;

% this function will check for some bad casses (like concat of concats etc)
% it will also check that all the component scans of an average have the
% same stim image and warn if they do not. It will then replace the stim cell
% array for the average with a single stim file, so that processing
% can continue as normal for pRFFit

% if not a cell, then ok, return
if ~iscell(stim),return,end

% first check for bad shiftList or refverseLIst
p = viewGet(v,'params',scanNum,groupNum);
if isfield(p,'params') && isfield(p.params,'shiftList') && any(p.params.shiftList~=0)
    disp(sprintf('(pRF_somatoFit) Component scan %s:%i has a shiftList that is non-zero (%s). pRFFit does not handle non-zero shifts in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
    keyboard
end
if isfield(p,'params') && isfield(p.params,'reverseList') && any(p.params.reverseList~=0)
    disp(sprintf('(pRF_somatoFit) Component scan %s:%i has a reverseList that is non-zero (%s). pRFFit does not handle time-reversed time series in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
    keyboard
end

% if is a cell, check to see if this is a concat or not
if ~isempty(concatInfo) && (concatInfo.isConcat)
    % this is a concat, so check each one of the elements
    [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
    for i = 1:length(stim)
        % get concatInfo for original scan
        concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
        if ~isempty(concatInfo)
            disp(sprintf('(pRF_somatoFit:checkStimForAverages) Detected concatenation of concatenations. pRFFit not implemented yet to handle this'));
            stim = [];
            keyboard
            return;
        end
        % check this next scan
        [stim{i} ignoreMismatchStimfiles] = checkStimForAverages(v,originalScanNum(i),originalGroupNum(i),stim{i},concatInfo,stimImageDiffTolerance);
        % if user has accepted all then set stimImageDiffTOlerance to infinity
        if isinf(ignoreMismatchStimfiles),stimImageDiffTolerance = inf;end
        if isempty(stim{i}),stim = [];return,end
    end
else
    % this for orignals
    [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
    % if it is an original than check each element
    if ~isempty(originalScanNum)
        % check that this is not an average of a concat
        for i = 1:length(stim)
            % get concatInfo for original scan
            concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
            if ~isempty(concatInfo)
                disp(sprintf('(pRF_somatoFit:checkStimForAverages) Detected average of a concatenations. pRFFit not implemented yet to handle this'));
                keyboard
                stim = [];
                return;
            end
            % see if it is an average of an average
            originalOfOriginalScanNum = viewGet(v,'originalScanNum',originalScanNum(i),originalGroupNum(i));
            if length(originalOfOriginalScanNum) > 1
                disp(sprintf('(pRF_somatoFit:checkStimForAverages) Detected average of an average. pRFFit not implemented yet to handle this'));
                keyboard
                stim = [];
                return;
            end
        end
        % ok, not an average of a concatenation/average so check all the stim files
        % and warn if there are any inconsistencies
        for i = 1:length(stim)
            if ~isequalwithequalnans(stim{1}.im,stim{i}.im)
                dispHeader
                disp(sprintf('(pRF_somatoFit:checkStimForAverages) !!! Average for %s:%i component scan %i does not match stimulus for other scans. If you wish to continue then this will use the stimfile associated with the first scan in the average !!!',viewGet(v,'groupName',groupNum),scanNum,originalScanNum(i)));
                % display which volumes are different
                diffVols = [];
                for iVol = 1:size(stim{1}.im,3)
                    if ~isequalwithequalnans(stim{1}.im(:,:,iVol),stim{i}.im(:,:,iVol))
                        diffVols(end+1) = iVol;
                    end
                end
                disp(sprintf('(pRF_somatoFit) Stimulus files are different at %i of %i vols (%0.1f%%): %s',length(diffVols),size(stim{1}.im,3),100*length(diffVols)/size(stim{1}.im,3),num2str(diffVols)));
                if 100*(length(diffVols)/size(stim{1}.im,3)) < stimImageDiffTolerance
                    disp(sprintf('(pRF_somatoFit) This could be for minor timing inconsistencies, so igorning. Set stimImageDiffTolerance lower if you want to stop the code when this happens'));
                else
                    % ask user if they want to continue (only if there is a difference of more than 10 vols
                    ignoreMismatchStimfiles = askuser('Do you wish to continue',1);
                    if ~ignoreMismatchStimfiles
                        stim = [];
                        return;
                    end
                end
                dispHeader
            end
        end
        % if we passed the above, this is an average of identical
        % scans, so just keep the first stim image since they are all the same
        stim = stim{1};
    end
end

end

%%%%%%%%%%%%%%%%%
%%    getStim    %
%%%%%%%%%%%%%%%%%
function stim = getStim(v,scanNum,fitParams)

% get stimfile
stimfile = viewGet(v,'stimfile',scanNum);
% get volume to trigger ratio
volTrigRatio = viewGet(v,'auxParam','volTrigRatio',scanNum);
% check if global matches
groupNum = viewGet(v,'curGroup');
global gpRFFitStimImage
if (isfield(fitParams,'recomputeStimImage') && fitParams.recomputeStimImage) || isempty(gpRFFitStimImage) || (gpRFFitStimImage.scanNum ~= scanNum)  || (gpRFFitStimImage.groupNum ~= groupNum) || (gpRFFitStimImage.xFlip ~= fitParams.xFlipStimulus) || (gpRFFitStimImage.yFlip ~= fitParams.yFlipStimulus) || (gpRFFitStimImage.timeShift ~= fitParams.timeShiftStimulus)
    disp(sprintf('(pRF_somatoFit) Computing stim image'));
    % if no save stim then create one
    stim = pRFGetSomatoStimImageFromStimfile(stimfile,'volTrigRatio',volTrigRatio,'xFlip',fitParams.xFlipStimulus,'yFlip',fitParams.yFlipStimulus,'timeShift',fitParams.timeShiftStimulus,'verbose',fitParams.verbose,'saveStimImage',fitParams.saveStimImage,'recomputeStimImage',fitParams.recomputeStimImage);
    % check for averages
    stim = checkStimForAverages(v,scanNum,viewGet(v,'curGroup'),stim,fitParams.concatInfo,fitParams.stimImageDiffTolerance);
    if isempty(stim),return,end
    % make into cell array
    stim = cellArray(stim);
    % save stim image in global
    gpRFFitStimImage.scanNum = scanNum;
    gpRFFitStimImage.groupNum = groupNum;
    gpRFFitStimImage.xFlip = fitParams.xFlipStimulus;
    gpRFFitStimImage.yFlip = fitParams.yFlipStimulus;
    gpRFFitStimImage.timeShift = fitParams.timeShiftStimulus;
    gpRFFitStimImage.stim = stim;
else
    % otherwise load from global
    disp(sprintf('(pRF_somatoFit) Using precomputed stim image'));
    stim = gpRFFitStimImage.stim;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    applyConcatFiltering    %
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


end
%%%%%%%%%%%%%
%%   r2d   %%
%%%%%%%%%%%%%
% function degrees = r2d(angle)
%
% degrees = (angle/(2*pi))*360;
%
% % if larger than 360 degrees then subtract
% % 360 degrees
% while (sum(degrees>360))
%     degrees = degrees - (degrees>360)*360;
% end
%
% % if less than 360 degreees then add
% % 360 degrees
% while (sum(degrees<-360))
%     degrees = degrees + (degrees<-360)*360;
% end
