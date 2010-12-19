% getGlmStatistics.m
%
%      usage: [d, contrastBetas, T, F] = getGlmStatistics(d, params, verbose)
%         by: Julien Besle
%       date: 18/01/10
%    purpose: Fits a GLM model and computes T and F statistics. Replaces getr2.
%              $Id$
%
%       e.g.: d = getGlmStatistics(d)
%             d needs to have a stimulus convolution matrix
%             when params.nBootstrap > 1, residual bootstrapping is performed
%             and output in the extra dimension of n.ehdr, 
%             the first element of the sixth dimension of n.ehdr is always the non-bootstrapped result

function [d, contrastBetas, T, F] = getGlmStatistics(d, params, verbose, precision, computeEstimates, computeTtests)

approximate_value = 0;
T = [];
F = [];
contrastBetas = [];
%bootstrap_max_voxels_ = [60 80 100 120 140 160 180 200];
%bootstrap_max_voxels = round(120000/nFrames);  %empirical value that optimizes the speed of bootstrapping, probably depends on available memory
bootstrap_max_voxels = 120; %optimizes the bootsrap computing time
optimalXLength = 60;   %optimizes the OLS residual and MSS computing time


%--------------------------------------------------- DEFAULT VALUES ------------------------------------%
if ieNotDefined('verbose'),verbose = 1;end
if ieNotDefined('precision'),precision = mrGetPref('defaultPrecision');end
if ieNotDefined('computeTtests'),computeTtests = 1;end
if ieNotDefined('computeEstimates'),computeEstimates = 1;end
if ieNotDefined('params')
  params = struct;
end

%default parameters
if ~isfield(params,'nBootstrap')
   nBootstrap = 1;
else
   nBootstrap = params.nBootstrap;
end

% check that contrasts and fTests have the same number of columns as the
% design matrix
if ~isfield(params,'contrasts') || isempty(params.contrasts)
  params.contrasts = {};
else
  if size(d.scm,2)~=size(params.contrasts,2)*d.nHrfComponents
    mrErrorDlg( 'contrasts incompatible with number of EVs');
    return;
  end
end

if~isfield(params,'numberContrasts')
  params.numberContrasts = 0;
end
contrasts = cell(1,params.numberContrasts);
for iContrast = 1:params.numberContrasts %convert contrasts to the same format as restrictions
  contrasts{iContrast} = params.contrasts(iContrast,:);
end

if ~isfield(params,'tTestSide')
   params.tTestSide = 'Both';
end

if ~isfield(params,'fTestNames') || isempty(params.fTestNames)
   params.fTestNames = {};
end
  
if ~isfield(params,'restrictions') || isempty(params.restrictions)
   params.restrictions = {};
else
   if size(d.scm,2)~=size(params.restrictions{1},2)*d.nHrfComponents
      mrErrorDlg( 'F tests incompatible with number of EVs');
      return;
   end
end
restrictions = params.restrictions;

if ~isfield(d,'nHrfComponents')
  d.nHrfComponents = size(d.scm,2)/d.nhdr;
end
if ~isfield(params,'componentsToTest') || isempty(params.componentsToTest)
   params.componentsToTest = ones(1,d.nHrfComponents);
end

if ~isfield(params,'componentsCombination') || isempty(params.componentsCombination)
   params.componentsCombination = 'Add';
end
if strcmp(params.componentsCombination,'Or') && d.nHrfComponents==1
  params.componentsCombination = 'Add';  %make sure this is set to add if there is only one component to test
end

if ieNotDefined('params') || ~isfield(params,'covCorrection') || ~isfield(params,'correctionType') || strcmp(params.correctionType,'none')
   params.covCorrection = 0;
end
if params.covCorrection
   if  ~isfield(params,'correctionType') || isempty(params.correctionType)
      params.correctionType = 'preWhitening';
   end
   if  ~isfield(params,'covEstimation') || isempty(params.covEstimation)
      params.covEstimation = 'singleTukeyTapers';
   else
      params.covEstimation = params.covEstimation;
   end
   if  ~isfield(params,'covEstimationAreaSize') || isempty(params.covEstimationAreaSize)
      params.covEstimationAreaSize= 1;
   else
      params.covEstimationAreaSize = params.covEstimationAreaSize;
   end
   if  ~isfield(params,'covFactorization') || isempty(params.covFactorization)
      params.covFactorization= 'Cholesky';
   else
      params.covFactorization = params.covFactorization;
   end
else
   params.correctionType = 'none';
end

%do not do anything if nothing is asked
if isempty(restrictions) && isempty(contrasts) && ~computeEstimates && ~computeTtests
  return;
end


%--------------------------------------------------- INITIALIZATION ------------------------------------%
% precalculate the normal equation for Ordinary Least Squares
covEVs = d.scm'*d.scm; %covariance matrix of the EVs (do I need to correct for covariance ?)
covarianceMatrixRank = rank(covEVs);
% use normal equations if we have a full ranked covariance matrix
if covarianceMatrixRank == size(covEVs,1)
   invCovEVs = covEVs^-1; % 
   pinv_X = covEVs\d.scm'; 
        
% otherwise use pinv
else
  % note that if we need to use the pseudo inverse it means that there is ambiguity in the design
  % such that there are an infinite number of possible solutions. The pseudo-inverse solution
  % choses the solution with the minimum length (i.e. Euclidian norm)
  if verbose, mrWarnDlg(sprintf('(getGlmStatistics) Design covariance matrix (%ix%i) is rank %i. Using pseudo-inverse to invert.',size(covEVs,1),size(covEVs,2),covarianceMatrixRank));end
  pinv_X = pinv(d.scm);
  invCovEVs = pinv(d.scm'*d.scm);
end


nFrames = size(d.volumes,2);

residualForming = eye(nFrames) - d.scm*pinv_X; %this matrix computes the residuals directly
d.rdf = nFrames-size(d.scm,2)-1; %degrees of freedom for the residual          %SHOULD BE nFrames-size(d.scm,2)

%initialize variables for covariance matrix estimation
if params.covCorrection 
   if params.covEstimationAreaSize>1
      sliceAverager3D = ones(1,params.covEstimationAreaSize,params.covEstimationAreaSize)/params.covEstimationAreaSize^2;
      if ~isfield(d,'roiPositionInBox') %this is in case we only compute for voxels in the loaded ROIs
         longMargin = ceil(params.covEstimationAreaSize/2); %we exclude voxels that are not surrounded by voxels with valid data
         shortMargin = floor(params.covEstimationAreaSize/2);
      end
   end
   switch params.covEstimation
      case 'singleTukeyTapers'
         tukeyM = 10;
         tukeyM = min(tukeyM,d.dim(4)-1);
   end
end

if strcmp(params.componentsCombination,'Or') %in this case, components of a given EV are tested independently from each other in an F-test
  params.componentsToTest = logical(diag(params.componentsToTest));
  %for contrasts, this amounts to  testing several contrats at once using an f test
  %So we have to change contrasts into f-tests (with the restriction that they have to be two-sided; this should be controlled for by the parameters)
  restrictions = [contrasts restrictions];
  contrasts = {};
end

%expand restriction matrices using kronecker products
for iR = 1:length(restrictions)
  restrictions{iR} = kron(restrictions{iR},params.componentsToTest);
end
for iContrast = 1:length(contrasts)
  contrasts{iContrast} = kron(contrasts{iContrast},params.componentsToTest);
end


if ~isempty(restrictions)
  
  %this has not been tested, but this is only for generalized F-tests which do not give the expected results anyway
  if strcmp(params.correctionType,'generalizedFTest')
    complementaryRestriction = cell(1,size(restrictions,1));
    baseRestriction =  kron(logical(eye(d.nhdr)),logical(params.componentsToTest)); 
    for iR = 1:length(restrictions)
      % not sure at all about these lines
      complementaryRestriction{iR} = baseRestriction - restrictions{iR};
      complementaryRestriction{iR} = complementaryRestriction{iR}(any(complementaryRestriction{iR},2),:);  
      % it actually doesn't make sense anymore now that F-tests can be made of any contrast
      % The logic relied on the fact that there were only elements on the diagonal of the restrictions matrix
      % There must be a generalization but I don't have time for this now
    end
    residualForming_f_tests = NaN(d.dim(4),d.dim(4),length(restrictions));
    residualForming_h = NaN(d.dim(4),d.dim(4),length(restrictions));
    traceRsV = NaN([d.dim(1:3) nBootstrap length(restrictions)]);
    effective_rdf = NaN(d.dim(1:3));
    effective_mdf = NaN([d.dim(1:3) nBootstrap length(restrictions)]);
    if ~approximate_value
      traceRV = NaN(d.dim(1:3));
    end
    scm_h = cell(1:length(restrictions));
    eig_residualForming_f_tests = cell(1:length(restrictions));
    for iR = 1:length(restrictions)
      scm_h{iR} = d.scm*complementaryRestriction{iR}';           
      residualForming_h(:,:,iR) = eye(nFrames) - scm_h{iR}*((scm_h{iR}'*scm_h{iR})^-1)*scm_h{iR}';       %TO REMOVE EVENTUALLY
      residualForming_f_tests(:,:,iR) = residualForming_h(:,:,iR) - residualForming; 
      [V,dump] = eig(residualForming_f_tests(:,:,iR));
      eig_residualForming_f_tests{iR} = V(:,1:size(restrictions{iR},1));
    end
  end
  
  mss = NaN(d.dim(1),d.dim(2),d.dim(3),nBootstrap,length(restrictions),precision);
  R_invCovEV_Rp = cell(1:length(restrictions));
  for iR = 1:length(restrictions)         %the degrees of freedom for the F-tests are the number of contrasts
    restrictions{iR} = restrictions{iR}(any(restrictions{iR},2),:); %remove lines of 0
    d.mdf(iR) = size(restrictions{iR},1); %this is provided that contrasts are independent and that they are no more than (regressors -1)
    %inv_R_invCovEV_Rp{iR} = (restrictions{iR}*invCovEVs*restrictions{iR}')^-1; 
    %I replaced all the precomputed inverses (previous line) by ml/mrDivide with the non-inverted matrix (faster and more accurate):
    R_invCovEV_Rp{iR} = restrictions{iR}*invCovEVs*restrictions{iR}';
  end
                         
end

% check roi
xvals = 1:d.dim(1);
yvals = 1:d.dim(2);
slices = 1:d.dim(3);
  


%--------------------------------------------------- PRE-ALLOCATION ------------------------------------%
ehdr = NaN(d.nhdr*d.nHrfComponents,d.dim(1),d.dim(2),d.dim(3),nBootstrap,precision);    %JB: replaces: d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);

if computeEstimates ||computeTtests
  ehdrste = NaN(d.nhdr*d.nHrfComponents,d.dim(1),d.dim(2),d.dim(3),precision); %JB: replaces: d.ehdrste = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
  contrastBetaSte = NaN([d.dim(1) d.dim(2) d.dim(3) nBootstrap params.numberContrasts],precision);
end
if ~isempty(restrictions) || computeTtests
  s2 = NaN(d.dim(1),d.dim(2),d.dim(3),precision); %residual variance
  rss = NaN(d.dim(1),d.dim(2),d.dim(3),precision); %JB: this is to store the sum of square of the residual error term
  tss = NaN(d.dim(1),d.dim(2),d.dim(3),precision); %JB: this is to store the total sum of square
end

% turn off warnings to avoid divide by zero warning
warning('off','MATLAB:divideByZero');



%--------------------------------------------------- MAIN LOOP: PARAMETER ESTIMATION ------------------------------------%
% display string
if verbose,hWaitBar = mrWaitBar(-1/(d.dim(2)*d.dim(3)),'(getGlmStatistics) Estimating model parameters');end
% cycle through images calculating the estimated hdr and r^2s of the estimate.
for z = slices
  residuals = NaN(d.dim(4),d.dim(1),d.dim(2),precision);
  %permute timeseries frames in first dimension to compute residuals
  timeseries = permute(d.data(:,:,z,d.volumes),[4 1 2 3]); 
  % subtract off column means
  colmeans = repmat(mean(timeseries,1),[nFrames 1 1]);
  timeseries = timeseries - colmeans;
  % convert to percent signal change
  timeseries = 100*timeseries./colmeans;
  clear('colmeans');

  %--------------------------------------------------- MAIN LOOP: ORDINARY LEAST SQUARES ------------------------------------%
  % the following section has been optimized to run faster by
  % eliminating the x loops. in addition, the x dimension is divided in chunks of 60 points.
  % (which makes a difference only when system swaps ?)
  for y = yvals
    % get OLS residuals
    for iX = 1:ceil(d.dim(1)/optimalXLength)
       xValues = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
       residuals(:,xValues,y) = residualForming*timeseries(:,xValues,y); %JB resplaces residuals = timeseries-d.scm*d.ehdr(:,:,y,z);
    end
    %residuals(:,:,y) = residualForming*timeseries(:,:,y); %This is too long if X dim is large (when calling with all data on the first dimension)
  end
  if computeEstimates || computeTtests || ~isempty(restrictions)
   %put x and y back in the 1st and 2nd dimensions
   rss(:,:,z) = permute(sum(residuals.^2),[2 3 1]);    %JB: replaces: sumOfSquaresResidual = sum((timeseries-d.scm*ehdr{y,z}).^2);
  end

  

  if params.covCorrection
    if params.covEstimationAreaSize>1     
      %average the residuals for covariance matrix estimation
      if isfield(d,'roiPositionInBox') %if the data are not spatially organized, we need to temporarily put them in a volume
        %first put them on a single line whose length is the product of the volume dimensions
        roiPositionInBox = reshape(d.roiPositionInBox,numel(d.roiPositionInBox),1);
        averaged_residuals = NaN(d.dim(4),size(roiPositionInBox,1),precision); 
        averaged_residuals(:,roiPositionInBox>0) = residuals;  
        %then reshape into a volume
        averaged_residuals = reshape(averaged_residuals,[d.dim(4) size(d.roiPositionInBox)]);
        %average
        averaged_residuals = convn(averaged_residuals,sliceAverager3D,'same');
        %put the data back in the original array size
        averaged_residuals = reshape(averaged_residuals,d.dim(4), numel(d.roiPositionInBox));
        averaged_residuals = averaged_residuals(:,d.roiPositionInBox>0);
        clear('roiPositionInBox');
      else
        averaged_residuals = NaN(size(residuals),precision);
        %averaged_residuals(isnan(residuals)) = 0;            %convert NaNs to 0... actually no, don't
        %average residuals over space
        averaged_residuals(:,longMargin:end-shortMargin,longMargin:end-shortMargin) = convn(residuals,sliceAverager3D,'valid');
      end
    else
      averaged_residuals = residuals;
    end
  end

  for y = yvals  %for F, noise covariance correction and residual bootstrapping we need to loop over y (and x)
      
%--------------------------------------------------- MAIN LOOP: NOISE COVARIANCE CORRECTION ------------------------------------%
      %if verbose,mrWaitBar( ((y-min(yvals)) + d.dim(2)*(z-min(slices))) / (d.dim(2)*d.dim(3)),hWaitBar);end
      
      % Covariance matrix estimation 
      if params.covCorrection
         autoCorrelation = [ones(1,d.dim(1),'double') ; zeros(d.dim(4)-1,d.dim(1),'double')]; %you don't want to make those single
         %because inversion operations are faster on doubles
         switch(params.covEstimation)
            %Single tapers
            case 'singleTukeyTapers'      %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
               this_s2 = var(averaged_residuals(:,:,y),1,1);   
               for tau=1:tukeyM-1   %NB as long as assignments are made in loops, the autocorrelation matrix stays in double format even if the data is single
                  autoCorrelation(tau+1,:) = .5*(1+cos(pi*tau/tukeyM)) * (sum(averaged_residuals(1:end-tau,:,y).*averaged_residuals(tau+1:end,:,y),1))./this_s2/(d.dim(4)-tau);
               end

            case 'nonParametricResiduals'
                %not implemented
            case 'nonParametricTimeSeries' %Can we do this on the timeseries instead of the residuals, like Wicker et al ?
                %not implemented
            case 'dampenedOscillator' %see Kruggel et al. (2002) Medical image Analysis, 6, p65
               %check if invertible/definite positive (might have negative values ?)
               this_s2 = var(averaged_residuals(:,:,y),1,1);  
               M = d.dim(4)-1;
               %M=30;
               for tau=1:M
                  autoCorrelation(tau+1,:) = (sum(averaged_residuals(1:end-tau,:,y).*averaged_residuals(tau+1:end,:,y),1))./this_s2/(d.dim(4)-tau);
               end
               %initialParameters = [0 0 0];
               initialParameters = rand(1,3);
               options = optimset('Display','iter');
               for x = xvals
                  [dampenedOscillatorParams(x,:), dummy, exitflag] = fminsearch(@(dummy) minimizeDampenedOscillator(dummy,(1:M)',autoCorrelation(1+(1:M),x)), initialParameters,options);
               end
               %not finished ...
         end
      end
      
      switch(params.correctionType)
         case 'none'
            if verbose,mrWaitBar( ((y-min(yvals)) + d.dim(2)*(z-min(slices)))/ prod(d.dim(2:3)), hWaitBar);end
            % get OLS estimates
            ehdr(:,:,y,z) = pinv_X*timeseries(:,:,y);    %JB: replaces: ehdr{y,z} = pinv_X*timeseries;
            % compute model sum of squares for F-tests
            %from Burock and Dale 2000 
            for iR = 1:length(restrictions)
               ss_beta = restrictions{iR}*ehdr(:,:,y,z,1);
               %computing time seems to be the fastest on my machine when computing 60 values at a time
               %but this might be because it was swapping. anyway, I'll leave it this way because if it doesn't swap
               %it doesn't seem to make much difference
%                optimalXLength = 10:200;                             %DEBUG/OPTIMIZATION
%                computingTimes = zeros(size(optimalXLength));        %DEBUG/OPTIMIZATION
%                for j = 1:length(optimalXLength)                     %DEBUG/OPTIMIZATION
%                   tic
                  for iX = 1:ceil(d.dim(1)/optimalXLength)
                     xValues = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
                     %mss(xValues,y,z,1,iR) = diag( ss_beta(:,xValues)' * inv_R_invCovEV_Rp{iR} * ss_beta(:,xValues) ); 
                     mss(xValues,y,z,1,iR) = diag( ss_beta(:,xValues)' / R_invCovEV_Rp{iR} * ss_beta(:,xValues) );  
                  end
%                   computingTimes(j) = toc;                          %DEBUG/OPTIMIZATION
%                end                                                  %DEBUG/OPTIMIZATION
%                figure;plot(optimalXLength,computingTimes);          %DEBUG/OPTIMIZATION

               %mss(:,y,z,1,iR) = diag( ss_beta' * inv_R_invCovEV_Rp{iR} * ss_beta ); %%%%%this is too slow when Xdim is large
               %mss(:,y,z,1,iR) = diag( ss_beta' / R_invCovEV_Rp{iR} * ss_beta ); 
               
            end
            

         case 'generalizedLeastSquares' %see Wicker & Fonlupt (2003) NeuroImage, 18, p589 and Burock & Dale (2000) Human Brain Mapping, 11, p249
            for x = xvals
               if verbose,mrWaitBar( (((y-min(yvals)) + d.dim(2)*(z-min(slices)))*d.dim(1)+x-min(xvals))/ prod(d.dim(1:3)), hWaitBar);end
               if ~any(isnan(timeseries(:,x,y))) && ~any(isnan(autoCorrelation(:,x)))
                  covResiduals = toeplitz(autoCorrelation(:,x));
                  
                  correctedCovEV = d.scm' / covResiduals * d.scm;
                  corrected_pinv_X = correctedCovEV \ d.scm';
                  ehdr(:,x,y,z) = corrected_pinv_X  / covResiduals * timeseries(:,x,y);
                  if computeEstimates || computeTtests || ~isempty(restrictions) || nBootstrap>1
                    residuals(:,x,y) = timeseries(:,x,y) - d.scm*ehdr(:,x,y,z);
                    if computeEstimates || computeTtests || ~isempty(restrictions)
                      s2(x,y,z) = residuals(:,x,y)' / covResiduals * residuals(:,x,y)/d.rdf;    %Wicker & Fonlupt (2003) 
                      %or
                      %preFilter = chol(invCovResiduals);                        %Burock & Dale (2000)
                      %s2(x,y,z) = sum((preFilter*residuals(:,x,y)).^2,1)/d.rdf;     %Slightly slower and give identical results
                      for iR = 1:length(restrictions)
                         corrected_R_invCovEV_Rp = restrictions{iR} / correctedCovEV * restrictions{iR}';
                         ss_beta = restrictions{iR}*ehdr(:,x,y,z,1);
                         mss(x,y,z,1,iR) = ss_beta' / corrected_R_invCovEV_Rp * ss_beta; 
                      end

    %                   %Wicker & Fonlupt (2003) about twice slower (and not generalizable to contrasts ?) - although didn't try replacing inv by mldivide/mrdivide
    %                   A = d.scm * corrected_pinv_X  * invCovResiduals;
    %                   residuals(:,x,y) = (eye(d.dim(4)) - A)*timeseries(:,x,y);
    %                   for iR = 1:length(restrictions)                                                                        
    %                      Ar = scm_h{iR} * (scm_h{iR}'* invCovResiduals * scm_h{iR})^-1 * scm_h{iR}' * invCovResiduals;
    %                      mss(x,y,z,1,iR) = timeseries(:,x,y)' * (A - Ar)' * invCovResiduals * (A - Ar) * timeseries(:,x,y);  
    %                   end                                                                                                   

                      %compute the corrected inverse covariance matrix
                      if computeEstimates || (isempty(contrasts) && computeTtests)
                        invCorrectedCovEV = inv(correctedCovEV);
                        % compute the std deviation of the beta estimates
                        if computeEstimates
                          ehdrste(:,x,y,z) = sqrt(s2(x,y,z) * diag(invCorrectedCovEV));        
                          % see for example Seber & Lee, Linear Regression Analysis, 1977, p67
                        end
                        if computeTtests
                        %this is from Woolrich et al. (2001) and extends the method to T-tests and therefore allows one-sided tests
                          for iContrast = 1:length(contrasts)
                            %k_eff(x,y,z,nBootstrap,iContrast) = contrasts{iContrast} / correctedCovEV * contrasts{iContrast}';
                            contrastBetaSte(x,y,z,nBootstrap,iContrast) = sqrt(s2(x,y,z)*contrasts{iContrast} * invCorrectedCovEV * contrasts{iContrast}');
                          end
                        end
                      end
                    end
                  end
               end
            end
            
         case 'varianceCorrection' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
            ehdr(:,:,y,z) = pinv_X*timeseries(:,:,y);    % get OLS hdr
            if computeTtests || ~isempty(restrictions) || computeEstimates
              for x = xvals
                if verbose,mrWaitBar( (((y-min(yvals)) + d.dim(2)*(z-min(slices)))*d.dim(1)+x-min(xvals))/ prod(d.dim(1:3)), hWaitBar);end
                if ~any(isnan(timeseries(:,x,y))) && ~any(isnan(autoCorrelation(:,x)))
                  covResiduals = (toeplitz(autoCorrelation(:,x)));
                  s2(x,y,z,1) = rss(x,y,z)/trace(residualForming * covResiduals);
                  corrected_invCovEV = pinv_X * covResiduals * pinv_X';
                  % compute the std error (variance of beta estimates)
                  ehdrste(:,x,y,z) = sqrt(s2(x,y,z) * diag(corrected_invCovEV));        
                  % see for example Seber & Lee, Linear Regression Analysis, 1977, p67
                  if computeTtests || ~isempty(restrictions)
                    if ~isempty(contrasts)
                       for iContrast = 1:length(contrasts)
                         %k_eff(x,y,z,nBootstrap,iContrast) = contrasts{iContrast} * corrected_invCovEV * contrasts{iContrast}';
                         contrastBetaSte(x,y,z,nBootstrap,iContrast) = sqrt(s2(x,y,z)*contrasts{iContrast} * corrected_invCovEV * contrasts{iContrast}');
                       end
                    end
                    %this is from Burock & Dale and extends the method to F-tests (and F-tests extended to contrasts)
                    for iR = 1:length(restrictions)
                       corrected_R_invCovEV_Rp = restrictions{iR}*corrected_invCovEV*restrictions{iR}';
                       ss_beta = restrictions{iR}*ehdr(:,x,y,z,1);
                       mss(x,y,z,1,iR) = ss_beta' / corrected_R_invCovEV_Rp * ss_beta;
                    end
                  end
                end
              end
            end
            

         case 'preWhitening' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
            for x = xvals
               if verbose,mrWaitBar( (((y-min(yvals)) + d.dim(2)*(z-min(slices)))*d.dim(1)+x-min(xvals))/ prod(d.dim(1:3)), hWaitBar);end
               if ~any(isnan(timeseries(:,x,y))) && ~any(isnan(autoCorrelation(:,x)))
                  tic;
                  covResiduals = toeplitz(autoCorrelation(:,x));
                  
                  invCovResiduals = inv(covResiduals);
                  % if this doesn't work, do pinv
                  if sum(isnan(invCovResiduals(:))) == length(invCovResiduals(:))
                     invCovResiduals = pinv(toeplitz(autoCorrelation(:,x)));
                     oneTimeWarning('pseudoInverseWarning','Covariance matrix is nearly singular, using pinv');
                  end
                  %compute inv(Cov)^1/2
                  switch(params.covFactorization)
                     case 'Cholesky'
                       try
                        preFilter = chol(invCovResiduals);
                       catch exception
                         mrDisp(sprintf('(geGlmStatistics) Cannot factorize inverse noise covariance matrix because it is not positive definite (voxel %d %d %d)\n',x,y,z));
                         preFilter = [];
                       end
                  end
                  if ~isempty(preFilter)
                    white_scm = preFilter * d.scm;                        %%%% IS IT POSSIBLE TO SKIP THE FACTORIZATION BY ONLY USING (COV)^-1 LATER ON ??
                    timeseries(:,x,y) = preFilter * timeseries(:,x,y);    %%%% NOT ACCORDING TO Kruggel et al. 2002...
                    white_pinv_X = (white_scm'*white_scm)\white_scm';
                    white_residualForming = eye(nFrames) - white_scm*white_pinv_X;
                    ehdr(:,x,y,z) = white_pinv_X * timeseries(:,x,y);
                    if computeEstimates || computeTtests || ~isempty(restrictions) || nBootstrap>1
                      residuals(:,x,y) = white_residualForming * timeseries(:,x,y);
                      if computeEstimates || computeTtests || ~isempty(restrictions)
                        rss(x,y,z) = sum(residuals(:,x,y).^2)'; 
                        s2(x,y,z) = rss(x,y,z)/trace(white_residualForming * preFilter * covResiduals * preFilter');
                        corrected_invCovEV = inv(d.scm' * invCovResiduals * d.scm);
                        if computeEstimates
                          % compute the std deviation of the beta estimates
                          ehdrste(:,x,y,z) = sqrt(s2(x,y,z) * diag(corrected_invCovEV));        
                          % see for example Seber & Lee, Linear Regression Analysis, 1977, p67
                        end
                        if computeTtests || ~isemtpy(restrictions)
                          if ~isempty(contrasts) && computeTtests
                            for iContrast = 1:length(contrasts)
                              %k_eff(x,y,z,nBootstrap,iContrast) = contrasts{iContrast} * corrected_invCovEV * contrasts{iContrast}';
                              contrastBetaSte(x,y,z,nBootstrap,iContrast) = sqrt(s2(x,y,z)*contrasts{iContrast} * corrected_invCovEV * contrasts{iContrast}');
                            end
                          end

                        %this is from Burock & Dale and extends method to F-tests (and F-tests extended to contrasts)
                        %corrected_invCovEV is equivalent to (white_scm'*white_scm)^-1 
                              %[this is because prefilter is chosen such that prefilter'*prefilter = invCovResiduals
                              % and then corrected_invCovEV = (d.scm'* invCovResiduals * d.scm)^-1 = (d.scm'* preFilter' * preFilter * d.scm)^-1 = (white_scm'* white_scm)^-1 ]
                          for iR = 1:length(restrictions)
                             corrected_inv_R_invCovEV_Rp = (restrictions{iR}*corrected_invCovEV*restrictions{iR}')^-1;
                             ss_beta = restrictions{iR}*ehdr(:,x,y,z,1);
                             mss(x,y,z,1,iR) = ss_beta' * corrected_inv_R_invCovEV_Rp * ss_beta;
                          end
                        end
                      end
                    end
                   end
               end
            end

         case 'generalizedFTest' %see Kruggel et al. (2002) Medical image Analysis, 6, p65
           %Does not give the expected result, but didn't have time to debug
           %furthermore, less accurate and not obvious that it is much faster than other methods, so not worth the trouble
            ehdr(:,:,y,z) = pinv_X*timeseries(:,:,y);    % get OLS hdr
            if computeTtests || ~isempty(restrictions) || computeEstimates
              for x = xvals
                if verbose,mrWaitBar( (((y-min(yvals)) + d.dim(2)*(z-min(slices)))*d.dim(1)+x-min(xvals))/ prod(d.dim(1:3)), hWaitBar);end
                if ~any(isnan(timeseries(:,x,y))) && ~any(isnan(autoCorrelation(:,x)))
                 covResiduals = toeplitz(autoCorrelation(:,x));
                  if approximate_value 
                     effective_rdf(x,y,z) = d.dim(4)*(d.dim(4)-d.nhdr*d.nHrfComponents)/trace(covResiduals*covResiduals);
                  else %this is the exact value (longer and not very different for large N, which is usually the case for fMRI data)
                     RV = residualForming*covResiduals;           
                     traceRV(x,y,z) = trace(RV);       
                     effective_rdf(x,y,z) = traceRV(x,y,z)^2/trace(RV*RV);  
                  end

                  for iR = 1:length(restrictions)
                     if approximate_value
                        temp = eig_residualForming_f_tests{iR}'*covResiduals*eig_residualForming_f_tests{iR};   %SEE IF FASTER WHEN LOOP INSTEAD OF TRACE
                        effective_mdf(x,y,z,1,iR) = trace(temp)^2 / trace(temp.^2);                                      %SEE IF FASTER WHEN LOOP INSTEAD OF TRACE
                        traceRsV(x,y,z,1,iR) = trace(residualForming_f_tests(:,:,iR)*covResiduals);
                     else
                        RsV = residualForming_f_tests(:,:,iR)*covResiduals;  %%this is the exact value (longer)
                        traceRsV(x,y,z,1,iR) = trace(RsV);
                        effective_mdf(x,y,z,1,iR) = traceRsV(x,y,z,1,iR)^2/trace(RsV*RsV);    
                     end
                  end
                end
              end
            end
            %from Burock and Dale 2000 
            for iR = 1:length(restrictions)
               ss_beta = restrictions{iR}*ehdr(:,:,y,z,1);
                  for iX = 1:ceil(d.dim(1)/optimalXLength)
                     xValues = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
                     %mss(xValues,y,z,1,iR) = diag( ss_beta(:,xValues)' * inv_R_invCovEV_Rp{iR} * ss_beta(:,xValues) ); 
                     mss(xValues,y,z,1,iR) = diag( ss_beta(:,xValues)' / R_invCovEV_Rp{iR} * ss_beta(:,xValues) ); 
                  end
%                mss(:,y,z,1,iR) = diag( ss_beta' * inv_R_invCovEV_Rp{iR} * ss_beta );%%%%%%%%%%%%%%%%%%%%%%%%%%THIS IS LIKELY TO BE SLOW WHEN X IS LARGE
                %mss(:,y,z,1,iR) = diag( ss_beta' / R_invCovEV_Rp{iR} * ss_beta ); 
            end
      end

 
%--------------------------------------------------- MAIN LOOP: BOOSTRAP CONFIDENCE INTERVAL ESTIMATION ------------------------------------%
      if nBootstrap>1 %BOOTSTRAP ONLY IMPLEMENTED FOR ORDINARY LEAST SQUARES
         
         %recompute the parameter estimates after bootstrapping the residuals to construct bootstrap timeseries
         if isempty(contrasts) && isempty(restrictions)
            bootstrap_ehdr=bootstrp(nBootstrap-1,@(bootstrap_residuals)pinv_X*(d.scm*ehdr(:,:,y,z)+bootstrap_residuals),residuals(:,:,y));
            ehdr(:,:,y,z,2:nBootstrap) = permute(reshape(bootstrap_ehdr,nBootstrap-1,size(ehdr,1),size(ehdr,2)),[2 3 4 5 1]);
         else
         %if we need to bootstrap estimate contrasts or fTests, we need more info, but it's slower
            bootstrap_index = ceil(size(residuals,1)*rand(nBootstrap-1,size(residuals,1)));
            %h_wait = waitbar(0,['Bootstrapping residuals for y=' num2str(y) ' and z=' num2str(z)]);
            for i = 1:ceil(size(residuals,2)/bootstrap_max_voxels)
               %waitbar(i/ceil(size(residuals,2)/bootstrap_max_voxels),h_wait);
               x_subset = (i-1)*bootstrap_max_voxels+1:min(i*bootstrap_max_voxels,size(residuals,2)); 
               res_temp = residuals(:,x_subset,y);
               for i_boot = 1:nBootstrap-1 
                  bootstrap_res_temp = res_temp(bootstrap_index(i_boot,:),:);
                  bootstrap_ts = d.scm*ehdr(:,x_subset,y,z) + bootstrap_res_temp;
                  ehdr(:,x_subset,y,z,i_boot+1) = pinv_X*bootstrap_ts; 
                  s2(x_subset,y,z,i_boot+1) = squeeze(sum((bootstrap_ts-d.scm*ehdr(:,x_subset,y,z,i_boot+1)).^2,1)/d.rdf);
                  for iR = 1:length(restrictions)
                     ss_beta = restrictions{iR}*ehdr(:,x_subset,y,z,i_boot+1);
                     %mss(x_subset,y,z,i_boot+1,iR) = diag( ss_beta' * inv_R_invCovEV_Rp{iR} * ss_beta );
                     mss(x_subset,y,z,i_boot+1,iR) = diag( ss_beta' / R_invCovEV_Rp{iR} * ss_beta ); 
                 end
               end
            end
            %close(h_wait);
         end
      end
   end
   clear('averaged_residuals','residuals');
   if computeEstimates
    tss(:,:,z) = permute(sum(timeseries.^2,1),[2 3 1]);    %JB 
   end
   clear('timeseries');
end
% disp([num2str(bootstrap_max_voxels) ' ' num2str(t) ';']);   %DEBUG

clear('residualForming');
oneTimeWarning('pseudoInverseWarning',0);

%--------------------------------------------------- STATISTICS ------------------------------------%


%---------------------------------- VARIANCE OF BETAS (STANDARD ERROR)------------------------------%
if strcmp(params.correctionType,'none')
  %compute residual variance for OLS
  s2(:,:,:,1) = rss/d.rdf;   %JB: replaces: S2 = sumOfSquaresResidual/(length(d.volumes)-size(d.scm,2));
  
  if computeEstimates
    % now distribute the error to each one of the points    
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    % (= compute the variance of the beta estimates, i.e. the standard error)
    ehdrste = sqrt(repmat(diag(invCovEVs),[1 d.dim(1:3) nBootstrap]).*repmat(permute(s2,[5 1 2 3 4]),[d.nhdr*d.nHrfComponents 1 1 1 1]));        %JB: replaces: ehdrste{y,z} = sqrt(diag(pinv(d.scm'*d.scm))*S2);
    % see for example Seber & Lee, Linear Regression Analysis, 1977, p42
  end
end


%---------------------------------- F-TESTS --------------------------------------------------------%
% basically, the F value (only for the EVs of interest) is 
%     ((MSS)/ number of betas of interest) / (RSS/(number of samples - number of betas ?-1?))
%     where MSS is the difference between the RSS without the EVs of interest and the RSS of the whole model
% MSS has been computed as betas'*R'*(R*(X'*X)^-1*R')^-1 * R * betas, where R is a restrictions matrix that isolates the EVs of interest (such that R*X = Xinterest)
% computed this way, the f-test can also be extended to contrasts and is equivalent to a two-sided T-test, if R describes a linear combination of EVs
if ~isempty(restrictions)
   
   switch(params.correctionType)
      case {'none','varianceCorrection','generalizedLeastSquares','preWhitening'}
         F = NaN(size(mss),precision);
         for iR = 1:length(restrictions)
            F(:,:,:,:,iR) = (mss(:,:,:,:,iR)/d.mdf(iR)) ./ s2;
         end
      case 'generalizedFTest'
         if approximate_value
            F = repmat(effective_rdf.^2,[1 1 1 1 length(restrictions)]) .* traceRsV .* mss ./ effective_mdf.^2 ./ repmat(rss,[1 1 1 1 length(restrictions)]) / (d.rdf+1);
         else
            F = repmat(effective_rdf.^2,[1 1 1 1 length(restrictions)]) .* traceRsV .* mss ./ effective_mdf.^2 ./ repmat(rss,[1 1 1 1 length(restrictions)]) ./ repmat(traceRV,[1 1 1 1 length(restrictions)]);
         end
   end

   F = permute(F,[1 2 3 5 4]); %put nBootstrap dim at the end
end
%put EV dimensions last so that contrasts are easier to perform
ehdr = permute(ehdr,[2 3 4 5 1]); 

%--------------------------------------------------- T-TESTS ------------------------------------%
if ~isempty(contrasts)
  contrastBetas = zeros([d.dim(1) d.dim(2) d.dim(3) nBootstrap length(contrasts)],precision);
  if strcmp(params.correctionType,'none')
    contrastBetaSte = NaN([d.dim(1) d.dim(2) d.dim(3) nBootstrap length(contrasts)],precision);
    for iContrast = 1:length(contrasts)
      %k_eff(:,:,:,:,iContrast) = repmat(contrasts{iContrast}*invCovEVs*contrasts{iContrast}',[d.dim(1) d.dim(2) d.dim(3) nBootstrap]);
      contrastBetaSte(:,:,:,:,iContrast) = sqrt(s2.*(contrasts{iContrast}*invCovEVs*contrasts{iContrast}'));
    end
  end
  for iContrast = 1:length(contrasts)
    for i_ev=1:d.nhdr*d.nHrfComponents
       contrastBetas(:,:,:,:,iContrast) = contrastBetas(:,:,:,:,iContrast)+contrasts{iContrast}(i_ev)*ehdr(:,:,:,:,i_ev);
    end
  end
  if computeTtests
    T = contrastBetas ./ contrastBetaSte;   
%     T = zeros([d.dim(1) d.dim(2) d.dim(3) nBootstrap length(contrasts)],precision);
%     for iContrast = 1:length(contrasts)
%       T(:,:,:,:,iContrast) = contrastBetas(:,:,:,:,iContrast) ./ (k_eff(:,:,:,:,iContrast).*s2).^(1/2);   
%     end
    T = permute(T,[1 2 3 5 4]); %put nBootstrap dim at the end
    switch(params.tTestSide)
     case 'Both'
        T = abs(T);
     case 'Left'
        T = -1 *T;
    end
  end
  contrastBetas = permute(contrastBetas,[1 2 3 5 4]); %put nBootstrap dim at the end
  contrastBetaSte = permute(contrastBetaSte,[1 2 3 5 4]); %put nBootstrap dim at the end
end

%Now in the case we computed contrasts on several components using option 'Or',
%we have to convert the appropriate F values into T values
if strcmp(params.componentsCombination,'Or') && ~isempty(params.contrasts)
  %T values are the square roots of the numberContrasts first F values 
  T = sqrt(F(:,:,:,1:params.numberContrasts,:));
  %remove T values form F values array
  F = F(:,:,:,params.numberContrasts+1:end,:);
  d.mdf=d.mdf(params.numberContrasts+1:end);
end

%reshape to the right dimensions
if computeEstimates
  % calculate variance accounted for by the estimated hdr
  d.r2 = 1 - rss./tss;     %JB: replaces: r2{y,z} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  ehdrste = permute(reshape(ehdrste, d.nHrfComponents, d.nhdr, d.dim(1),d.dim(2),d.dim(3)),[3 4 5 2 1]);  %JB
  d.ehdrste = ehdrste;
  %reshape nhdr and nhdrlen on two dimensions, put bootstrap dimension last, and swap nHrfComponents and nhdr order
  ehdr = permute(reshape(ehdr, d.dim(1), d.dim(2), d.dim(3), nBootstrap, d.nHrfComponents, d.nhdr),[1 2 3 6 5 4]);    
  d.ehdr = ehdr;
  d.rss = rss;
  d.tss = tss;
  d.contrastSte = contrastBetaSte;
end

if verbose,mrCloseDlg(hWaitBar);end



%this function computes a new model from scm, ehdr and residuals 
%and then fits the model described by scm amd pinv_X 
%to compute new estimates and residuals, that are concatenated in a single output matrix
%so that it can be used by bootstrp (not used finally)
function new_residuals_ehdr = bootstrapResiduals(pinv_X, scm, ehdr, residuals)

bootstrap_ts = scm*ehdr+residuals;
ehdr = pinv_X*bootstrap_ts; 
residuals = bootstrap_ts-scm*ehdr;
new_residuals_ehdr = cat(1,residuals,ehdr);

%this function computes the sum of squared errors between the dampened oscillator
%model (for xdata) and the sample autocorrelation function (ydata)
function sse = minimizeDampenedOscillator(params, xdata,ydata)
  FittedCurve = params(1)^2 - exp(params(2) * xdata) .* cos(params(3)*xdata);
  ErrorVector = FittedCurve - ydata;
  sse = sum(ErrorVector.^2);
