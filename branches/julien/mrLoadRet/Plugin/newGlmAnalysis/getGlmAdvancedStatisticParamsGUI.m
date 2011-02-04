% getGlmAdvancedStatisticParamsGUI.m
%
%        $Id: getGlmAdvancedStatisticParamsGUI.m 2013 2011-01-21 17:43:19Z julien $
%      usage: params = getGlmAdvancedStatisticParamsGUI(params,useDefault)
%         by: julien besle, 
%       date: 19/01/2011
%    purpose: returns additional statistical test parameters for GLM analysis
%

function params = getGlmAdvancedStatisticParamsGUI(thisView,params,useDefault)

keepAsking = 1;

while keepAsking
  %Default params
  if fieldIsNotDefined(params, 'testOutput')
      params.testOutput = 'Z value';
  end
  if fieldIsNotDefined(params,'outputStatistic')
    params.outputStatistic = 0;
  end

  if fieldIsNotDefined(params, 'bootstrapIntervals')
      params.bootstrapIntervals = 0;
  end
  if fieldIsNotDefined(params, 'alphaConfidenceIntervals')
      params.alphaConfidenceIntervals = 0.05;
  end
  
  if fieldIsNotDefined(params, 'resampleFweMethod')
    params.resampleFweMethod = 'Adaptive Step-down';
  end
  
  if fieldIsNotDefined(params, 'fweAdjustment')
      params.fweAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'fweMethod')
      params.fweMethod = 'Adaptive Step-down';
  end
  
  if fieldIsNotDefined(params, 'fdrAdjustment')
      params.fdrAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'fdrAssumption')
      params.fdrAssumption = 'Independence/Positive dependence';
  end
  if fieldIsNotDefined(params, 'fdrMethod')
      params.fdrMethod = 'Adaptive Step-up';
  end
  
  if fieldIsNotDefined(params, 'trueNullsEstimationMethod')
      params.trueNullsEstimationMethod = 'Least Squares';
  end
  if fieldIsNotDefined(params, 'trueNullsEstimationThreshold')
      params.trueNullsEstimationThreshold = .05;
  end
  if fieldIsNotDefined(params, 'maskContrastOverlay')
      params.maskContrastOverlay = 1;
  end
  if fieldIsNotDefined(params, 'noDoubleBootstrap')
      params.noDoubleBootstrap = 0;
  end
  
  
  testOutputMenu = putOnTopOfList(params.testOutput,{'P value','Z value','-log10(P) value'});
  resampleFweMethodMenu = putOnTopOfList(params.resampleFweMethod,{'Single-step','Adaptive Single-step','Step-down','Adaptive Step-down'});
  fdrAssumptionMenu = putOnTopOfList(params.fdrAssumption,{'Independence/Positive dependence','None'});
  fweMethodMenu = putOnTopOfList(params.fweMethod,{'Single-step (Bonferroni)','Adaptive Single-step','Step-down (Holm)','Adaptive Step-down','Step-up (Hochberg)','Hommel'});
  fdrMethodMenu = putOnTopOfList(params.fdrMethod,{'Step-up','Adaptive Step-up','Two-stage Step-up','Multiple-stage Step-down','Multiple-stage Step-up'});
  trueNullsEstimationMethodMenu = putOnTopOfList(params.trueNullsEstimationMethod,{'Lowest Slope','Raw P-values Cut-off','Least Squares'});
  
  paramsInfo = {...
      {'testOutput', testOutputMenu,'type=popupmenu', 'Type of statistics for output overlay.  P: outputs the probability value associated with the statistic. p-values less than 1e-16 will be replaced by 0; Z: outputs standard normal values associated with probability p. Z values with a probability less than 1e-16 will be replaced by +/-8.209536145151493'},...
      {'outputStatistic', params.outputStatistic,'type=checkbox', 'If you want to output the statistic (T/F) as an overlay as well'},...
      {'bootstrapIntervals', params.bootstrapIntervals,'type=checkbox', 'Whether to compute Bootstrap confidence intervals for contrasts and parameter estimates'},...
      {'alphaConfidenceIntervals', params.alphaConfidenceIntervals,'contingent=bootstrapIntervals', 'minmax=[0 1]', 'Confidence Intervals will be computed as fractiles (alpha/2) and (1-alpha/2) of the bootstrap estimated null distribution'},...
      {'resampleFweMethod', resampleFweMethodMenu,'type=popupmenu', 'bootstrap/permutation-based family-wise error control method. Adaptive methods estimate the number of true Null hypotheses form the raw p-values. Adaptive Step-down is the most powerful'},...
      {'fweAdjustment', params.fweAdjustment,'type=checkbox', 'Whether to perform adjustments for Bonferroni-type familywise error rate (FWER) control'},...
      {'fweMethod', fweMethodMenu,'contingent=fweAdjustment','type=popupmenu', 'Bonferroni-type FWER control method. Hommel method is the most powerful of non-adaptive methods. Adaptive methods estimate the number of true Null hypotheses from the raw p-values, making the procedure more powerful, but are only implemented for single-step and step-up methods.'},...
      {'fdrAdjustment', params.fdrAdjustment,'type=checkbox', 'Whether to perform adjustments for false discovery rate (FDR) control'},...
      {'fdrAssumption', fdrAssumptionMenu,'contingent=fdrAdjustment','type=popupmenu', 'Distributional assumption for the FDR adjustment methods. Most FDR methods assume that tests are independent or positively correlated, although some have a correcting factor in case no such assumption is made'},...
      {'fdrMethod', fdrMethodMenu,'contingent=fdrAdjustment','type=popupmenu', 'Type of FDR control adjustment method. The multistep adaptive method is supposedly the most powerful'},...
      {'trueNullsEstimationMethod', trueNullsEstimationMethodMenu,'type=popupmenu', 'Method to estimate the number of true null outcomes for adaptive adjustment methods'},...
      {'trueNullsEstimationThreshold', params.trueNullsEstimationThreshold, 'Raw p-values cutoff for true nulls estimation'},...
      {'maskContrastOverlay', params.maskContrastOverlay,'type=checkbox', 'Whether to mask the contrast overlay(s) with the results of the test outcome'},...
      {'noDoubleBootstrap', params.noDoubleBootstrap,'type=checkbox', 'Prevents boostrap-based FWER control adjustment on bootstrap P-value to spare memory'},...
       };

  if ~params.showAdvancedStatisticMenu || useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo,'Advanced Statistics Menu');
  end

  % user hit cancel
  if isempty(tempParams)
    params = tempParams;
    return;
  end
  
  params = mrParamsCopyFields(tempParams,params);
  
  if (params.computeTtests || params.numberFtests) && params.permutationTests && ...
      ischar(params.scanParams{params.scanNum(1)}.stimDuration) && strcmp(params.scanParams{params.scanNum(1)}.stimDuration,'fromFile')
    mrWarnDlg('(getTestParamsGUI) Permutation tests are not compatible with stimulus duration from log file','Yes');
  elseif params.TFCE && params.parametricTests && (params.bootstrapFweAdjustment || params.permutationFweAdjustment) && ismember(params.resampleFweMethod,{'Step-down','Adaptive Step-down'})
    mrWarnDlg('(getTestParamsGUI) Step-down resample-based FWE adjustment is not implemented for TFCE-transformed data','Yes');
  elseif params.TFCE && params.parametricTests && (params.bootstrapFweAdjustment || params.permutationFweAdjustment) && ismember(params.resampleFweMethod,{'Adaptive Single-step'})
    mrWarnDlg('(getTestParamsGUI) Adaptive resample-based FWE adjustment is not implemented for TFCE-transformed data','Yes');
  elseif params.fdrAdjustment && strcmp(params.fdrAssumption,'None') && ismember(params.fdrMethod,{'Adaptive Step-down','Multiple-stage Adaptive Step-up'})
    mrWarnDlg('(getTestParamsGUI) Multi-stage and Step-down adaptive FDR adjustments require the assumption that voxel are independent or positively correlated','Yes');
  else
    keepAsking = 0;
  end
  if keepAsking && useDefault
    params = [];
    return;
  end

end


