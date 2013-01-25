% getGlmTestParamsGUI.m
%
%        $Id$
%      usage: params = getGlmTestParamsGUI(params,useDefault)
%         by: julien besle, 
%       date: 09/11/2010
%    purpose: return statistical test parameters for GLM analysis
%

function params = getGlmTestParamsGUI(thisView,params,useDefault)

keepAsking = 1;

while keepAsking
  %Default params
  if isfield(params, 'contrasts') && strcmp(params.contrasts,'all')
    params.numberContrasts = params.numberEVs;
  end
  if fieldIsNotDefined(params, 'contrasts') || ~isnumeric(params.contrasts) || ...
    ~isequal(size(params.contrasts,1),params.numberContrasts) || ~isequal(size(params.contrasts,2),params.numberEVs)
      params.contrasts=eye(params.numberContrasts,params.numberEVs);
  end
  if fieldIsNotDefined(params,'tTestSide')
    params.tTestSide = 'Both';
  end

  if fieldIsNotDefined(params, 'fTestNames') || ~isequal(length(params.fTestNames),params.numberFtests) 
    if params.numberFtests
      params.fTestNames=cellstr(reshape(sprintf('fTest%2d',1:params.numberFtests),7,params.numberFtests)');
    else
      params.fTestNames={};
    end
  end
  if fieldIsNotDefined(params, 'restrictions') || ...
    ~isequal(length(params.restrictions),params.numberFtests) || ~isequal(size(params.restrictions{1},2),params.numberEVs)
    params.restrictions = {};
    for iFtest = 1:params.numberFtests
        params.restrictions{iFtest}=zeros(params.numberEVs,params.numberEVs);
    end
  end

  %create model HRF to get the number of components per EV
  %here we assume that all scans in this group have the same sampling parameters 
  %(which are needed to determine the number of components in the deconvolution case) 
  framePeriod = viewGet(thisView,'framePeriod',params.scanNum(1),viewGet(thisView,'groupNum',params.groupName));
  if ~fieldIsNotDefined(params.scanParams{params.scanNum(1)},'estimationSupersampling')
    estimationSupersampling = params.scanParams{params.scanNum(1)}.estimationSupersampling;
  else
    estimationSupersampling = 1;
  end
  [hrfParams,hrf] = feval(params.hrfModel, params.hrfParams,framePeriod/estimationSupersampling,[],1);
  nComponents = size(hrf,2);
  if fieldIsNotDefined(params, 'componentsToTest') || ~isequal(nComponents,length(params.componentsToTest));
    params.componentsToTest = ones(1,nComponents);
  end
  if fieldIsNotDefined(params, 'componentsCombination')
    params.componentsCombination = 'Or';
  end

  if params.computeTtests || params.numberFtests
    if fieldIsNotDefined(params,'parametricTests')
      params.parametricTests = 1;
    end
    if fieldIsNotDefined(params,'bootstrapTests')
      params.bootstrapTests = 0;
    end
    if fieldIsNotDefined(params,'permutationTests')
      params.permutationTests = 0;
    end
    if fieldIsNotDefined(params, 'bootstrapFweAdjustment')
      params.bootstrapFweAdjustment = 0;
    end
    if fieldIsNotDefined(params, 'permutationFweAdjustment')
      params.permutationFweAdjustment = 0;
    end
    if fieldIsNotDefined(params, 'fweAdjustment')
      params.fweAdjustment = 1;
    end
    if fieldIsNotDefined(params, 'fdrAdjustment')
      params.fdrAdjustment = 1;
    end
    if fieldIsNotDefined(params, 'TFCE')
      params.TFCE = 0;
    end
    if fieldIsNotDefined(params, 'nResamples')
      params.nResamples = 10000;
    end
    if fieldIsNotDefined(params, 'showAdvancedStatisticMenu') 
      params.showAdvancedStatisticMenu = 0;
    end
  else
    params.parametricTests = 0;
    params.bootstrapTests = 0;
    params.permutationTests = 0;
    params.bootstrapFweAdjustment = 0;
    params.permutationFweAdjustment = 0;
    params.fweAdjustment = 0;
    params.fdrAdjustment = 0;
    params.TFCE = 0;
    params.nResamples = 1;
    params.showAdvancedStatisticMenu = 0;
  end
  
  tTestSideMenu = putOnTopOfList(params.tTestSide,{'Both','Right','Left'});
  componentsCombinationMenu = putOnTopOfList(params.componentsCombination,{'Add','Or'});
  
  contrastOptionsVisible = 'visible=0';
  tTestOptionsVisible = 'visible=0';
  fTestOptionsVisible = 'visible=0';
  componentOptionsVisible = 'visible=0';
  testOptionsVisible = 'visible=0';
  tfceOptionVisible = 'visible=0';
  
  if params.numberContrasts
    contrastOptionsVisible = 'visible=1';
  end
  if params.computeTtests
    tTestOptionsVisible = 'visible=1';
  end
  if params.numberFtests
    fTestOptionsVisible = 'visible=1';
  end
  if params.computeTtests || params.numberFtests
    if nComponents>1
      componentOptionsVisible = 'visible=1';
    end
    testOptionsVisible = 'visible=1';
  end
  if strcmp(mrGetPref('fslPath'),'FSL not installed')
    params.TFCE = 0;
    if params.computeTtests || params.numberFtests
      tfceOptionVisible = 'enable=0';
    end
  else
      tfceOptionVisible = testOptionsVisible;
  end
  
  currentNumberContrasts = params.numberContrasts;
  currentNumberFtests= params.numberFtests;
  
  paramsInfo = {...
      {'numberContrasts',params.numberContrasts,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus', 'Number of contrasts that will be ouput as an overlay. If modifying the value, press OK to redraw the dialog with the new number of contrasts'},...
      {'numberFtests',params.numberFtests,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of F-tests to be computed and output as overlays'},...
      {'EVnames', params.EVnames, 'type=stringarray','editable=0','Name of the EVs'},...
      {'contrasts', params.contrasts,contrastOptionsVisible,...
            'incdec=[-1 1]','incdecType=plusMinus',...
            'Matrix defining the contrasts of interest that will be output as overlays. Each row defines a contrast, which is a linear combination of EVs.'},...
      {'tTestSide', tTestSideMenu,tTestOptionsVisible,'type=popupmenu', 'Sidedness of contrast T-tests (Both = two-sided, Right = one-sided positive, Left = one-sided negative)'},...
       };
  for iFtest = 1:params.numberFtests
    paramsInfo{end+1} = {fixBadChars(sprintf('fTest%2d',iFtest)), params.fTestNames{iFtest},fTestOptionsVisible ,['Name of F-test ' num2str(iFtest)]};
    paramsInfo{end+1} = {fixBadChars(sprintf('restriction%2d',iFtest)), params.restrictions{iFtest},fTestOptionsVisible,...
            'incdec=[-1 1]','incdecType=plusMinus',...
            ['Restriction matrix defining F-test ' num2str(iFtest)]};
  end
  paramsInfo = [paramsInfo {...
      {'componentsToTest', params.componentsToTest,componentOptionsVisible,...
            'Vector defining which EV components are tested. Put zeros to exclude components or a weight to include them. '},...
      {'componentsCombination', componentsCombinationMenu,componentOptionsVisible,'type=popupmenu', 'How to combine EV components. ''Add'' adds the weighted components into a single EV for contrasts/F-test. ''Or'' ignores the weights and tests contrasts/F-tests at any component that is not 0. Note that ''Or'' is not compatible with one-sided T-tests'}...
      {'TFCE', params.TFCE,tfceOptionVisible,'type=checkbox', 'Performs Threshold Free Cluster Enhancement on T/F maps using fslmaths. This option is only enabled if a path is specified for FSL by running mrSetPref(''fslPath'',''yourpath''). In addition it can only be used in conjonction with bootstrap/permutation tests or resample-based FWER control adjustment.'},...
      {'parametricTests', params.parametricTests,testOptionsVisible,'type=checkbox', 'Performs parametric tests on contrasts/F values'},...
      {'bootstrapTests', params.bootstrapTests,testOptionsVisible,'type=checkbox', 'Performs non-parametric residual bootstrap tests on T/F values. Bootstrapping consists in resampling the residuals with replacement after OLS/GLS fit, using these bootstrapped residuals as the new time-series for OLS/GLS fitting and estimating the null-hypothesis distributions for nResamples resamplings.'},...
      {'permutationTests', params.permutationTests,testOptionsVisible,'type=checkbox', 'Performs non-parametric permutation tests on contrasts/F values. Permutations constist in shuffling stimulus event labels, re-fitting the GLM using OLS/GLS and estimating the null-hypothesis distribution for each statistic for nResamples permutations.'},...
      {'nResamples', params.nResamples,testOptionsVisible, 'minmax=[10 inf]', 'Number of iterations for bootstrap/permutation tests.'},...
      {'bootstrapFweAdjustment', params.bootstrapFweAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts P-values for familywise error rate control by residuals bootstrapping (see advanced menu to change the method; default is the step-down method)'},...
      {'permutationFweAdjustment', params.permutationFweAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts P-values for familywise error rate control by permutations (see advanced menu to change the method; default is the adaptive step-down method)'},...
      {'fweAdjustment', params.fweAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts P-values for familywise error rate control by Bonferroni-type methods (see advanced menu to change the method; default is the adaptive step-down method)'},...
      {'fdrAdjustment', params.fdrAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts P-values for false discovery rate control (see advanced menu to change the method; default is the adaptive step-up method)'},...
      {'showAdvancedStatisticMenu', params.showAdvancedStatisticMenu,testOptionsVisible,'type=checkbox', 'If you want to change adjustment methods or their parameters, compute bootstrap confidence intervals, change the default output of hypothesis tests...'},...
       }];

  if useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo,'Define Contrasts, T-tests and F-tests');
  end

  % user hit cancel
  if isempty(tempParams)
    params = tempParams;
    return;
  end
  
  params = mrParamsCopyFields(tempParams,params);
  %check that contrasts are not empty
  actualNumberContrasts=0;
  if isempty(params.contrasts) %mrParamsDialog returns an empty string instead of an empty array
    params.contrasts = [];
  else
    for iContrast = 1:currentNumberContrasts
      if ~any(params.contrasts(iContrast,:))
        mrWarnDlg('(getGlmTestParamsGUI) Discarding empty contrast');
      else
        actualNumberContrasts = actualNumberContrasts+1;
        params.contrasts(actualNumberContrasts,:) = params.contrasts(iContrast,:);
      end
    end
    params.contrasts = params.contrasts(1:actualNumberContrasts,:);
  end
  allContrasts = params.contrasts;
  
  %check that F-tests are not empty
  restrictions = {};
  fTestNames={};
  actualNumberFtests=0;
  for iFtest = 1:currentNumberFtests
    thisRestriction=params.(fixBadChars(sprintf('restriction%2d',iFtest)));
    params = mrParamsRemoveField(params,fixBadChars(sprintf('restriction%2d',iFtest)));
    if ~any(any(thisRestriction))
      mrWarnDlg('(getGlmTestParamsGUI) Discarding F-test with empty restriction matrix');
    else
      actualNumberFtests = actualNumberFtests+1;
      fTestNames{actualNumberFtests} = params.(fixBadChars(sprintf('fTest%2d',iFtest)));
      restrictions{actualNumberFtests,1} = thisRestriction;
      allContrasts = [allContrasts;thisRestriction];
    end
    params = mrParamsRemoveField(params,fixBadChars(sprintf('fTest%2d',iFtest)));
  end
  params = mrParamsCopyFields(mrParamsDefault({...
  {'fTestNames',fTestNames,'Self-explanatory'},...
  {'restrictions',restrictions,'Restrictions matrix defining F-tests. Each line of each matrix defines a contrast of EVs'},...
  }),params);
  params.restrictions = restrictions; %copy this again, as mrParamsDefault changes the value of these two parameters
  params.fTestNames = fTestNames; %but we use it to create a paramInfo entry for them, which gives acces to the hlep info later on
  

  %this is because of the incoherent behaviour of mrParamsGet that empties disabled params fields
  if isempty(params.TFCE)
    params.TFCE = 0;
  end
  
  if params.numberContrasts && params.computeTtests  && ~strcmp(params.tTestSide,'Both') && ...
      nnz(params.componentsToTest)>1 && strcmp(params.componentsCombination,'Or')
    mrWarnDlg('(getTestParamsGUI) One-sided T-tests on several EV components with ''Or'' combination are not implemented','Yes');
  elseif params.bootstrapTests  && ~ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
    mrWarnDlg('(getTestParamsGUI) Bootstrap tests are currently only allowed for ROI(s) analyses','Yes');
  elseif params.bootstrapFweAdjustment && params.bootstrapTests  && ~ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
    mrWarnDlg('(getTestParamsGUI) Resampled-based FWE adjustments are currently only allowed for ROI(s) analyses','Yes');
  elseif params.permutationTests && ~all( (sum(logical(allContrasts),2)==2 & sum(allContrasts,2)==0) | ~any(allContrasts,2))
    mrWarnDlg('(getTestParamsGUI) Randomization tests can only be run if all contrasts of all T/F tests are 1 to 1 comparisons','Yes');
  elseif params.TFCE && ~(params.bootstrapTests || params.bootstrapFweAdjustment || params.permutationTests || params.permutationFweAdjustment)
    mrWarnDlg('(getTestParamsGUI) TFCE requires resampling tests (bootstrap, permutation or resample-based FWE control)','Yes');
  elseif params.permutationTests && params.permutationFweAdjustment && ~params.parametricTests
    mrWarnDlg('(getTestParamsGUI) permutation-based adjustment for permutation tests is not implemented','Yes');
  elseif params.bootstrapFweAdjustment && ~(params.bootstrapTests ||params.parametricTests)
    mrWarnDlg('(getTestParamsGUI) bootstrap-based FWE adjustments must be performed on bootsrap or parametric tests','Yes');
  elseif params.permutationFweAdjustment && ~(params.permutationTests ||params.parametricTests)
    mrWarnDlg('(getTestParamsGUI) permutation-based FWE adjustments must be performed on permutation or parametric tests','Yes');
  elseif (params.fweAdjustment || params.fdrAdjustment) && ~(params.permutationTests ||params.parametricTests||params.bootstrapTests)
    mrWarnDlg('(getTestParamsGUI) FWE/FDR adjustments must be performed on bootsrap/permutation/parametric tests','Yes');
  elseif params.numberContrasts~=currentNumberContrasts || params.numberFtests~=currentNumberFtests 
    %check if the number of contrasts and F-tests has changed and change the contrasts/restrictions matrices accordingly
    if params.numberContrasts<=size(params.contrasts,1)
      params.contrasts = params.contrasts(1:params.numberContrasts,:);
    else
      params.contrasts = [params.contrasts;zeros(params.numberContrasts-size(params.contrasts,1),size(params.contrasts,2))];
    end
    if params.numberFtests<=length(params.restrictions)
      params.restrictions = params.restrictions(1:params.numberFtests);
      params.fTestNames=params.fTestNames(1:params.numberFtests);
    else
      for iFtest = length(params.restrictions)+1:params.numberFtests
        params.restrictions{iFtest} = zeros(params.numberEVs,params.numberEVs);
        params.fTestNames{iFtest}=sprintf('fTest%2d',iFtest);
      end
    end
    disp('User changed the number of contrasts/F-tests')
  else
    params.numberFtests = actualNumberFtests;
    params.numberContrasts = actualNumberContrasts;
    tempParams = getGlmAdvancedStatisticParamsGUI(thisView,params,useDefault);
    % if empty, user hit cancel, go back
    if isempty(tempParams)
      if size(tempParams,2) %if the top close button has been pressed
        params=tempParams;
        return
      else
        keepAsking = 1;
      end
    else
      params = tempParams;
      keepAsking = 0;
    end
  end
  if keepAsking && useDefault
    params = [];
    return;
  end

end



