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
  if ~isfield(params, 'contrasts') || isempty(params.contrasts) || ...
    ~isequal(size(params.contrasts,1),params.numberContrasts) || ~isequal(size(params.contrasts,2),params.numberEVs)
      params.contrasts=eye(params.numberContrasts,params.numberEVs);
  end
  if ~isfield(params,'tTestSide') || isempty(params.tTestSide)
    params.tTestSide = 'Both';
  end

  if ~isfield(params, 'fTestNames') || isempty(params.fTestNames) || ...
    ~isequal(length(params.fTestNames),params.numberFtests) 
    if params.numberFtests
      params.fTestNames=cellstr(reshape(sprintf('fTest%2d',1:params.numberFtests),7,params.numberFtests)');
    else
      params.fTestNames={};
    end
  end
  if ~isfield(params, 'restrictions') || isempty(params.restrictions) || ...
    ~isequal(length(params.restrictions),params.numberFtests) || ~isequal(size(params.restrictions{1},2),params.numberEVs)
    params.restrictions = {};
    for iFtest = 1:params.numberFtests
        params.restrictions{iFtest}=zeros(params.numberEVs,params.numberEVs);
    end
  end

  %create model HRF
  %here we assume that all scans in this group have the same framePeriod
  [hrfParams,hrf] = feval(params.hrfModel, params.hrfParams,...
    viewGet(thisView,'framePeriod',1,viewGet(thisView,'groupNum',params.groupName))/params.scanParams{params.scanNum(1)}.estimationSupersampling,0,1);
  nComponents = size(hrf,2);
  if ~isfield(params, 'componentsToTest') || isempty(params.componentsToTest) || ~isequal(nComponents,length(params.componentsToTest));
    params.componentsToTest = ones(1,nComponents);
  end
  if ~isfield(params, 'componentsCombination') || isempty(params.componentsCombination)
    params.componentsCombination = 'Or';
  end

  if ~isfield(params,'parametricTests') || isempty(params.parametricTests)
    params.parametricTests = 1;
  end
  if ~isfield(params,'parametricTestOutput')
    params.parametricTestOutput = 'Z value';
  end
  if ~isfield(params, 'TFCE') || isempty(params.TFCE)
      params.TFCE = 0;
  end

  if ~isfield(params,'randomizationTests') || isempty(params.randomizationTests)
    params.randomizationTests = 0;
  end
  if ~isfield(params,'randomizationTestOutput') || isempty(params.randomizationTestOutput)
    params.randomizationTestOutput = 'Z value';
  end
  if ~isfield(params, 'nRand') || isempty(params.nRand)
      params.nRand = 10000;
  end

  if ~isfield(params,'bootstrapStatistics') || isempty(params.bootstrapStatistics)
    params.bootstrapStatistics = 0;
  end
  if ~isfield(params,'bootstrapTestOutput') || isempty(params.bootstrapTestOutput)
    params.bootstrapTestOutput = 'Z value';
  end
  if ~isfield(params, 'nBootstrap') || isempty(params.nBootstrap)
      params.nBootstrap = 10000;
  end
  if ~isfield(params, 'nBootstrap') || isempty(params.nBootstrap)
      params.nBootstrap = 10000;
  end
  if ~isfield(params, 'bootstrapIntervals') || isempty(params.bootstrapIntervals)
      params.bootstrapIntervals = 0;
  end
  if ~isfield(params, 'alphaConfidenceIntervals') || isempty(params.alphaConfidenceIntervals)
      params.alphaConfidenceIntervals = 0.05;
  end

  tTestSideMenu = putOnTopOfList(params.tTestSide,{'Both','Right','Left'});
  componentsCombinationMenu = putOnTopOfList(params.componentsCombination,{'Add','Or'});
  parametricTestOutputMenu = putOnTopOfList(params.parametricTestOutput,{'T/F value','P value','Z value'});
  randomizationTestOutputMenu = putOnTopOfList(params.randomizationTestOutput,{'P value','Z value'});
  bootstrapTestOutputMenu = putOnTopOfList(params.bootstrapTestOutput,{'P value','Z value'});
  
  contrastOptionsVisible = 'visible=0';
  tTestOptionsVisible = 'visible=0';
  fTestOptionsVisible = 'visible=0';
  componentOptionsVisible = 'visible=0';
  testOptionsVisible = 'visible=0';
  tfceOptionVisible = 'visible=0';
  tfceContigency = 'contingent=parametricTests';
  
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
      tfceContigency = '';
    end
  else
      tfceOptionVisible = testOptionsVisible;
  end
  
  paramsInfo = {...
      {'EVnames', params.EVnames, 'type=stringarray','editable=0','Name of the EVs'},...
      {'contrasts', params.contrasts,contrastOptionsVisible,...
            'incdec=[-1 1]','incdecType=plusMinus',...
            'Matrix defining the contrasts of interest that will be output as overlays. Each row defines a contrast, which is a linear combination of EVss. Contrasts are computed after running the preprocessing function, so the number of colums should match the the number of EVs after preprocessing'},...
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
      {'parametricTests', params.parametricTests,testOptionsVisible,'type=checkbox', 'Performs parametric tests on contrasts/F values'},...
      {'parametricTestOutput', parametricTestOutputMenu,testOptionsVisible,'contingent=parametricTests','type=popupmenu', 'Type of statistics for output overlay. T/F: outputs the value of the statistic (T for contrasts and F for F-tests); P: outputs the probability value associated with the statistic. p-values less than 1e-16 will be replaced by 0; Z: outputs standard normal values associated with probability p. Z values with a probability less than 1e-16 will be replaced by +/-8.209536145151493'},...
      {'TFCE', params.TFCE,tfceContigency,tfceOptionVisible,'type=checkbox', 'Performs Threshold Free Cluster Enhancement on T/F maps using fslmaths. This option is only enabled if a path is specified for FSL by running mrSetPref(''fslPath'',''yourpath'')'},...
      {'randomizationTests', params.randomizationTests,testOptionsVisible,'type=checkbox', 'Performs non-parametric randomization tests on contrasts/F values'},...
      {'nRand', params.nRand,testOptionsVisible,'contingent=randomizationTests', 'minmax=[10 inf]', 'Number of randomizations for randomization tests. Randomizations are implemented by shuffling event labels before running the pre-processing function'},...
      {'randomizationTestOutput', randomizationTestOutputMenu,testOptionsVisible,'contingent=randomizationTests','type=popupmenu', 'p: Outputs the frequency of randomizations giving Contrast/F values greater than the actual value; Z: outputs values for the standard normal distribution associated with the probability p.'},...
      {'bootstrapStatistics', params.bootstrapStatistics,testOptionsVisible,'type=checkbox', 'Performs non-parametric residual bootstrap tests on T/F values or confidence intervals on contrasts and parameter estimates. Bootstrapping consists in resampling the residuals with replacement after OLS/GLS fit and estimating the null-hypothesis distributions using the bootstrapped residuals as the new time-series.'},...
      {'nBootstrap', params.nBootstrap,testOptionsVisible,'contingent=bootstrapStatistics', 'minmax=[10 inf]', 'Number of resamplings for bootstrap tests.'},...
      {'bootstrapTestOutput', bootstrapTestOutputMenu,testOptionsVisible,'contingent=bootstrapStatistics','type=popupmenu', 'p: Outputs the frequency of bootstrap resamples giving Contrast/F values greater than the actual value; Z: outputs values for the standard normal distribution associated with the probability p.'},...
      {'bootstrapIntervals', params.bootstrapIntervals,testOptionsVisible,'contingent=bootstrapStatistics','type=checkbox', 'Whether to compute Bootstrap confidence intervals for contrasts and parameter estimates'},...
      {'alphaConfidenceIntervals', params.alphaConfidenceIntervals,testOptionsVisible,'contingent=bootstrapIntervals', 'minmax=[0 1]', 'Confidence Intervals will be computed as fractiles (alpha/2) and (1-alpha/2) of the bootstrap estimated null distribution'},...
       }];

  if useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo,'Define Contrasts, T-tests and F-tests');
    %buttonWidth = min(5/params.numberEVs,.4);
    %tempParams = mrParamsDialog(paramsInfo,'Define Explanatory Variables, Contrasts and F-tests','buttonWidth',buttonWidth);
  end

  % user hit cancel
  if isempty(tempParams)
    params = tempParams;
    return;
  end
  
  params = mrParamsCopyFields(tempParams,params);
  %check that contrasts are not empty
  actualNumberContrasts=0;
  for iContrast = 1:size(params.contrasts,1)
    if ~any(params.contrasts(iContrast,:))
      mrWarnDlg('(getGlmTestParamsGUI) Discarding empty contrast');
    else
      actualNumberContrasts = actualNumberContrasts+1;
      params.contrasts(actualNumberContrasts,:) = params.contrasts(iContrast,:);
    end
  end
  params.contrasts = params.contrasts(1:actualNumberContrasts,:);
  
  %check that F-tests are not empty
  restrictions = {};
  fTestNames={};
  actualNumberFtests=0;
  for iFtest = 1:params.numberFtests
    thisRestriction=params.(fixBadChars(sprintf('restriction%2d',iFtest)));
    params = mrParamsRemoveField(params,fixBadChars(sprintf('restriction%2d',iFtest)));
    if ~any(any(thisRestriction))
      mrWarnDlg('(getGlmTestParamsGUI) Discarding F-test with empty restriction matrix');
    else
      actualNumberFtests = actualNumberFtests+1;
      fTestNames{actualNumberFtests} = params.(fixBadChars(sprintf('fTest%2d',iFtest)));
      restrictions{actualNumberFtests} = thisRestriction;
    end
    params = mrParamsRemoveField(params,fixBadChars(sprintf('fTest%2d',iFtest)));
  end
  params = mrParamsCopyFields(mrParamsDefault({...
  {'fTestNames',fTestNames,'Self-explanatory'},...
  {'restrictions',restrictions,'Restrictions matrix defining F-tests. Each line of each matrix defines a contrast of EVs'},...
  }),params);
  params.restrictions = restrictions; %copy this again, as mrParamsDefault changes the value of thses two parameters
  params.fTestNames = fTestNames; %but we use it to create a paramInfo entry for them, which gives acces to the hlep info later on

  %this is because of the incoherent behaviour of mrParamsGet that empties disabled params fields
  if isempty(params.TFCE)
    params.TFCE = 0;
  end
    
  if (params.computeTtests || params.numberFtests) && ...
      params.randomizationTests && ischar(params.scanParams{params.scanNum(1)}.stimDuration) && strcmp(params.scanParams{params.scanNum(1)}.stimDuration,'fromFile')
  elseif params.numberContrasts && params.computeTtests  && ~strcmp(params.tTestSide,'Both') && ...
      length(params.componentsToTest)>1 && strcmp(params.componentsCombination,'Or')
     mrWarnDlg('(getTestParamsGUI) One-sided T-tests on several EV components with ''Or'' combination are not implemented','Yes');
  else
     keepAsking = 0;
  end
  if keepAsking && useDefault
    params = [];
    return;
  end

end


