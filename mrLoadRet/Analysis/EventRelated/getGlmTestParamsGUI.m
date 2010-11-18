% getGlmTestParamsGUI.m
%
%        $Id$
%      usage: testParams = getGlmTestParamsGUI(params,useDefault)
%         by: julien besle, 
%       date: 09/11/2010
%    purpose: return statistical test parameters for GLM analysis
%

function testParams = getGlmTestParamsGUI(thisView,params,useDefault)

keepAsking = 1;
if isfield(params,'testParams')
   testParams = params.testParams;
else
   % make the output as long as the number of scans
   testParams = struct;
end

while keepAsking
  keepAsking = 0;
  %Default params
  if ~isfield(testParams, 'stimToEVmatrix') || isempty(testParams.contrasts) || ...
    ~isequal(size(testParams.stimToEVmatrix,1),params.numberEvents) || ~isequal(size(testParams.stimToEVmatrix,2),params.numberEVs)
      testParams.stimToEVmatrix = eye(params.numberEvents,params.numberEVs);
  end
  if ~isfield(testParams, 'EVnames') || isempty(testParams.EVnames) || ...
    ~isequal(length(testParams.EVnames),size(testParams.stimToEVmatrix,2)) 
      testParams.EVnames = repmat({' '},1,params.numberEVs);
      if params.numberEVs < params.numberEvents
        testParams.EVnames = params.stimNames(1:params.numberEVs);
      else
        testParams.EVnames(1:params.numberEVs) = params.stimNames;
      end
  end
  if ~isfield(testParams, 'contrasts') || isempty(testParams.contrasts) || ...
    ~isequal(size(testParams.contrasts,1),params.numberContrasts) || ~isequal(size(testParams.contrasts,2),params.numberEVs)
      testParams.contrasts=eye(params.numberContrasts,params.numberEVs);
  end
  if ~isfield(testParams,'tTestSide') || isempty(testParams.tTestSide)
    testParams.tTestSide = 'Both';
  end

  if ~isfield(testParams, 'fTestNames') || isempty(testParams.fTestNames) || ...
    ~isequal(length(testParams.fTestNames),params.numberFtests) 
      testParams.fTestNames=cellstr(reshape(sprintf('fTest%2d',1:params.numberFtests),7,params.numberFtests)');
  end
  if ~isfield(testParams, 'restrictions') || isempty(testParams.restrictions) || ...
    ~isequal(length(testParams.restrictions),params.numberFtests) || ~isequal(size(testParams.restrictions{1},2),params.numberEVs)
    for iFtest = 1:params.numberFtests
        testParams.restrictions{iFtest}=zeros(params.numberEVs,params.numberEVs);
    end
  end

  %create model HRF
  %here we assume that all scans in this group have the same framePeriod
  [hrfParams,hrf] = feval(params.hrfModel, params.hrfParams,viewGet(thisView,'framePeriod',1,viewGet(thisView,'groupNum',params.groupName)),0,1);
  nComponents = size(hrf,2);
  if ~isfield(testParams, 'componentsToTest') || isempty(testParams.componentsToTest) || ~isequal(nComponents,length(testParams.componentsToTest));
    testParams.componentsToTest = ones(1,nComponents);
  end
  if ~isfield(testParams, 'componentsCombination') || isempty(testParams.componentsCombination)
    testParams.componentsCombination = 'Or';
  end

  if ~isfield(testParams,'parametricTests') || isempty(testParams.parametricTests)
    testParams.parametricTests = 1;
  end
  if ~isfield(testParams,'parametricTestOutput')
    testParams.parametricTestOutput = 'Z value';
  end
  if ~isfield(testParams, 'TFCE') || isempty(testParams.TFCE)
      testParams.TFCE = 0;
  end

  if ~isfield(testParams,'randomizationTests') || isempty(testParams.randomizationTests)
    testParams.randomizationTests = 0;
  end
  if ~isfield(testParams,'randomizationTestOutput') || isempty(testParams.randomizationTestOutput)
    testParams.randomizationTestOutput = 'Z value';
  end
  if ~isfield(testParams, 'nRand') || isempty(testParams.nRand)
      testParams.nRand = 10000;
  end


  tTestSideMenu = putOnTopOfList(testParams.tTestSide,{'Both','Right','Left'});
  componentsCombinationMenu = putOnTopOfList(testParams.componentsCombination,{'Add','Or'});
  parametricTestOutputMenu = putOnTopOfList(testParams.parametricTestOutput,{'T/F value','P value','Z value'});
  randomizationTestOutputMenu = putOnTopOfList(testParams.randomizationTestOutput,{'P value','Z value'});
  
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
    testParams.TFCE = 0;
    if params.computeTtests || params.numberFtests
      tfceOptionVisible = 'enable=0';
      tfceContigency = '';
    end
  else
      tfceOptionVisible = testOptionsVisible;
  end
  
  paramsInfo = cell(1,params.numberEvents);
  for iEvent = 1:params.numberEvents
    paramsInfo{iEvent} = {fixBadChars(params.stimNames{iEvent}),...
      testParams.stimToEVmatrix(iEvent,:),...
      ['How much stimulus ' fixBadChars(params.stimNames{iEvent}) ' contributes to each EV']};
  end
  paramsInfo = [paramsInfo {...
      {'EVnames', testParams.EVnames, 'self explanatory'},...
      {'contrasts', testParams.contrasts,contrastOptionsVisible, 'Matrix defining the contrasts of interest that will be output as overlays. Each row defines a contrast, which is a linear combination of EVss. Contrasts are computed after running the preprocessing function, so the number of colums should match the the number of EVs after preprocessing'},...
      {'tTestSide', tTestSideMenu,tTestOptionsVisible,'type=popupmenu', 'Sidedness of contrast T-tests (Both = two-sided, Right = one-sided positive, Left = one-sided negative)'},...
       }];
  for iFtest = 1:params.numberFtests
    paramsInfo{end+1} = {fixBadChars(sprintf('fTest%2d',iFtest)), testParams.fTestNames{iFtest},fTestOptionsVisible ,['Name of F-test ' num2str(iFtest)]};
    paramsInfo{end+1} = {fixBadChars(sprintf('restriction%2d',iFtest)), testParams.restrictions{iFtest},fTestOptionsVisible, ['Restriction matrix defining F-test ' num2str(iFtest)]};
  end
  paramsInfo = [paramsInfo {...
      {'componentsToTest', testParams.componentsToTest,componentOptionsVisible, 'Vector defining which EV components are tested. Put zeros to exclude components or a weight to include them. '},...
      {'componentsCombination', componentsCombinationMenu,componentOptionsVisible,'type=popupmenu', 'How to combine EV components. ''Add'' adds the weighted components into a single EV for contrasts/F-test. ''Or'' ignores the weights and tests contrasts/F-tests at any component that is not 0. Note that ''Or'' is not compatible with one-sided T-tests'}...
      {'parametricTests', testParams.parametricTests,testOptionsVisible,'type=checkbox', 'Performs parametric tests on contrasts/F values'},...
      {'parametricTestOutput', parametricTestOutputMenu,testOptionsVisible,'contingent=parametricTests','type=popupmenu', 'Type of statistics for output overlay. T/F: outputs the value of the statistic (T for contrasts and F for F-tests); P: outputs the probability value associated with the statistic. p-values less than 1e-16 will be replaced by 0; Z: outputs standard normal values associated with probability p. Z values with a probability less than 1e-16 will be replaced by +/-8.209536145151493'},...
      {'TFCE', testParams.TFCE,tfceContigency,tfceOptionVisible,'type=checkbox', 'Performs Threshold Free Cluster Enhancement on T/F maps using fslmaths. This option is only enabled if a path is specified for FSL by running mrSetPref(''fslPath'',''yourpath'')'},...
      {'randomizationTests', testParams.randomizationTests,testOptionsVisible,'type=checkbox', 'Performs non-parametric randomization tests on contrasts/F values'},...
      {'nRand', testParams.nRand,testOptionsVisible,'contingent=randomizationTests', 'minmax=[1 inf]', 'Number of randomizations for randomization tests. Randomizations are implemented by shuffling event labels before running the pre-processing function'},...
      {'randomizationTestOutput', randomizationTestOutputMenu,testOptionsVisible,'contingent=randomizationTests','type=popupmenu', 'p: Outputs the frequency of randomizations giving Contrast/F values greater than the actual value; Z: outputs values for the standard normal distribution associated with the probability p.'},...
       }];

  if useDefault
    testParams = mrParamsDefault(paramsInfo);
  else
    buttonWidth = min(500/params.numberEVs,50);
    testParams = mrParamsDialog(paramsInfo,'Define Explanatory Variables, Contrasts and F-tests','buttonWidth',buttonWidth);
  end

  % user hit cancel
  if isempty(testParams)
     return;
  end
  
  tTestSideMenu = putOnTopOfList(testParams.tTestSide,tTestSideMenu);
  componentsCombinationMenu = putOnTopOfList(testParams.componentsCombination,componentsCombinationMenu);
  parametricTestOutputMenu = putOnTopOfList(testParams.parametricTestOutput,parametricTestOutputMenu);
  randomizationTestOutputMenu = putOnTopOfList(testParams.randomizationTestOutput,randomizationTestOutputMenu);

  
  for iEvent = 1:params.numberEvents
    testParams.stimToEVmatrix(iEvent,:) = testParams.(fixBadChars(params.stimNames{iEvent}));
    testParams = rmfield(testParams,fixBadChars(params.stimNames{iEvent}));
  end
  for iFtest = 1:params.numberFtests
    testParams.fTestNames{iFtest} = testParams.(fixBadChars(sprintf('fTest%2d',iFtest)));
    testParams.restrictions{iFtest} = testParams.(fixBadChars(sprintf('restriction%2d',iFtest)));
    testParams.restrictions{iFtest} = testParams.restrictions{iFtest}(any(testParams.restrictions{iFtest},2),:); %remove lines of 0
  end
    
  if (params.computeTtests || params.numberFtests) && ...
      testParams.randomizationTests && ischar(params.scanParams{1}.stimDuration) && strcmp(params.scanParams{1}.stimDuration,'fromFile')
     mrWarnDlg('(getTestParamsGUI) Randomization tests are not (yet) compatible with stim duration from files');
     if ~useDefault
       keepAsking = 1;
     else
       testParams = [];
     end
  end
  
  if params.numberContrasts && params.computeTtests  && ~strcmp(testParams.tTestSide,'Both') && ...
      length(testParams.componentsToTest)>1 && strcmp(testParams.componentsCombination,'Or')
     mrWarnDlg('(getTestParamsGUI) One-sided T-tests on several EV components with ''Or'' combination are not implemented','Yes');
     if ~useDefault
       keepAsking = 1;
     else
       testParams = [];
     end
  end

end

