% eventRelatedGlmGUI.m
%
%      usage: params = eventRelatedGlmGUI(varargin)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 08/01/2010
%    purpose: return parameters for GLM analysis
%           $Id$

function params = eventRelatedGlmGUI(varargin)   

% get the arguments
eval(evalargs(varargin));

if ieNotDefined('useDefault'),useDefault = 0;end
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('params'),params = struct;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;

%convert old params
params = convertOldParams(params);


% put default parameters if not already set
if ~isfield(params,'hrfParams')
  params.hrfParams = struct;
end
if ~isfield(params,'groupName') || isempty(params.groupName)
  params.groupName = groupName;
end
if ~isfield(params,'saveName') || isempty(params.saveName)
  params.saveName = 'GLM';
end
if ~isfield(params,'hrfModel') || isempty(params.hrfModel)
  params.hrfModel ='hrfDiffGamma';
end
if ~isfield(params,'analysisVolume') || isempty(params.analysisVolume)
  params.analysisVolume ='Whole volume';
end

if ~isfield(params,'numberEVs') || isempty(params.numberEVs)
  params.numberEVs = 0;
end
if ~isfield(params,'numberContrasts') || isempty(params.numberContrasts)
  params.numberContrasts = 0;
end
if ~isfield(params,'computeTtests') || isempty(params.computeTtests)
  params.computeTtests = 0;
end
if ~isfield(params,'numberFtests') || isempty(params.numberFtests)
  params.numberFtests = 0;
end

if ~isfield(params,'covCorrection') || isempty(params.covCorrection)
   params.covCorrection = 0;
end
if ~isfield(params,'correctionType') || isempty(params.correctionType)
   params.correctionType = 'generalizedLeastSquares';
end
if ~isfield(params,'covEstimation') || isempty(params.covEstimation)
   params.covEstimation = 'singleTukeyTapers';
end
if ~isfield(params,'covEstimationAreaSize') || isempty(params.covEstimationAreaSize) %this is necessary so that it doesn't become a string is it is empty
   params.covEstimationAreaSize = 20;
end
if ~isfield(params,'covFactorization') || isempty(params.covFactorization)
   params.covFactorization = 'Cholesky';
end

if ~isfield(params,'nonLinearityCorrection') || isempty(params.nonLinearityCorrection)
   params.nonLinearityCorrection = 0;
end
if ~isfield(params,'saturationThreshold') || isempty(params.saturationThreshold)
   params.saturationThreshold = 2;
end

askForParams = 1;
% put group name on top of list to make it the default
groupNames = putOnTopOfList(params.groupName,viewGet(thisView,'groupNames'));
hrfModelMenu = putOnTopOfList(params.hrfModel,{'hrfDiffGamma','hrfDeconvolution'});
analysisVolumeMenu = {'Whole volume','Subset box'};
%if viewGet(thisView,'numberofROIs') %This doesn't work if called from recompute because the call does not include the view (thisView)
   analysisVolumeMenu{end+1} = 'Loaded ROI(s)'; 
%end
analysisVolumeMenu = putOnTopOfList(params.analysisVolume,analysisVolumeMenu);
correctionTypeMenu = putOnTopOfList(params.correctionType,{'varianceCorrection','preWhitening','generalizedLeastSquares'});%, 'generalizedFTest'});
covEstimationMenu = putOnTopOfList(params.covEstimation,{'singleTukeyTapers','dampenedOscillator'});
covFactorizationMenu = putOnTopOfList(params.covFactorization,{'Cholesky'});




while askForParams
  % set the parameter string
  paramsInfo = {...
    {'groupName',groupNames,'type=popupmenu','Name of the group from which to do the analysis'},...
    {'saveName',params.saveName,'File name to try to save as'},...
    {'hrfModel',hrfModelMenu,'type=popupmenu','Name of the function that defines the hrf used in glm',},...
    {'analysisVolume',analysisVolumeMenu,'type=popupmenu','Which voxels perform the analysis on. If subsetbox, the dimensions of the box can be specific to each scan and have to include a margin if covEstimationAreaSize>1',},...
    {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','bla'},...
    {'numberContrasts',params.numberContrasts,'minmax=[0 inf]','incdec=[-1 1]','bla'},...
    {'computeTtests', params.computeTtests,'type=checkbox', 'Whether contrasts are tested for significance'},...
    {'numberFtests',params.numberFtests,'minmax=[0 inf]','incdec=[-1 1]','bla'},...
    {'covCorrection',params.covCorrection,'type=checkbox','Correction for correlated noise'},...
    {'correctionType',correctionTypeMenu,'type=popupmenu','contingent=covCorrection','Type of correction (see Wiki for the different algorithms and computing time comparison)'},...
    {'covEstimation',covEstimationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of Estimation of the noise covariance matrix'},...
    {'covEstimationAreaSize',params.covEstimationAreaSize, 'minmax=[1 inf]','contingent=covCorrection','round=1','dimensions in voxels of the spatial window on which the covariance matrix is estimated (in the X and Y dimensions)'},...
    {'covFactorization',covFactorizationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of factorization of the covariance matrix for the computation of the pre-whitening filter'},...
    {'nonLinearityCorrection',params.nonLinearityCorrection,'visible=0','type=checkbox','Correction for non linearity. if Yes, the response saturates when it reaches some value defined below'},...
    {'saturationThreshold',params.saturationThreshold,'visible=0', 'minmax=[1 inf]','contingent=nonLinearityCorrection','Saturation threshold, expressed in terms of the maximum amplitude of the model HRF (e.g. 2 means that the response saturates when it reaches twice the maximum of the model HRF '},...
   };

  % Get parameter values
  if useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo, 'GLM parameters');
  end


  % if empty user hit cancel
  if isempty(tempParams)
    params = [];
    return;
  end
  %replace params in the structure that has been passed in
  hrfModelMenu = putOnTopOfList(tempParams.hrfModel,hrfModelMenu);
  analysisVolumeMenu = putOnTopOfList(tempParams.analysisVolume,analysisVolumeMenu);
  correctionTypeMenu = putOnTopOfList(tempParams.correctionType,correctionTypeMenu);
  covEstimationMenu = putOnTopOfList(tempParams.covEstimation,covEstimationMenu);
  covFactorizationMenu = putOnTopOfList(tempParams.covFactorization,covFactorizationMenu);
  fieldNames=fieldnames(tempParams);
  for i_field = 1:length(fieldNames)
    params.(fieldNames{i_field})=eval(['tempParams.' fieldNames{i_field}]);
  end

  while askForParams       % get hrf model parameters
    %here we assume that all scans in this group have the same framePeriod
    hrfParams = feval(params.hrfModel,params.hrfParams,viewGet(thisView,'framePeriod',1,viewGet(thisView,'groupNum',params.groupName)),1,useDefault);
    % if empty user hit cancel, go back
    if isempty(hrfParams)
      break;
    else
      params.hrfParams = hrfParams;
      if ~isfield(params.hrfParams, 'description')
         params.hrfParams.description = params.hrfModel;
      end

      while askForParams    % get the scan params
        % get scans
        %thisView = viewSet(thisView,'groupName',params.groupName);
        groupNum = viewGet(thisView,'groupnum',params.groupName);
        if ~ieNotDefined('scanList')
           params.scanNum = scanList;
        elseif useDefault
           %params.scanNum = 1:viewGet(thisView,'nScans');
           params.scanNum = 1:viewGet(thisView,'nScans',groupNum);
        elseif viewGet(thisView,'nScans',groupNum) >1
           %params.scanNum = selectScans(thisView);
           params.scanNum = selectScans(thisView,[],groupNum);
           if isempty(params.scanNum)
              askForParams = 1;
              break;
           end
        else
           params.scanNum = 1;
        end

        % get the parameters for each scan
        [scanParams, params] = getGlmScanParamsGUI(thisView,params,useDefault);
        if isempty(scanParams)
          askForParams = 1;
          break;
        else
          params.scanParams = scanParams;
          
          if ~params.numberEVs
            params.numberEVs = params.numberEvents;
            defaultTestParams = 1;
          else
            defaultTestParams = useDefault;
          end
          while askForParams           %get testParams
            testParams = getGlmTestParamsGUI(thisView,params,defaultTestParams);
            % if empty, user hit cancel, go back
            if isempty(testParams)
              askForParams = 1;
              break;
            else
              params.testParams = testParams;
              askForParams = 0;
            end
          end
        end
      end
    end
  end
end

% set the scan number
for i = 1:length(params.scanNum)
  params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
end



function params = convertOldParams(params)

if isfield(params,'contrasts')
  if ischar(params.contrasts)
    params.testParams.contrasts = eval(params.contrasts);
  else
    params.testParams.contrasts = params.contrasts;
  end
  params = rmfield(params,'contrasts');
end
if isfield(params,'f_tests')
  if ischar(params.f_tests)
    params.testParams.fTests = eval(params.f_tests);
  else
    params.testParams.fTests = params.f_tests;
  end
  params = rmfield(params,'f_tests');
end
if isfield(params,'fTests')
  if ischar(params.fTests)
    params.testParams.fTests = eval(params.fTests);
  else
    params.testParams.fTests = params.fTests;
  end
  params = rmfield(params,'fTests');
end
if isfield(params,'n_rand')
  params.testParams.nRand = params.n_rand;
  params = rmfield(params,'n_rand');
end
if isfield(params,'nRand')
  params.testParams.nRand = params.nRand;
  params = rmfield(params,'nRand');
end
if isfield(params,'outputStatistic')
  params = rmfield(params,'outputStatistic');
end
if isfield(params,'outputZStatistic')
  params = rmfield(params,'outputZStatistic');
end
if isfield(params,'outputPValue')
  params = rmfield(params,'outputPValue');
end
if isfield(params,'TFCE')
  params.testParams.TFCE = params.TFCE;
  params = rmfield(params,'TFCE');
end
if isfield(params,'tTestSide')
  params.testParams.tTestSide = params.tTestSide;
  params = rmfield(params,'tTestSide');
end

