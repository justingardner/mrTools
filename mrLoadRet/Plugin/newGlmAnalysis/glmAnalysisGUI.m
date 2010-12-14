% glmAnalysisGUI.m
%
%      usage: params = glmAnalysisGUI(varargin)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 08/01/2010
%    purpose: return parameters for GLM analysis
%           $Id: glmAnalysisGUI.m 1921 2010-12-13 10:22:30Z julien $

function params = glmAnalysisGUI(varargin)   

% get the arguments
eval(evalargs(varargin));

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('params'),params = struct;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;

%convert old params
params = convertOldGlmParams(params);

% rm the tseriesFile field because we're creating a new analysis. It will be merged later if required
if isfield(params,'tseriesFile')
  params = mrParamsRemoveField(params,'tseriesFile');
end

% put default parameters if not already set
if ~isfield(params,'hrfParams')
  params.hrfParams = struct;
end
if ~isfield(params,'groupName') || isempty(params.groupName)
  params.groupName = groupName;
end
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = [];
end
if ~isfield(params,'saveName') || isempty(params.saveName)
  params.saveName = 'GLM';
else
  params.saveName = viewGet(thisView,'analysisName');
end
if ~isfield(params,'hrfModel') || isempty(params.hrfModel)
  params.hrfModel ='hrfDoubleGamma';
end
if ~isfield(params,'analysisVolume') || isempty(params.analysisVolume)
  params.analysisVolume ='Whole volume';
end

if ~isfield(params,'numberEVs') || isempty(params.numberEVs)
  params.numberEVs = 0;
end
if ~isfield(params,'numberContrasts') || isempty(params.numberContrasts)
  if isfield(params,'contrasts')
    params.numberContrasts = size(params.contrasts,1);
  else
    params.numberContrasts = 0;
  end
end
if ~isfield(params,'computeTtests') || isempty(params.computeTtests)
  params.computeTtests = 0;
end
if ~isfield(params,'numberFtests') || isempty(params.numberFtests)
  if isfield(params,'params') && isfield(params,'restrictions')
    params.numberFtests = length(params.restrictions);
  else
    params.numberFtests = 0;
  end
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
hrfModelMenu = putOnTopOfList(params.hrfModel,{'hrfDoubleGamma','hrfDeconvolution'});
analysisVolumeMenu = {'Whole volume','Subset box','Loaded ROI(s)'};
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
    {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','bla'},...
    {'numberContrasts',params.numberContrasts,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','bla'},...
    {'computeTtests', params.computeTtests,'type=checkbox', 'Whether contrasts are tested for significance'},...
    {'numberFtests',params.numberFtests,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','bla'},...
    {'covCorrection',params.covCorrection,'type=checkbox','Correction for correlated noise'},...
    {'correctionType',correctionTypeMenu,'type=popupmenu','contingent=covCorrection','Type of correction (see Wiki for the different algorithms and computing time comparison)'},...
    {'covEstimation',covEstimationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of Estimation of the noise covariance matrix'},...
    {'covEstimationAreaSize',params.covEstimationAreaSize, 'minmax=[1 inf]','contingent=covCorrection','round=1','dimensions in voxels of the spatial window on which the covariance matrix is estimated (in the X and Y dimensions)'},...
    {'covFactorization',covFactorizationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of factorization of the covariance matrix for the computation of the pre-whitening filter'},...
    {'nonLinearityCorrection',params.nonLinearityCorrection,'visible=0','type=checkbox','Correction for non linearity. if Yes, the response saturates when it reaches some value defined below'},...
    {'saturationThreshold',params.saturationThreshold,'visible=0', 'minmax=[1 inf]','contingent=nonLinearityCorrection','Saturation threshold, expressed in terms of the maximum amplitude of the model HRF (e.g. 2 means that the response saturates when it reaches twice the maximum of the model HRF '},...
   };

  % Get parameter values
  if defaultParams
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
%   fieldNames=fieldnames(tempParams);
%   for i_field = 1:length(fieldNames)
%     params.(fieldNames{i_field})=eval(['tempParams.' fieldNames{i_field}]);
%   end
  params = copyFields(tempParams,params);


  %perform some checks
  if 0
    %perform check on parameters here if needed
  elseif 0
    %perform check on parameters here if needed
  else
    while askForParams       % get hrf model parameters
      %here we assume that all scans in this group have the same framePeriod
      hrfParams = feval(params.hrfModel,params.hrfParams,viewGet(thisView,'framePeriod',1,viewGet(thisView,'groupNum',params.groupName)),1,defaultParams);
      % if empty user hit cancel, go back
      if isempty(hrfParams)
        break;
      end
      params.hrfParams = hrfParams;
      if ~isfield(params.hrfParams, 'description')
         params.hrfParams.description = params.hrfModel;
      end

      while askForParams    % get the scan numbers
        groupNum = viewGet(thisView,'groupnum',params.groupName);
        nScans=viewGet(thisView,'nScans',groupNum);
        if nScans == 1
          params.scanNum = 1;
        elseif ~ieNotDefined('scanList')
          params.scanNum = scanList;
        elseif defaultParams
          params.scanNum = 1:nScans;
        elseif viewGet(thisView,'nScans',groupNum) >1
          scanNums = selectScans(thisView,[],groupNum,params.scanNum);
          if isempty(scanNums)
            askForParams = 1;
            break;
          else
            params.scanNum = scanNums;
          end
        end

        while askForParams    % get the parameters for each scan
          [scanParams, params] = getGlmScanParamsGUI(thisView,params,defaultParams);
          if isempty(scanParams)
            askForParams = 1;
            break;
          end
          params.scanParams = scanParams;

          while askForParams    %get the stim to EV matrices for each scan
            if ~params.numberEVs
              thisDefaultParams =1;
            else
              thisDefaultParams = defaultParams;
            end
            [scanParams, params] = getGlmEVParamsGUI(thisView,params,thisDefaultParams);
            if isempty(scanParams)
              askForParams = 1;
              break;
            end
            params.scanParams = scanParams;
            if ~params.numberContrasts && ~params.numberFtests
              defaultParams =1;
            end

            while askForParams           %get params
              tempParams = getGlmTestParamsGUI(thisView,params,defaultParams);
              % if empty, user hit cancel, go back
              if isempty(tempParams)
                askForParams = 1;
                break;
              else
                params = tempParams;
                %update the number of tests in case they've changed
                params.numberContrasts = size(params.contrasts,1);
                params.numberFtests = length(params.restrictions);
                askForParams = 0;
              end
            end
          end
        end
        if nScans == 1 || ~ieNotDefined('scanList') || defaultParams
          break;
        end
      end
    end
  end
end

% set the scan number
for i = 1:length(params.scanNum)
  params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
end

