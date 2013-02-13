% glmAnalysisGUI.m
%
%      usage: params = glmAnalysisGUI(varargin)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 08/01/2010
%    purpose: return parameters for GLM analysis
%           $Id$

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
nRois = viewGet(thisView,'numberOfRois');
nVisibleRois = length(viewGet(thisView,'visibleRois'));
if ~isfield(params,'analysisVolume') || isempty(params.analysisVolume)
  params.analysisVolume ='Whole volume';
else
  if strcmp(params.analysisVolume,'Loaded ROI(s)')
    if ~nRois
      mrWarnDlg('(glmAnalysisGUI) No ROI loaded, switching to Whole Volume');
      params.analysisVolume ='Whole volume';
    end
  end
  if strcmp(params.analysisVolume,'Visible ROI(s)')
    if ~nVisibleRois
      mrWarnDlg('(glmAnalysisGUI) No ROI is visible, switching to Whole Volume');
      params.analysisVolume ='Whole volume';
    end
  end
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

if fieldIsNotDefined(params,'spatialSmoothing')
   params.spatialSmoothing = 0;
end

if fieldIsNotDefined(params,'covCorrection')
   params.covCorrection = 0;
end
if fieldIsNotDefined(params,'correctionType')
   params.correctionType = 'generalizedLeastSquares';
end
if fieldIsNotDefined(params,'covEstimation')
   params.covEstimation = 'singleTukeyTapers';
end
if fieldIsNotDefined(params,'covEstimationAreaSize') %this is necessary so that it doesn't become a string is it is empty
   params.covEstimationAreaSize = 20;
end
if fieldIsNotDefined(params,'covFactorization')
   params.covFactorization = 'Cholesky';
end
roiNames = viewGet(thisView,'roinames');
if fieldIsNotDefined(params,'covEstimationBrainMask')
   params.covEstimationBrainMask = 'None';
elseif ~ismember(params.covEstimationBrainMask,roiNames)
   mrWarnDlg(['(glmAnalysisGUI) covariance estimation ROI mask ' params.covEstimationBrainMask ' is not loaded, switching to no brain mask']);
   params.covEstimationBrainMask = 'None';
end

if fieldIsNotDefined(params,'nonLinearityCorrection')
   params.nonLinearityCorrection = 0;
end
if fieldIsNotDefined(params,'saturationThreshold')
   params.saturationThreshold = 2;
end

askForParams = 1;
% put group name on top of list to make it the default
groupNames = putOnTopOfList(params.groupName,viewGet(thisView,'groupNames'));
hrfModelMenu = putOnTopOfList(params.hrfModel,{'hrfDoubleGamma','hrfFslFlobs','hrfDeconvolution','hrfBoxcar'});
analysisVolumeMenu = {'Whole volume','Subset box'};
if nRois
  analysisVolumeMenu{end+1} = 'Loaded ROI(s)';
end
if nVisibleRois
  analysisVolumeMenu{end+1} = 'Visible ROI(s)';
end
analysisVolumeMenu = putOnTopOfList(params.analysisVolume,analysisVolumeMenu);
correctionTypeMenu = putOnTopOfList(params.correctionType,{'varianceCorrection','preWhitening','generalizedLeastSquares'});%, 'generalizedFTest'});
covEstimationMenu = putOnTopOfList(params.covEstimation,{'singleTukeyTapers','dampenedOscillator'});
covFactorizationMenu = putOnTopOfList(params.covFactorization,{'Cholesky'});
covEstimationBrainMaskMenu = putOnTopOfList(params.covEstimationBrainMask,[{'None'} roiNames]);


while askForParams
  % set the parameter string
  paramsInfo = {...
    {'groupName',groupNames,'type=popupmenu','Name of the group from which to do the analysis'},...
    {'saveName',params.saveName,'File name to try to save as'},...
    {'hrfModel',hrfModelMenu,'type=popupmenu','Name of the function that defines the Hemodynamic Response Function model that will be convolved with the design matrix',},...
    {'analysisVolume',analysisVolumeMenu,'type=popupmenu','Which voxels to perform the analysis on. If subsetbox, the dimensions of the box can be specific to each scan and have to include a margin if covEstimationAreaSize>1',},...
    {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of Explanatory Variables in the model = number of columns in the design matrix. if 0, the number of EV will be set to the number of stimulus type in the stimulsu file'},...
    {'numberContrasts',params.numberContrasts,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of contrasts that will be ouput as an overlay'},...
    {'computeTtests', params.computeTtests,'type=checkbox', 'Whether contrasts are tested for significance using T-tests'},...
    {'numberFtests',params.numberFtests,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of F-tests to be computed and output as overlays'},...
    {'spatialSmoothing',params.spatialSmoothing, 'minmax=[0 inf]','Width at half-maximum in voxels of a 2D gaussian kernel that will be convolved with each slice at each time-point. If 0, no spatial smoothing is applied'},...
    {'covCorrection',params.covCorrection,'type=checkbox','Correction for correlated noise'},...
    {'correctionType',correctionTypeMenu,'type=popupmenu','contingent=covCorrection','Type of correction (see Wiki for the different algorithms and computing time comparison)'},...
    {'covEstimation',covEstimationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of Estimation of the noise covariance matrix'},...
    {'covEstimationAreaSize',params.covEstimationAreaSize, 'minmax=[1 inf]','contingent=covCorrection','round=1','dimensions in voxels of the spatial window on which the covariance matrix is estimated (in the X and Y dimensions)'},...
    {'covFactorization',covFactorizationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of factorization of the covariance matrix for the computation of the pre-whitening filter'},...
    {'covEstimationBrainMask',covEstimationBrainMaskMenu,'contingent=covCorrection','ROI mask for the estimation of the noise covariance matrix.'},...
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
  % remove description of HRF model if the type of model has changed
  if ~strcmp(params.hrfModel,tempParams.hrfModel)
    params.hrfParams.description='';
  end
  params = copyFields(tempParams,params);
  
  %perform some checks
  if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None') && ~ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
    mrWarnDlg('(glmAnalysisGUI) Noise covariance estimation mask can only be used with ROI(s) analysis');
  elseif 0
    %perform check on parameters here if needed
  else
    while askForParams       % get hrf model parameters
      %first get some parameters that will be used to plot the HRF model
      %if scans have been chosen already, take the first chosen scan, otherwise take the first scan of the group
      if ~fieldIsNotDefined(params,'scanNum')
        scanNum = params.scanNum(1);
      else
        scanNum=1;
      end
      %we assume that all (chosen) scans in this group have the same framePeriod
      framePeriod = viewGet(thisView,'framePeriod',scanNum,viewGet(thisView,'groupNum',params.groupName));
      %and that the acquisition time is half the frame period
      hrfParams = feval(params.hrfModel,params.hrfParams,framePeriod,framePeriod/2,defaultParams);
      % if empty user hit cancel, go back
      if isempty(hrfParams)
        if size(hrfParams,2) %if the top close button has been pressed
          params=[];
          return
        else
          break;
        end
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
          scanNums = selectInList(thisView,'scans','Select scans to analyse',params.scanNum,groupNum);
          if isempty(scanNums)
            if size(scanNums,2) %if the top close button has been pressed
              params=[];
              return
            else
              askForParams = 1;
              break;
            end
          else
            %check that scan frame Period are all identical
            framePeriod = viewGet(thisView,'framePeriod',scanNums(1));
            for iScan = scanNums(2:end);
              if viewGet(thisView,'framePeriod',iScan)~=framePeriod
                mrWarnDlg('(glmAnalysisGUI) GLM analysis cannot be performed on scans with different frame periods. Copy scans in different groups.');
                params=[];
                return
              end
            end

            params.scanNum = scanNums;
          end
        end

        while askForParams    % get the parameters for each scan
          [scanParams, params] = getGlmScanParamsGUI(thisView,params,defaultParams);
          if isempty(scanParams)
            if size(scanParams,2) %if the top close button has been pressed
              params=[];
              return
            else
              askForParams = 1;
              break;
            end
          end
          params.scanParams = scanParams;


          while askForParams    %get the stim to EV matrices for each scan
            [scanParams, params] = getGlmEVParamsGUI(thisView,params,defaultParams);
            if isempty(scanParams)
              if size(scanParams,2) %if the top close button has been pressed
                params=[];
                return
              else
                askForParams = 1;
                break;
              end
            end
            params.scanParams = scanParams;
            if ~params.numberContrasts && ~params.numberFtests
              defaultParams =1;
            end

            while askForParams           %get params
              tempParams = getGlmTestParamsGUI(thisView,params,defaultParams);
              % if empty, user hit cancel, go back
              if isempty(tempParams)
                if size(tempParams,2) %if the top close button has been pressed
                  params=[];
                  return
                else
                  askForParams = 1;
                  break;
                end
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

