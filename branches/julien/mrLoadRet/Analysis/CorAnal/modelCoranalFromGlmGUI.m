% modelCoranalFromGlmGUI
%
%        $Id$
%      usage: params = modelCoranalFromGlmGUI('thisView',thisView,'groupName',groupName,'params',params)
%         by: julien besle
%       date: 2010-10-24
%     inputs: 
%    outputs: 
% 
%    purpose: parameters for modelCoranalFromGlm.m
%        e.g:
%

function params = modelCoranalFromGlmGUI(varargin)

% get the arguments
eval(evalargs(varargin));
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;

if ieNotDefined('params')
  params = struct;
  if isempty(thisView.analyses) || ...
      ~ismember(thisView.analyses{thisView.curAnalysis}.type,{'glmAnalStats','glmAnal','erAnal','deconvAnal'})
    mrWarnDlg('You must load an event-related/GLM Analysis');
    params = [];
    return;
  end
  
  analysisParams = viewGet(thisView,'analysisParams');
  params.originalAnalysis = analysisParams.saveName;
  params.originalGroup = analysisParams.groupName;

  d = thisView.analyses{thisView.curAnalysis}.d{1};
  if ~isfield(params,'newGroupName')
    params.newGroupName = ['Average modeled from ' groupName ' GLM'];
  end
  if ~isfield(params,'saveName')
    params.saveName = 'corAnalFromGLM';
  end

  if ~isfield(params,'ncycles')
    params.ncycles = 12;
  end
  if ~isfield(params,'stimDurationS')
    params.stimDurationS = 2*d.tr;
  end
  if ~isfield(params,'cycleLengthTR')
    params.cycleLengthTR = 10;
  end
  if ~isfield(params,'dummyStimsTR')
    params.dummyStimsTR = 6;
  end
  if ~isfield(params,'contrasts')
    params.contrasts = diag(ones(d.nhdr,1));
    %shift the order of contrasts to account for dummies
    params.contrasts = circshift(params.contrasts,round(params.dummyStimsTR/params.cycleLengthTR*size(params.contrasts,1)));
    params.contrasts = mat2str(params.contrasts);
  end
  if ~isfield(params,'nScans')
    %this is the maximum number of scans we can have
    params.nScans = floor(size(d.scm,1)/(params.ncycles*params.cycleLengthTR+params.dummyStimsTR));
  end
  if ~isfield(params,'reverseScans')
    params.reverseScans = mat2str(2:2:params.nScans);
  end
  if ~isfield(params,'delayTR')
    params.delayTR = -2;
  end
  if ~isfield(params,'subsetBox')
    params.subsetBox = mat2str([1 size(d.ehdr,1);1 size(d.ehdr,2);1 size(d.ehdr,3)]);
  end
  if ~isfield(params,'TRCorrection')
    params.TRCorrection = 'Yes';
  end
  if ~isfield(params,'noiseMode')
    params.noiseMode = 'Residuals';
  end
  if ~isfield(params,'detrend')
    params.detrend = 'Highpass';
  end
  if ~isfield(params,'spatialnorm')
    params.spatialnorm = 'None';
  end
  if ~isfield(params,'trigonometricFunction')
    params.trigonometricFunction = 'Cosine';
  end
elseif isfield(params,'scanParams')
  params = params.scanParams{viewGet(thisView,'curScan')};
  %if params is not empty, we assume that they are all present and correct
  if ~isfield(params,'trigonometricFunction') %except for this one, which has been added later
    params.spatialnorm = 'Cosine';
  end
end

groupNames = putOnTopOfList(params.originalGroup,viewGet(thisView,'groupNames'));
TRCorrectionMenu = putOnTopOfList(params.TRCorrection,{'Yes','No'});
noiseModeMenu = putOnTopOfList(params.noiseMode,{'Residuals','No noise'});
detrendMenu = putOnTopOfList(params.detrend,{'Highpass','None','Linear','Quadratic'});
spatialnormMenu = putOnTopOfList(params.spatialnorm,{'Divide by mean','None'});
trigonometricFunctionMenu = putOnTopOfList(params.trigonometricFunction,{'Sine','Cosine'});

askForParams = 1;
while askForParams
  paramsInfo = {...
      {'originalGroup',groupNames,'type=popupmenu','Name of the Group from which the GLM anlaysis is taken'},...
      {'originalAnalysis',params.originalAnalysis,'Name of the GLM analysis from which the estimates are taken'},...
      {'newGroupName',params.newGroupName,'Name of the group in which the created scans/analysis will be install'},...
      {'saveName',params.saveName,'Name of the analysis to save. If identical to an existing analysis a scan will be added and overlays will be added to the existed analysis'},...
      {'contrasts',params.contrasts,'linear combination of estimates to model the different stimulations'},...
      {'nScans',params.nScans,'to do'},...
      {'reverseScans', params.reverseScans, 'to do'},...
      {'ncycles', params.ncycles, 'to do'},...
      {'cycleLengthTR', params.cycleLengthTR, 'to do'},...
      {'stimDurationS', params.stimDurationS, 'to do'},...
      {'dummyStimsTR', params.dummyStimsTR,  'to do'},...
      {'delayTR', params.delayTR, 'to do'},...
      {'TRCorrection', TRCorrectionMenu, 'type=popupmenu', 'to do'},...
      {'noiseMode', noiseModeMenu, 'type=popupmenu', 'to do'},...
      {'detrend', detrendMenu, 'type=popupmenu', 'to do'},...
      {'spatialnorm', params.spatialnorm, 'editable=0', 'Spatial normalization is disabled because time-series mean might be negative'},...
      {'trigonometricFunction', trigonometricFunctionMenu, 'type=popupmenu', 'to do'},...
      {'subsetBox', params.subsetBox, 'subset of voxels, of the form [X1 X2;Y1 Y2;Z1 Z2] (Zs are optional)'},...
    };
%      {'spatialnorm', spatialnormMenu, 'type=popupmenu','editable=0', 'Spatial normalization is disabled because time-series mean might be negative'},...

  params = mrParamsDialog(paramsInfo, 'Correlation from GLM parameters');
  % Abort if params empty
  if ieNotDefined('params'),return,end

  groupNames = putOnTopOfList(params.originalGroup,groupNames);
  TRCorrectionMenu = putOnTopOfList(params.TRCorrection,TRCorrectionMenu);
  noiseModeMenu = putOnTopOfList(params.noiseMode,noiseModeMenu);
  detrendMenu = putOnTopOfList(params.detrend,detrendMenu);
  spatialnormMenu = putOnTopOfList(params.spatialnorm,spatialnormMenu);
  trigonometricFunctionMenu = putOnTopOfList(params.trigonometricFunction,trigonometricFunctionMenu);

  askForParams=0;
  if ~ismember(params.originalAnalysis,viewGet(thisView,'analysisNames',viewGet(thisView,'groupNum',groupNames{1})))
    mrWarnDlg(['(modelCoranalFromGlmGUI) There is no analysis ' params.originalAnalysis ' in group ' groupNames{1}] );
    askForParams=1;
  end
  if max(eval(params.reverseScans))>params.nScans
    mrWarnDlg('(modelCoranalFromGlmGUI) Max(reverseScan) must be less than number of scans');
    askForParams=1;
  end
  if params.nScans/length(eval(params.reverseScans))~=2
    mrWarnDlg('(modelCoranalFromGlmGUI) The number of reverse scans must equal half the total number of scans');
    askForParams=1;
  end
  if ~ieNotDefined('d')
    if strcmp(params.noiseMode,'Residuals') && params.nScans*(params.ncycles*params.cycleLengthTR+params.dummyStimsTR) >size(d.scm,1)
      mrWarnDlg('(modelCoranalFromGlmGUI) If residuals are used as noise, the total number of simulated TRs cannot be greater than the number of TRs in the GLm analysis');
      askForParams=1;
    end
    if rem(params.stimDurationS,d.tr) || rem(params.cycleLengthTR,size(eval(params.contrasts),1))
      mrWarnDlg('(modelCoranalFromGlmGUI) Function not implemented for stimulus durations that are not multiples of TR or cycle lengths that are not multiples of the number of stimuli/contrasts');
      askForParams=1;
    end
  end
            
end
