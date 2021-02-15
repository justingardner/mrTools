% pRF_gaussian.m
%
%        $Id:$ 
%      usage: pRF_gaussianhdr(varargin)
%         by: akshay jagadeesh
%       date: 09/02/16
%    purpose: Template file to create new models
%
%             - Allows you to create a new model type 
%             - Simply add the model type to pRFGUI and create a 
%               file like this one for your model with the appropriate 
%               methods filled in.
%
%     This model template is set up for the standard Gaussian model.

function output = pRF_gaussianhdr(varargin)

if nargin < 2
  disp(sprintf('Not enough arguments'));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%pRF_gaussianhdr('getModelResponse', fitParams, rfModel, hrf, i)%%%%%
%%%%                Called from getModelResidual                   %%%%%
if strcmp(varargin{1}, 'getModelResponse')

  fitParams = varargin{2};
  rfModel = varargin{3};
  hrf = varargin{4};
  p = varargin{5};
  i = varargin{6};

  nFrames = fitParams.concatInfo.runTransition(i,2);
  thisModelResponse = convolveModelWithStimulus(rfModel,fitParams.stim{i},nFrames);

  % and convolve in time.
  thisModelResponse = convolveModelResponseWithHRF(thisModelResponse,hrf);

  % drop junk frames here
  thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);

  % return the calculated model response
  output = thisModelResponse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%pRF_gaussian('setParams', fitParams)%%%%%
%%%%        Called from setParams       %%%%%

elseif strcmp(varargin{1}, 'setParams')

  fitParams = varargin{2};

  fitParams.paramNames = {'x','y','rfWidth', 'timelag', 'tau'};
  fitParams.paramDescriptions = {'RF x position','RF y position','RF width (std of gaussian)', 'Time before start of rise of hemodynamic function', 'Width of the hemodynamic function (tau parameter of gamma)'};
  fitParams.paramIncDec = [1 1 1 0.1 0.5];
  fitParams.paramMin = [-inf -inf 0 0 0];
  fitParams.paramMax = [inf inf inf inf inf];
  % set min/max and init
  fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0 0 0];
  fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf 3 inf];
  fitParams.initParams = [0 0 4 fitParams.timelag fitParams.tau];
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

  % return fitParams with modified values
  output = fitParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0%%%%%
%%%%  pRF_gaussianhdr(command, fitParams, params)%%%%%
%%%%         Called from getFitParams            %%%%% 

elseif strcmp(varargin{1}, 'getFitParams')

  fitParams = varargin{2};
  params = varargin{3};

  p.rfType = fitParams.rfType;

  % Define your parameters here
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
  output = struct(p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
function modelResponse = convolveModelWithStimulus(rfModel,stim,nFrames)

% get number of frames
nStimFrames = size(stim.im,3);

% preallocate memory
modelResponse = zeros(1,nFrames);

for frameNum = 1:nStimFrames
  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
  modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

