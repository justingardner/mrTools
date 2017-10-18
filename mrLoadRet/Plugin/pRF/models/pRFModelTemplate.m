% pRFModelTemplate.m
%
%        $Id:$ 
%      usage: pRFModelTemplate(varargin)
%         by: akshay jagadeesh
%       date: 08/23/16
%    purpose: Template file to create new models
%
%             - Allows you to create a new model type 
%             - Simply add the model type to pRFGUI and create a 
%               file like this one for your model with the appropriate 
%               methods filled in.
%
%     This model template is set up for the standard Gaussian model.

function output = pRFModelTemplate(varargin)

if nargin <=2
  disp(sprintf('Not enough arguments'));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%pRFModelTemplate(command, fitParams, rfModel, hrf, i)%%%%%
%%%%          Called from getModelResidual               %%%%%
if strcmp(varargin{1}, 'getModelResponse')

  fitParams = varargin{2};
  rfModel = varargin{3};
  hrf = varargin{4};
  i = varargin{5};

  nFrames = fitParams.concatInfo.runTransition(i,2);
  thisModelResponse = convolveModelWithStimulus(rfModel,fitParams.stim{i},nFrames);

  % and convolve in time.
  thisModelResponse = convolveModelResponseWithHRF(thisModelResponse,hrf);

  % drop junk frames here
  thisModelResponse = thisModelResponse(fitParams.concatInfo.totalJunkedFrames(i)+1:end);

  % return the calculated model response
  output = thisModelResponse;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%pRFModelTemplate(command, fitParams)%%%%%
%%%%    Called from setFitParams    %%%%%

elseif strcmp(varargin{1}, 'setParams')

  fitParams = varargin{2};

  fitParams.paramNames = {'x','y','rfWidth'};
  fitParams.paramDescriptions = {'RF x position','RF y position','RF width (std of gaussian)'};
  fitParams.paramIncDec = [1 1 1];
  fitParams.paramMin = [-inf -inf 0];
  fitParams.paramMax = [inf inf inf];
  % set min/max and init
  fitParams.minParams = [fitParams.stimExtents(1) fitParams.stimExtents(2) 0];
  fitParams.maxParams = [fitParams.stimExtents(3) fitParams.stimExtents(4) inf];
  fitParams.initParams = [0 0 4];

  % return fitParams with modified values
  output = fitParams;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%pRFModelTemplate(command, fitParams, params)%%%%%
%%%%         Called from getFitParams           %%%%% 

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
  p.canonical.timelag = fitParams.timelag;
  p.canonical.tau = fitParams.tau;
  p.canonical.exponent = fitParams.exponent;
  p.canonical.offset = 0;
  p.canonical.diffOfGamma = fitParams.diffOfGamma;
  p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
  p.canonical.timelag2 = fitParams.timelag2;
  p.canonical.tau2 = fitParams.tau2;
  p.canonical.exponent2 = fitParams.exponent2;
  p.canonical.offset2 = 0;

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

