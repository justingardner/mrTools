% hrfDoubleGamma.m
%
%        $Id$
%      usage: [params,hrf] = hrfDoubleGamma(params, designSampling, justGetParams, defaultParams)
%         by: farshad moradi, modified by julien besle
%       date: 14/06/07, 09/02/2010
%    purpose: returns a canonical hrf that's a difference of two gamma distribution function
%
function [params, hrf] = hrfDoubleGamma(params, designSampling, justGetParams, defaultParams)%, varargin )

threshold = 1e-3; %threshold for removing trailing zeros at the end of the model

if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('designSampling'),designSampling = 1;end

if ieNotDefined('params')
  params = struct;
end
if ~isfield(params,'description')
  params.description = 'Double Gamma Function';
end
if ~isfield(params,'x')
  params.x = 6;
end
if ~isfield(params,'y')
  params.y = 16;
end
if ~isfield(params,'z')
  params.z = 6;
end
if ~isfield(params,'includeDerivative')
  params.includeDerivative = 0;
end

% figure handle if one is up
global modelHRFFig;
modelHRFFig = [];

paramsInfo = {...
    {'description', params.description, 'Comment describing the hdr model'},...
    {'x', params.x, 'Shape parameter of the positive Gamma distribution function; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'y', params.y, 'Shape parameter of the negative Gamma distribution function; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'z', params.z, 'Scaling factor between the positive and negative gamma componenets; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'includeDerivative',params.includeDerivative,'type=checkbox','Includes derivative of the hrf in the model'},...
    {'displayHRF',0,'type=pushbutton','callback',@dispModelHRF,'buttonString=Display HRF','passParams=1','callbackArg',{threshold designSampling},'Display the hrf with current parameters'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set model HRF parameters');
end

% close figure that displays HRF if it is up
if ~isempty(modelHRFFig)
  close(modelHRFFig);
end

if justGetParams
   return
end

% get the model HRF
if ~isempty(params)
  [params hrf] = getModelHrf(params,threshold,designSampling);
end


%%%%%%%%%%%%%%%%%%%%%
%    getModelHrf    %
%%%%%%%%%%%%%%%%%%%%%
function [params hrf t] = getModelHrf(params,threshold,designSampling)

tmax = max(params.y*3, 20); %min length of the hrf model in seconds

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.05;

t = 0:dt:tmax;
warning('off', 'MATLAB:log:logOfZero');
modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
warning('on', 'MATLAB:log:logOfZero');

if shift<0
  modelHrf = [zeros(1, ceil(-shift/dt)), modelHrf];
elseif shift>0
  modelHrf = modelHrf( ceil(shift/dt):end );
end


if params.includeDerivative
  % take the derivative
  modelHrfDerivative = [diff(modelHrf), 0];
  % orthogonalize
  modelHrfDerivative = modelHrfDerivative - modelHrf*(modelHrfDerivative/modelHrf);
  % remove mean
  modelHrfDerivative = modelHrfDerivative - mean(modelHrfDerivative);
  % normalize so that it's norm equals the Hrf norm
  modelHrfDerivative = modelHrfDerivative / norm(modelHrfDerivative)*norm(modelHrf);
  %concatenate
  modelHrf = [modelHrf; modelHrfDerivative];
end

%remove trailing zeros
modelHrf = modelHrf(1:end-find(flipud(max(abs(modelHrf),[],2))>threshold,1,'first')+1,:);
%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));
    
%downsample with constant integral
hrf = downsample(modelHrf', round(designSampling/dt));


params.maxModelHrf = designSampling/dt * max(modelHrf'); %output the max amplitude of the actual model HRF

% return actual time
t = t(1:round(designSampling/dt):length(t));
% make sure t is same length as hrf
if length(t) < length(hrf)
  t = [t nan(1,length(hrf)-length(t))];
elseif length(t) > length(hrf)
  t = t(1:length(hrf));
end

%%%%%%%%%%%%%%%%%%%%%%
%    dispModelHRF    %
%%%%%%%%%%%%%%%%%%%%%%
function retval = dispModelHRF(callbackArg,params)

global modelHRFFig;

[params hrf t] = getModelHrf(params,callbackArg{1},callbackArg{2});
modelHRFFig = mlrSmartfig('hrfDoubleGamma','reuse');clf;
plot(t,hrf);
title('Model HRF');
xlabel('Time (sec)');
ylabel('Response magnitude');
