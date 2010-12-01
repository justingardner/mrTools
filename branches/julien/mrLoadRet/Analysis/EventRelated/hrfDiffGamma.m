% hrfDiffGamma.m
%
%        $Id$
%      usage: [params,hrf] = hrfDiffGamma(params, tr, stimDuration )
%         by: farshad moradi, modified by julien besle
%       date: 14/06/07, 09/02/2010
%    purpose: returns a canonical hrf
%
function [params, hrf] = hrfDiffGamma(params, tr, justGetParams, defaultParams )

if ~any(nargin == [1 2 3 4])
  help hrfDiffGamma
  return
end

threshold = 1e-3; %threshold for removing trailing zeros at the end of the model

if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

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
  
paramsInfo = {...
    {'description', params.description, 'comment describing the hdr model'},...
    {'x', params.x, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'y', params.y, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'z', params.z, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'includeDerivative',params.includeDerivative,'type=checkbox','include derivative of the hrf in the model?'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set model HRF parameters');
end

if justGetParams
   return
end

tmax = max(params.y*3, 20); %min length of the hrf model in seconds

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.01;

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

%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));
    
%downsample with constant integral
hrf = downsample(modelHrf', round(tr/dt));
%remove trailing zeros
hrf = hrf(1:end-find(flipud(max(abs(hrf),[],2))>threshold,1,'first')+1,:);


params.maxModelHrf = tr/dt * max(modelHrf'); %output the max amplitude of the actual model HRF

