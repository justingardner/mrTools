% hrfDeconvolution.m
%
%        $Id$
%      usage: [params,hrf] = hrfDeconvolution(params, sampleDuration, notUsed, defaultParams)
%         by: julien besle
%       date: 13/04/2010
%    purpose: returns a deconvolution matrix given the design sampling period and the estimation sampling period
%
function [params,hrf] = hrfDeconvolution(params, designSamplingPeriod, ~, defaultParams, estimationSamplingPeriod)

if ~any(nargin == [1 2 3 4 5])
  help hrfDeconvolution
  return
end

%estimationSampling = varargin{1};

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('estimationSamplingPeriod'),estimationSamplingPeriod = designSamplingPeriod;end

if ieNotDefined('params')
  params = struct;
end
if fieldIsNotDefined(params,'description')
  params.description = 'Deconvolution';
end
if fieldIsNotDefined(params,'hdrlenS')
  params.hdrlenS = 25;
end
   
paramsInfo = {...
    {'description', params.description, 'comment describing the hdr model'},...
    {'hdrlenS',params.hdrlenS,'length of the HDR to estimate in seconds'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set Deconvolution parameters');
end

if nargout==1
   return
end

hrf = eye(round(params.hdrlenS/designSamplingPeriod));

% if the estimation sampling period is a multiple of the design sampling period
downsamplingFactor = round(estimationSamplingPeriod/designSamplingPeriod); % downsampling factor from design to estimation sampling rates
hrf = hrf(:,1:downsamplingFactor:end);                                     % (although this is rounded, we generally assume that the estimation sampling period is a multiple of the design sampling period)
