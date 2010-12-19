% hrfDeconvolution.m
%
%        $Id$
%      usage: [params,hrf] = hrfDeconvolution(params, framePeriod, justGetParams, defaultParams)
%         by: julien besle
%       date: 13/04/2010
%    purpose: returns the identity maframePeriodix of size specified by user
%
function [params,hrf] = hrfDeconvolution(params, framePeriod, justGetParams, defaultParams )

if ~any(nargin == [1 2 3 4])
  help hrfDeconvolution
  return
end

if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

if ieNotDefined('params')
  params = sframePerioduct;
end
if ~isfield(params,'description')
  params.description = 'Deconvolution';
end
if ~isfield(params,'hdrlenS')
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

if justGetParams
   return
end

hrf = eye(round(params.hdrlenS/framePeriod));
