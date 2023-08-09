% hrfDeconvolution.m
%
%        $Id$
%      usage: [params,hrf] = hrfCustom(params, sampleDuration, sampleDelay, defaultParams)
%         by: julien besle
%       date: 13/04/2010
%    purpose: returns the HRF specified as values in the parameters. If a time vector is specified
%             checks that times correspond to actual TR and acquistion time (sampleDuration and sampleDelay)
%             otherwise, assumes that HRF sample times correspond to those TR and acquisition times
%
function [params,hrf] = hrfCustom(params, sampleDuration, sampleDelay, defaultParams, ~)

if ~any(nargin == [1 2 3 4])% 5])
  help hrfCustom
  return
end

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('sampleDelay')
  sampleDelay=sampleDuration/2;
end

if ieNotDefined('params')
  params = struct;
end
if fieldIsNotDefined(params,'description')
  params.description = 'Custom HRF';
end
if fieldIsNotDefined(params,'hrf')
  [~, params.hrf] = hrfDoubleGamma([],sampleDuration,sampleDelay,1);
  params.hrf = params.hrf';
end
if fieldIsNotDefined(params,'hrfTimes')
  params.hrfTimes = sampleDelay+sampleDuration*(0:length(params.hrf)-1);
end
   
paramsInfo = {...
    {'description', params.description, 'comment describing the hdr model'},...
    {'hrf',params.hrf,'values of the the HRF'},...
    {'hrfTimes',params.hrfTimes,'Times of the HRF samples'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set Custom HRF parameters');
end

if nargout==1
   return
end

%check that the times correspond to 
if ~isequal(sampleDelay+sampleDuration*(0:length(params.hrf)-1), params.hrfTimes)
  mrWarnDlg('(hrfCustom) HRF times are not compatible with TR and acquisition time'); 
  keyoard
else
  if size(params.hrf,1)==1
    hrf = params.hrf';
  else
    hrf = params.hrf;
  end
end