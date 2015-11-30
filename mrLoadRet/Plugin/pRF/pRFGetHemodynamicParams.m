% pRFGetHemodynamicParams.m
%
%        $Id:$ 
%      usage: hemoParams = pRFGetHemodynamicParams(fitParams)
%         by: justin gardner
%       date: 11/25/15
%    purpose: Get initial hemodynamic parameters for fitParams
%
function hemo = pRFGetHemodynamicParams(fitParams)

% check arguments
if ~any(nargin == [1])
  help pRFGetHemodynamicParams
  return
end

% parameters for single gamma
hemo.names = {'timelag','tau'};
hemo.initParams = [fitParams.timelag fitParams.tau];
hemp.descriptions = {'Time before start of rise of hemodynamic function','Width of the hemodynamic function (tau parameter of gamma)'};
hemo.incDec = [0.1 0.5];
hemo.min = [0 0];
hemo.max = [inf inf];

% parameters for difference of gammas
if fitParams.diffOfGamma
  hemo.names = {hemo.names{:} 'amp2' 'timelag2','tau2'};
  hemo.descriptions = {hemo.paramDescriptions{:} 'Amplitude of second gamma for HDR' 'Timelag for second gamma for HDR','tau for second gamma for HDR'};
  hemo.incDec = [hemo.incDec(:) 0.1 0.1 0.5];
  hemo.min = [hemo.min(:) 0 0 0];
  hemo.max = [hemo.max(:) inf 6 inf];
  % set init params
  hemo.initParams = [hemo.initParams fitParams.amplitudeRatio fitParams.timelag2 fitParams.tau2];
end
  
