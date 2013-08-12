% mrParamsEnable.m
%
%        $Id$ 
%      usage: mrParamsEnable(paramname,enable)
%         by: justin gardner
%       date: 07/23/09
%    purpose: enable or disable a a parameter in an mrParamsDialog by name. Enable should be 1 or 0
%
function retval = mrParamsEnable(paramname,enable)

% check arguments
if ~any(nargin == [2])
  help mrParamsEnable
  return
end

if isequal(enable,1)
  enable = 'on';
elseif isequal(enable,0)
  enable = 'off';
end

global gParams;
if isempty(gParams)
  disp('(mrParamsEnable) mrParamsDialog is not running');
  return
end

for i = 1:length(gParams.varinfo)
  if strcmp(gParams.varinfo{i}.name,paramname)
    for j = 1:length(gParams.ui.varentry{i})
      set(gParams.ui.varentry{i}(j),'enable',enable);
    end
    return
  end
end

disp(sprintf('(mrParamsEnable) Could not find parameter %s',paramname));


