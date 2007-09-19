% mrParamsDefault.m
%
%      usage: mrParamsDefault(paramsinfo)
%         by: justin gardner
%       date: 03/13/07
%    purpose: returns default params given paramInfo
%             see wiki for details
%
function params = mrParamsDefault(paramInfo)

% check arguments
if ~any(nargin == [1])
  help mrParamsDefault
  return
end

% parse the input parameter string
[vars varinfo] = mrParamsParse(paramInfo);
for i = 1:length(vars)
  % make sure it is not a contingent value that has been shut
  % off, first get value it is contingent on
  if isfield(varinfo{i},'contingent')
    contingentValue = varinfo{varinfo{i}.contingentOn}.value;
    if isstr(contingentValue),contingentValue = str2num(contingentValue);,end
    % if it has been shut down, give the parameter an empty
    % value and continue on
    if isequal(contingentValue,0)
      params.(varinfo{i}.name) = [];
      continue
    end
  end
  % otherwise set the parameter to the default value
  if ~iscell(varinfo{i}.value)
    params.(varinfo{i}.name) = varinfo{i}.value;
  else
    params.(varinfo{i}.name) = varinfo{i}.value{1};
  end
  % change to numeric if need be
  if strcmp(varinfo{i}.type,'numeric') || strcmp(varinfo{i}.type,'checkbox')
    if isstr(params.(varinfo{i}.name))
      params.(varinfo{i}.name) = str2num(params.(varinfo{i}.name));
    end
  end
end
params.paramInfo = paramInfo;

