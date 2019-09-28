% mlrFixDescriptionInParams.m
%
%      usage: mlrFixDescriptionInParams()
%         by: justin gardner
%       date: 08/12/19
%    purpose: Pulling out function that converts descriptions in params that have [x...x] to make a pretty
%             description string like min - max or [1 3 5] etc.
%
function params = mlrFixDescriptionInParams(params)

% check arguments
if ~any(nargin == [1])
  help mlrFixDescriptionInParams
  params = [];
  return
end

if isfield(params,'description') && ~isempty(strfind(params.description,'[x...x]'))
  swaploc = strfind(params.description,'[x...x]');
  swaploc = swaploc(1);
  % get the group name
  if isfield(params,'groupName')
    groupName = params.groupName;
  else
    groupName = '';
  end
  % scan numbers
  if isfield(params,'scanList')
    scanNums = params.scanList;
  elseif isfield(params,'scanNum')
    scanNums = params.scanNum;
  else
    scanNums = [];
  end
  % convert to scanNames
  if length(scanNums) > 1
    % see if we can break down into smaller chunks of consecutive numbers
    scanNums = sort(scanNums);
    diffScanNums = diff(scanNums);
    % a break in consecutive numbers is one in which diff is more than one
    scanNumRunBreak = [find(diffScanNums > 1) length(scanNums)];
    % default scanNames
    scanNames = '';
    % start the last break at 1
    lastBreak = 1;
    for iBreak = 1:length(scanNumRunBreak)
      thisBreak = scanNumRunBreak(iBreak);
      if (thisBreak-lastBreak) > 0
	% if it is a run then put in the run
	scanNames = sprintf('%s,%i-%i',scanNames,scanNums(lastBreak),scanNums(thisBreak));
      else
	% otherwise it is just a singleton
	scanNames = sprintf('%s,%i',scanNames,scanNums(lastBreak));
      end
      lastBreak = thisBreak+1;
    end
    % remove starting ,
    if length(scanNames)>1
      scanNames = sprintf(' [%s]',scanNames(2:end));
    end
  elseif length(scanNums) == 1
    scanNames = sprintf('%i',scanNums);
  else
    scanNames = '';
  end
  % now make the description string
  if ~isempty(groupName)
    params.description = sprintf('%s%s:%s%s',params.description(1:swaploc-1),groupName,scanNames,params.description(swaploc+7:end));
  else
    params.description = sprintf('%s%s%s',params.description(1:swaploc-1),scanNames,params.description(swaploc+7:end));
  end
end
