% mrParamsReconcile.m
%
%      usage: mrParamsReconcile(params)
%         by: justin gardner
%       date: 03/13/07
%    purpose: A default params reconcile. This does not necessarily
%             have to be called for params created by mrDefaultParamsGUI,
%             but if it is, then there will be a paramInfo field (for format
%             see wiki), and variables will be validated against that structure
%
function [params data] = mrParamsReconcile(groupName,params,data)

% check arguments
if ~any(nargin == [1 2 3])
  help mrParamsReconcile
  return
end

% if this has the field paramInfo then it comes from'
% mrDefaultParamsGUI and we can check the parameters 
% for consistency
if isfield(params,'paramInfo')
  [vars varinfo] = mrParamsParse(params.paramInfo);
  for i = 1:length(varinfo)
    useDefault = 0;
    % make sure the field exists
    if ~isfield(params,varinfo{i}.name)
      useDefault = 1;
    % check if it is good
    else
      % is it numeric
      if strcmp(lower(varinfo{i}.type),'numeric')
	if ~isnumeric(params.(varinfo{i}.name))
	  useDefault = 1;
	end
	% also check if there is a minmax
	if isfield(varinfo{i},'minmax')
	  if params.(varinfo{i}.name) < varinfo{i}.minmax(1)
	    disp(sprintf('(mrParamsReconcile) %s out of range. Setting to min=%f',varinfo{i}.name,varinfo.minmax(1)));
	    params.(varinfo{i}.name) = varinfo{i}.minmax(1);
	  elseif params.(varinfo{i}.name) > varinfo{i}.minmax(2)
	    disp(sprintf('(mrParamsReconcile) %s out of range. Setting to max=%f',varinfo{i}.name,varinfo.minmax(2)));
	    params.(varinfo{i}.name) = varinfo{i}.minmax(2);
	  end
	end
      % or a checkbox
      elseif strcmp(lower(varinfo{i}.type),'checkbox')
	if ~any(params.(varinfo{i}.name) == [0 1])
	  useDefault = 1;
	end
      % or a list of items
      elseif strcmp(lower(varinfo{i}.type),'popupmenu')
	if ~any(strcmp(params.(varinfo{i}.name),varinfo{i}.value))
	  useDefault = 1;
	end
      end
    end      
    % see if we have to switch it to default
    if useDefault
      if ~iscell(varinfo{i}.value)
	params.(varinfo{i}.name) = varinfo{i}.value;
      else
	params.(varinfo{i}.name) = varinfo{i}.value{1};
      end
      if isnumeric(params.(varinfo{i}.name))
	disp(sprintf('(mrParamsReconcile) Using default value for %s (%s)',varinfo{i}.name,num2str(params.(varinfo{i}.name))));
      else
	disp(sprintf('(mrParamsReconcile) Using default value for %s (%s)',varinfo{i}.name,params.(varinfo{i}.name)));
      end
    end
  end
end

% look for a description, and see if it has something like [x...x], that should be replaced by the scan numbers selected in scanList and groupName 
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
    scanNames = num2str(params.scanList);
  elseif isfield(params,'scanNum')
    scanNames = num2str(params.scanNum);
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

% get scan numbers
if isfield(params,'scanList')
  scanNums = params.scanList;
elseif isfield(params,'scanNum')
  scanNums = params.scanNum;
else
  disp('(mrParamsReconcile) Could not find scan numbers');
  return
end
% get group number
if isfield(params,'groupName')
  groupName = params.groupName;
end
groupNum = viewGet([],'groupNum',groupName);

  % get the tseries name
if ~isfield(params,'tseriesFile')
  for iscan = 1:length(scanNums)
    params.tseriesFile{iscan} = viewGet([],'tseriesFile',scanNums(iscan),groupNum);
  end
else
  % see if the tseriesFile name matches
  for iscan = 1:length(scanNums)
    if ~strcmp(params.tseriesFile{iscan},viewGet([],'tseriesFile',scanNums(iscan),groupNum));
      disp(sprintf('(mrParamsReconcile) Scan %i has filename %s which does not match previous one %s',scanNums(iscan),viewGet([],'tseriesFile',scanNums(iscan),groupNum),params.tseriesFile{iscan}));
    end
  end
end  
