% mrParamsGet.m
%
%        $Id$ 
%      usage: mrParamsGet(<vars>)
%         by: justin gardner
%       date: 09/23/10
%    purpose: get the current settings of the mrParamsDialog
%
function params = mrParamsGet(vars)

% check arguments
if ~any(nargin == [0 1])
  help mrParamsGet
  return
end
  
global gParams;
if nargin == 0,vars = gParams.vars;end
% return the var entries
for i = 1:length(gParams.varinfo)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % for checkboxes, just return 0 or 1
  if strcmp(gParams.varinfo{i}.type,'checkbox')
    if isfield(gParams.varinfo{i},'group')
      for j = 1:length(gParams.varinfo{i}.allValues)
        % if this is the current one then use field val
        if gParams.varinfo{i}.oldControlVal == j
          params.(gParams.varinfo{i}.name)(j) = get(gParams.ui.varentry{i},'Value');
        else
          params.(gParams.varinfo{i}.name)(j) = gParams.varinfo{i}.allValues{j};
        end
      end
    else
      params.(gParams.varinfo{i}.name) = get(gParams.ui.varentry{i},'Value');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for arrays, have to get all values
  elseif ismember(gParams.varinfo{i}.type,{'array' 'stringarray'})
    if ~isfield(gParams.varinfo{i},'group')
      % if not grouped, just get the value from the gui
      for iRows = 1:size(gParams.ui.varentry{i},1)
        for iCols = 1:size(gParams.ui.varentry{i},2)
          value = get(gParams.ui.varentry{i}(iRows,iCols),'String');
          if ~isempty(value)
            if strcmp(gParams.varinfo{i}.type,'array')
              params.(gParams.varinfo{i}.name)(iRows,iCols) = mrStr2num(value);
            else
              params.(gParams.varinfo{i}.name){iRows,iCols} = value;
            end
          end
        end
      end
      % if grouped, either get value from gui or form allvalues
    else
      for j = 1:length(gParams.varinfo{i}.allValues)
        if gParams.varinfo{i}.oldControlVal == j
	  for iRows = 1:size(gParams.ui.varentry{i},1)
	    for iCols = 1:size(gParams.ui.varentry{i},2)
          if strcmp(gParams.varinfo{i}.type,'array')
            params.(gParams.varinfo{i}.name){j}(iRows,iCols) = mrStr2num(get(gParams.ui.varentry{i}(iRows,iCols),'String'));
          else
            params.(gParams.varinfo{i}.name){j}(iRows,iCols) = get(gParams.ui.varentry{i}(iRows,iCols),'String');
          end
	    end
	  end
        else
          params.(gParams.varinfo{i}.name){j} = gParams.varinfo{i}.allValues{j};
        end
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for a push button just return whatever is in the value field
  elseif strcmp(gParams.varinfo{i}.type,'pushbutton')
    params.(gParams.varinfo{i}.name) = gParams.varinfo{i}.value;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for pop up menus, get the value and look it up in the original list
  elseif strcmp(gParams.varinfo{i}.type,'popupmenu')
    % get the current value
    val = get(gParams.ui.varentry{i},'Value');
    list = get(gParams.ui.varentry{i},'String');
    % if this is a group return a cell array
    if isfield(gParams.varinfo{i},'group')
      for j = 1:length(gParams.varinfo{i}.allValues)
        % if this is the current list thenuse current val
        if gParams.varinfo{i}.oldControlVal == j
          params.(gParams.varinfo{i}.name){j} = list{val};
          % else get the list form allValues
        else
          params.(gParams.varinfo{i}.name){j} = gParams.varinfo{i}.allValues{j}{1};
        end
      end
    else
      params.(gParams.varinfo{i}.name) = list{val};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All other values
  else
    if isfield(gParams.varinfo{i},'group')
      for j = 1:length(gParams.varinfo{i}.allValues)
        % if this is the current one then use field val
        if gParams.varinfo{i}.oldControlVal == j
          params.(gParams.varinfo{i}.name){j} = get(gParams.ui.varentry{i},'String');
        else
          params.(gParams.varinfo{i}.name){j} = gParams.varinfo{i}.allValues{j};
        end
      end
    else
      params.(gParams.varinfo{i}.name) = get(gParams.ui.varentry{i},'String');
    end
    % check to see if this is an array that could not be displayed
    if isfield(gParams.varinfo{i},'tooBigArrayValue')
      % if so then set to the original array value
      params.(gParams.varinfo{i}.name) = gParams.varinfo{i}.tooBigArrayValue;
    end
  end
  % change numeric popupmenu to number
  if strcmp(gParams.varinfo{i}.type,'popupmenu') && strcmp(gParams.varinfo{i}.popuptype,'numeric')
    params.(gParams.varinfo{i}.name) = mrStr2num(params.(gParams.varinfo{i}.name));
  end
  % if non numeric then convert back to a number
  if ~any(strcmp(gParams.varinfo{i}.type,{'string' 'popupmenu' 'array' 'stringarray' 'checkbox' 'pushbutton','statictext'}))
    if isfield(gParams.varinfo{i},'group')
      for j = 1:length(gParams.varinfo{i}.allValues)
        % if this is the current one then use field val
        if isempty(params.(gParams.varinfo{i}.name){j})
          temp(j) = nan;
        elseif isstr(params.(gParams.varinfo{i}.name){j})
          temp(j) = mrStr2num(params.(gParams.varinfo{i}.name){j});
        else
          temp(j) = params.(gParams.varinfo{i}.name){j};
        end
      end
      params.(gParams.varinfo{i}.name) = temp;
    else
      params.(gParams.varinfo{i}.name) = mrStr2num(params.(gParams.varinfo{i}.name));
    end
  end
  % not enabled then set parameter to empty
  if strcmp(get(gParams.ui.varentry{i},'Enable'),'off'); %JB: it would make more sense to leave disabled parameters as they are
    params.(gParams.varinfo{i}.name) = [];                % and not set them to empty
  end                                                     % functions that use these parameters should know what to do
end
params.paramInfo = vars;

