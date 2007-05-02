% mrParamsParse.m
%
%      usage: mrParamsParse(vars)
%         by: justin gardner
%       date: 03/13/07
%    purpose: called by mrDefaultParamsGUI and mrDefaultParamsReconcile
%             parses the variable information, see wiki for details on
%             format of vairable string
%
%	$Id$	

function [vars varinfo nrows ncols] = mrParamsParse(vars)

% check arguments
if ~any(nargin == [1])
  help mrParamsParse
  return
end

nrows = 1;
ncols = 4;
% first parse the argument
for i = 1:length(vars)
  nrows = nrows+1;
  % if the variable is just a string, then
  % it got passed in without a default argument
  % so make it into a cell array with the second
  % element set to empty
  if ~iscell(vars{i})
    thisvars = vars{i};
    vars{i} = {};
    vars{i}{1} = thisvars;
    vars{i}{2} = '0';
    varinfo{i}.type = 'numeric';
  % no default argument
  elseif length(vars{i}) == 1
    vars{i}{2} = '0';
    varinfo{i}.type = 'numeric';
  % default arguments have to be strings so they
  % can be put into the text fields properly. here
  % we change them into strings, but remember what
  % type they were
  elseif length(vars{i}) >= 2
    if isnumeric(vars{i}{2})
      % check to see if it is an array
      if isscalar(vars{i}{2})
	vars{i}{2} = num2str(vars{i}{2});
	varinfo{i}.type = 'numeric';
      else
	varinfo{i}.type = 'array';
	nrows = nrows+size(vars{i}{2},1)-1;
	ncols = max(ncols,1+size(vars{i}{2},2));
      end
    elseif iscell(vars{i}{2})
      varinfo{i}.type = 'popupmenu';
    else
      varinfo{i}.type = 'string';
    end
  end
  % set info in varinfo
  varinfo{i}.name = vars{i}{1};
  varinfo{i}.value = vars{i}{2};
  varinfo{i}.description = '';
  % check for options
  if length(vars{i}) > 2
    for j = 3:length(vars{i})
      % if this looks like a description then save it as a
      % description
      if (length(strfind(vars{i}{j},'=')) ~= 1) || ...
	    ~isempty(strfind(vars{i}{j}(1:strfind(vars{i}{j},'=')),' '))
	varinfo{i}.description = vars{i}{j};
      else
	% we are going to call evalargs but we want the variables
	% set as a part of gParams (also do it quietly)--> that is
	% evalargs, will do the parsing of the variable=value strings
	setparam{1} = 'gVerbose = 0';
	setparam{2} = sprintf('varinfo{i}.%s',vars{i}{j});
        % set the argument
	eval(evalargs(setparam));
      end
    end
  end
  % make sure type is in lower case
  varinfo{i}.type = lower(varinfo{i}.type);
  % check for minmax violation
  if strcmp(varinfo{i}.type,'numeric') && isfield(varinfo{i},'minmax')
    if vars{i}{2} < varinfo{i}.minmax(1)
      vars{i}{2} = varinfo{i}.minmax(1);
    elseif vars{i}{2} > varinfo{i}.minmax(2)
      vars{i}{2} = varinfo{i}.minmax(2);
    end
  end
  % make editable by default
  if ~isfield(varinfo{i},'editable')
    varinfo{i}.editable = 1;
  end
end

% now check any contingencies
for i = 1:length(varinfo)
  if isfield(varinfo{i},'contingent')
    % go look for the control variable and set it to 
    % have a pointer to this variable
    foundControlVariable = 0;
    for j = 1:length(varinfo)
      if strcmp(varinfo{i}.contingent,varinfo{j}.name)
	foundControlVariable = 1;
	varinfo{i}.contingentOn = j;
	if ~isfield(varinfo{j},'controls')
	  varinfo{j}.controls = i;
	else
	  varinfo{j}.controls(end+1) = i;
	end
      end
    end
    if ~foundControlVariable
      disp(sprintf('Control variable for %s (%s) not found, ignoring contingency',varinfo{i}.name,varinfo{i}.contingent));
    end
  end
end

