% askuser.m
%
%      usage: askuser(question,<askOption>,<useDialog>)
%         by: justin gardner
%       date: 02/08/06
%    purpose: ask the user a yes/no question. Question is a string or a cell arry of strings with the question. This
%             function will return 1 for yes and 0 for no. If askOption is set to 1, then
%             'Yes to all' will be an option, which if selected will return inf
%
%             If askOption is set to -1 then there will be a default answer of Y
%             If askOption is set to -2 then there will be a default answer of N
%
function r = askuser(question,askOption,useDialog)

% check arguments
if ~any(nargin == [1 2 3])
  help askuser
  return
end

if ieNotDefined('askOption'),askOption = 0;,end
if ieNotDefined('useDialog'),useDialog=0;end

% if explicitly set to useDialog then use dialog, otherwise
% check verbose setting
if useDialog
  verbose = 1;
else
  verbose = mrGetPref('verbose');
  if strcmp(verbose,'Yes'),verbose = 1;else,verbose = 0;end
end

r = [];

question=cellstr(question); %convert question into a cell array
  

while isempty(r)
  % ask the question
  if ~verbose
    % not verbose, use text question
    %fprintf('\n');
    for iLine = 1:length(question)-1
      fprintf([question{iLine} '\n']);
    end
    if askOption==1
      % ask the question (with askOption for all)
      r = input([question{end} ' (y/n or a for Yes to all)? '],'s');
    elseif askOption==-1
      % ask the question (with default to Y)
      r = input([question{end} ' (y/n [default: y])? '],'s');
    elseif askOption==-2
      % ask the question (with default to Y)
      r = input([question{end} ' (y/n [default: n])? '],'s');
    else
      % ask question (without askOption for all)
      r = input([question{end} ' (y/n)? '],'s');
    end
  else
    if askOption==1
      % verbose, use dialog
      r = questdlg(question,'','Yes','No','All','Yes');
      r = lower(r(1));
    else
      r = questdlg(question,'','Yes','No','Yes');
      r = lower(r(1));
    end
  end
  % make sure we got a valid answer
  if (lower(r) == 'n')
    r = 0;
  elseif (lower(r) == 'y')
    r = 1;
  elseif (lower(r) == 'a') & (askOption==1)
    r = inf;
  else
    % if empty then can apply defaults if asked for
    if isempty(r)
      % default answer is yes
      if askOption==-1
	r = 1;
	% default answer is no
      elseif askOption==-2
	r = 0;
      end
    else
      r =[];
    end
  end
end


