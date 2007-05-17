% fixBadChars.m
%
%      usage: str = fixBadChars(str)
%         by: justin gardner
%       date: 04/20/07
%    purpose: takes a string and replaces bad characters not
%    allowed in variable names like space or * with variable name acceptable characters
%
function str = fixBadChars(str)

% check arguments
if ~any(nargin == [1])
  help fixBadChars
  return
end

% this is the list of what characters will map to what
swapchars = {{'-','_'},{' ','_'},{'*','star'},{'+','plus'},{'%','percent'},{'[',''},{']',''},{'(',''},{')',''},{'/','_div_'},{'=','_eq_'},{'^','_pow_'},{'.','_period_'},{':','_'},{'&','_and_'},{'!','_bang_'},{'#','_hash_'},{'$','_dollar_'},{'{',''},{'}',''},{'|','_bar_'},{'\','_backslash_'},{';','_'},{'?','_question_'},{',','_comma_'},{'<','_less_'},{'>','_greater_'},{'~','_tilde_'},{'`','_backtick_'}};

% now swap any occurrences of these characters
for i = 1:length(swapchars)
  % look for where we have a bad character
  swaplocs = strfind(str,swapchars{i}{1});
  % if any found replace them
  if ~isempty(swaplocs)
    newstr = '';
    swaplocs = [0 swaplocs];
    for j = 2:length(swaplocs)
      newstr = sprintf('%s%s%s',newstr,str(swaplocs(j-1)+1:swaplocs(j)-1),swapchars{i}{2});
    end
    str = sprintf('%s%s',newstr,str(swaplocs(end)+1:end));
  end
end

% check for non character at beginning
if ~isempty(regexp(str,'^[^a-zA-Z]'))
  str = sprintf('x%s',str);
end
