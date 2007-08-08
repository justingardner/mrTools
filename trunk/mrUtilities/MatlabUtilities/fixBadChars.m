% fixBadChars.m
%
%      usage: str = fixBadChars(str,<fixList>)
%         by: justin gardner
%       date: 04/20/07
%    purpose: takes a string and replaces bad characters not
%    allowed in variable names like space or * with variable name acceptable characters
%             you can also provide your own fixlist, i.e. pairs of
%             things that are the match and replacement, e.g.
%             fixBadChars('remove *this* and replace with that',{'*this*','that'})
function str = fixBadChars(str,fixList)

% check arguments
if ~any(nargin == [1 2])
  help fixBadChars
  return
end

% this is the list of what characters will map to what
if ieNotDefined('fixList')
  fixList = {{'-','_'},{' ','_'},{'*','star'},{'+','plus'},{'%','percent'},{'[',''},{']',''},{'(',''},{')',''},{'/','_div_'},{'=','_eq_'},{'^','_pow_'},{'.','_period_'},{':','_'},{'&','_and_'},{'!','_bang_'},{'#','_hash_'},{'$','_dollar_'},{'{',''},{'}',''},{'|','_bar_'},{'\','_backslash_'},{';','_'},{'?','_question_'},{',','_comma_'},{'<','_less_'},{'>','_greater_'},{'~','_tilde_'},{'`','_backtick_'}};
  userDefinedFixList = 0;
else
  fixList = cellArray(fixList,2);
  userDefinedFixList = 1;
end

% now swap any occurrences of these characters
for i = 1:length(fixList)
  % look for where we have a bad character
  swaplocs = strfind(str,fixList{i}{1});
  % if any found replace them
  if ~isempty(swaplocs)
    newstr = '';
    swaplocs = [-length(fixList{i}{1})+1 swaplocs];
    for j = 2:length(swaplocs)
      newstr = sprintf('%s%s%s',newstr,str((swaplocs(j-1)+length(fixList{i}{1})):swaplocs(j)-1),fixList{i}{2});
    end
    str = sprintf('%s%s',newstr,str((swaplocs(end)+length(fixList{i}{1})):end));
  end
end

% check for non character at beginning
if ~userDefinedFixList
  if ~isempty(regexp(str,'^[^a-zA-Z]'))
    str = sprintf('x%s',str);
  end
end
