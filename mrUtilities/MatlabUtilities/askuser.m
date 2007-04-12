% askuser.m
%
%      usage: askuser(question)
%         by: justin gardner
%       date: 02/08/06
%    purpose: ask the user something
%
function r = askuser(question,toall)

% check arguments
if ~any(nargin == [1 2])
  help askuser
  return
end
if ieNotDefined('toall'),toall = 0;,end
r = [];
while isempty(r)
  % ask the question
  if toall
    r = input(sprintf('%s (y/n or a for Yes to all)? ',question),'s');
  else
    r = input(sprintf('%s (y/n)? ',question),'s');
  end
  % make sure we got a valid answer
  if (lower(r) == 'n')
    r = 0;
  elseif (lower(r) == 'y')
    r = 1;
  elseif (lower(r) == 'a') & toall
    r = inf;
  else
    r =[];
  end
end


