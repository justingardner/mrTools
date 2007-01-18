% askuser.m
%
%      usage: askuser(question)
%         by: justin gardner
%       date: 02/08/06
%    purpose: ask the user something
%
function r = askuser(question)

% check arguments
if ~any(nargin == [1])
  help askuser
  return
end

r = [];
while isempty(r)
  % ask the question
  r = input(sprintf('%s (y/n)? ',question),'s');
  % make sure we got a valid answer
  if (lower(r) == 'n')
    r = 0;
  elseif (lower(r) == 'y')
    r = 1;
  else
    r =[];
  end
end


