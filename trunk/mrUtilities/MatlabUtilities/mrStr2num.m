% mrStr2num
%
%        $Id$
%      usage: num = mrStr2n(str)
%         by: justin gardner
%       date: 07/17/07
%    purpose: returns number or empty if string is not a number
%             matlab's str2num is *very* annoying since it
%             evaluates strings, so that if your string happens to
%             have, say a function name in it then that function
%             will be called, instead of just returning []
%
function retval = mrStr2num(str)

% check arguments
if ~any(nargin == [1])
  help mrStr2num
  return
end

% remove from the string any nan/inf for testing
% since these are valid strings to have.
teststr = fixBadChars(str,{{'nan',''},{'inf',''}});

% check if the string is a valid function or if it has
% any characters in it
if any(exist(teststr) == [2 3 5]) || any(regexp(teststr,'[a-zA-Z]')) 
  % then return empty
  retval = [];
else
  retval = str2num(str);
end
