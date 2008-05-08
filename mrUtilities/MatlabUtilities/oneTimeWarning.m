% oneTimeWarning.m
%
%        $Id$
%      usage: oneTimeWarning(warningName,warningText,<justPrint>)
%         by: justin gardner
%       date: 02/26/08
%    purpose: Brought over from Shani's viewGet function. This
%             will only warn once about something, you give it
%             a warningName and a warningText e.g.:
%
%             oneTimeWarning('someWarning','Something has happened');
%
%             The first time this happens it prints out the
%             warning, after that it won't do anything.
%
%             If the warning is actually just a message, and
%             not a warning (i.e. you just want it printed
%             as is to the matlab buffer), set the justPrint
%             argument to 1
%

function oneTimeWarning(fieldCheck,warnText,justPrint)

% check arguments
if ~any(nargin == [2 3])
  help oneTimeWarning
  return
end

global gMLRWarning
verbose = mrGetPref('verbose');
fieldCheck = fixBadChars(fieldCheck);
if ~isfield(gMLRWarning,fieldCheck)
  gMLRWarning.(fieldCheck) = 1;
  if ieNotDefined('justPrint')
    mrWarnDlg(warnText)
  else
    disp(warnText);
  end
end



