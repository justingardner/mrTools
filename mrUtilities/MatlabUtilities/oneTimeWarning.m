% oneTimeWarning.m
%
%        $Id$
%      usage: oneTimeWarning(warningName,warningText)
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
%
function oneTimeWarning(fieldCheck,warnText)

% check arguments
if ~any(nargin == [2])
  help oneTimeWarning
  return
end

global gMLRWarning
verbose = mrGetPref('verbose');
fieldCheck = fixBadChars(fieldCheck);
if ~isfield(gMLRWarning,fieldCheck)
  gMLRWarning.(fieldCheck) = 1;
  mrWarnDlg(warnText)
end



