% mrParamsClose.m
%
%      usage: mrParamsClose()
%         by: justin gardner
%       date: 09/04/11
%    purpose: Closes an open mrParamsDialog
%
function retval = mrParamsClose()

% check arguments
if ~any(nargin == [0])
  help mrParamsClose
  return
end

global gParams;
if isfield(gParams,'fignum')
  if ishandle(gParams.fignum)
    close(gParams.fignum);
  end
end

