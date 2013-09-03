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
  for i = 1:length(gParams.fignum)
    if ishandle(gParams.fignum(i))
      close(gParams.fignum(i));
    end
    if ishandle(gParams.fignum(i))
      delete(gParams.fignum(i));
    end
  end
end

