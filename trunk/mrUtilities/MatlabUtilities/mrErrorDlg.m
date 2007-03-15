function mrErrorDlg(errstr)
%
% function mrErrorDlg(errstr)
%
% Uses 'verbose' preference to either open a matlab errordlg or just signal
% an error.
%
% To set the 'verbose' preference:
%    setpref('mrLoadRet','verbose',0);
%    setpref('mrLoadRet','verbose',1);
%
% djh, 5/2005
%	$Id$	



if ispref('mrLoadRet','verbose')
    verbose = getpref('mrLoadRet','verbose');
else
    verbose = 0;
end

if verbose
  errordlg(errstr,'Error!');
  error(errstr);
else
    disp(errstr);
    keyboard
end
