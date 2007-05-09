function mrErrorDlg(errstr)
%
% function mrErrorDlg(errstr)
%
% Uses 'verbose' preference to either open a matlab errordlg or just signal
% an error.
%
% To set the 'verbose' preference:
%    mrSetPref('verbose',0);
%    mrSetPref('verbose',1);
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref

verbose = mrGetPref('verbose');

if verbose
  errordlg(errstr,'Error!');
  error(errstr);
else
  error(errstr);
end
