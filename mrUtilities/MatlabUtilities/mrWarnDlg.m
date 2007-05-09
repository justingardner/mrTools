function mrWarnDlg(warnstr)
%
% mrWarnDlg(warnstr)
%
% Calls Matlab's warndlg or warning depending on verbose preference
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
    warndlg(warnstr);
	drawnow;
else
    warning(warnstr);
end
