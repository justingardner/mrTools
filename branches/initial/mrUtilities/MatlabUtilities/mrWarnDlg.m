function mrWarnDlg(warnstr)
%
% mrWarnDlg(warnstr)
%
% Calls Matlab's warndlg or warning depending on verbose preference
%
% To set the 'verbose' preference:
%    setpref('mrLoadRet','verbose',0);
%    setpref('mrLoadRet','verbose',1);
%
% djh, 5/2005

if ispref('mrLoadRet','verbose')
    verbose = getpref('mrLoadRet','verbose');
else
    verbose = 0;
end

if verbose
    warndlg(warnstr);
	drawnow;
else
    warning(warnstr);
end
