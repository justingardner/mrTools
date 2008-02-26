function mrWarnDlg(warnstr)
%
% mrWarnDlg(warnstr)
%
% Calls Matlab's warndlg or warning depending on verbose preference
%
% To set the 'verbose' preference:
%    mrSetPref('verbose','Yes');
%    mrSetPref('verbose','No');
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref

verbose = mrGetPref('verbose');

% always display warning on command line
disp(sprintf('Warning: %s',warnstr));
if strcmp(verbose,'Yes')
  warndlg(warnstr);
  drawnow;
end
