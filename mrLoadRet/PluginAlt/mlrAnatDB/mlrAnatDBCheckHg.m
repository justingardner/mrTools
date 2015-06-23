% mlrAnatDBCheckHg.m
%
%        $Id:$ 
%      usage: tf = lrAnatDBCheckHg()
%         by: justin gardner
%       date: 06/22/15
%    purpose: Returns whether mercurial is correctly installed or not
%
function tf = mlrAnatDBCheckHg

tf = false;

% check for hg
[status,result] = system('which hg');
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) You do not have mercurial installed. You will need to install mercurial. Typicaly by going to the website: http://mercurial.selenic.com and following download instructions.'));
  return
end
% check here for config stuff
[status,result] = system('hg config');
if ~strfind(result,'extensions.largefiles')
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Your hg config needs to have the extension for largefiles enabled. This is done by running "hg config --edit" and adding the line "largefiles =" after the section header "[extensions]"'));
  return
end
[status,result] = system('hg config ui.username');
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Your hg config needs to have your username and email address for you to commit. Use "hg config --edit" to fix this.'));
  return
end
[status,result] = system('hg config web.cacerts');
if status ~= 0
  oneTimeWarning('cacerts',sprintf('(mlrAnatDBPlugin) Your hg config does not have web.cacerts specified - this is useful so that it will allow committing to a self-certified https:// site.'));
end
[status,result] = system('hg config auth');
if status ~= 0
  oneTimeWarning('auth',sprintf('(mlrAnatDBPlugin) Your hg config does not have any authorizations specified. You may want to add auth for the site so that you do not have to keep putting in your username password.'));
end
  
tf = true;

