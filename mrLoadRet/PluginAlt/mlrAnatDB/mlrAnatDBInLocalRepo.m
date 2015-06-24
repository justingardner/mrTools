% mlrAnatDBInLocalRepo.m
%
%        $Id:$ 
%      usage: tf = mlrAnatDBInLocalRepo(v)
%         by: justin gardner
%       date: 06/23/15
%    purpose: Tells whether the session is in the mlrAnatDB repo or not
%
%             v = newView;
%             tf = mlrAnatDBInLocalRepo(v);
%
function tf = mlrAnatDBInLocalRepo(v)

localRepoTop = mlrReplaceTilde(mrGetPref('mlrAnatDBLocalRepo'));
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));
% if they are not the same, then offer add session as a menu item,
% but nothing else.
if ~strncmp(localRepoTop,homeDir,length(localRepoTop))
  tf = false;
else
  tf = true;
end

