% mlrAnatDBPush.m
%
%        $Id:$ 
%      usage: mlrAnatDBPush(subjectID)
%         by: justin gardner
%       date: 07/05/15
%    purpose: Push the subjectID repo. 
%
function retval = mlrAnatDBPush(subjectID)

% check arguments
if ~any(nargin == [1])
  help mlrAnatDBPush
  return
end

% get the subjectID
subjectID = mlrAnatDBSubjectID(subjectID);
if isempty(subjectID),return,end

% find out if the repos exist
[localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID,'noPull=1');

% keep current path
curpwd = pwd;

% push them if they exist
if ~isempty(localRepo)
  disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s',localRepo));
  cd(localRepo)
  mysystem(sprintf('hg push --new-branch'));
  cd(curpwd);
  disppercent(inf);
end

% push them if they exist
if ~isempty(localRepoLargeFiles)
  disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s. This may take a few minutes',localRepoLargeFiles));
  cd(localRepoLargeFiles)
  mysystem(sprintf('hg push --new-branch'));
  cd(curpwd);
  disppercent(inf);
end


%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPut): %s',command));
[status,result] = system(command,'-echo');

