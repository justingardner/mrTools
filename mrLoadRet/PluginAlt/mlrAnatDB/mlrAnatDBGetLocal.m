% mlrAnatDBGetLocal.m
%
%        $Id:$ 
%      usage: [tf localRepoROI localRepoSession] = mlrAnatDBGetLocal(subjectID,ROIonly)
%         by: justin gardner
%       date: 06/22/15
%    purpose: Pass subject ID and will return a local repo for that subject. This works
%             by checking whether the local repo exists and if not pulling the desired
%             repo from the central repository. Note that it always tries to update
%             the local repo (pull request). In future, we may want to have a flag
%             to suppress this so that someone can be using a repo without keeping
%             it up to date. 
%
%             Returns:
% 
%             tf: indicates success
%             localRepoROI: string containing directory name where smaller structures like ROIs,
%                           sufaces and flat maps live (of general use)
%             localRepoSession: string containing direcotry name with larger structures like sessions
%                               and other raw data used to create ROIs/surfaces/flat maps (not necessarily of
%                               general use unless you want to check the original data or redraw rois)
%
%             e.g.
%             [tf localRepoROI localRepoSession] = mlrAnatDBGetLocal('s0025');
%
function [tf localRepoROI localRepoSession] = mlrAnatDBGetLocal(subjectID,ROIonly)

% default return arguments
tf = false;
localRepoROI = [];
localRepoSession = [];
if nargin < 2,ROIonly = false;end

% current password
curpwd = pwd;

% get where the anatomy database lives
localRepo = mlrReplaceTilde(mrGetPref('mlrAnatDBLocalRepo'));
centralRepo = mlrReplaceTilde(mrGetPref('mlrAnatDBCentralRepo'));

% check existence of local repo
if isempty(localRepo)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) mlrAnatDBLocalRepo must be set to the location that you want the local repo to be in. You can change this in File/Anat DB/Anat DB Preferences'));
  return
end

% check existence of central repo
if isempty(centralRepo)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) mlrAnatDBCentralRepo must be set to the location (typically an https address or a shared drive) that contains the Anat DB central repository. You can change this in File/Anat DB/Anat DB Preferences'));
  return
end

% make the local repo directory if it does not exist
if ~isdir(localRepo)
  mkdir(localRepo);
end

% check again, if directory exists - in case the mkdir failed so that we can 
% report failure
if ~isdir(localRepo)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not make mlrAnatDB directory %s. Permission problem?',localRepo));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Now get repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetLocal) Getting local repo for %s',subjectID));
localRepoROI = fullfile(localRepo,sprintf('%s',subjectID));
if isdir(localRepoROI)
  % update it
  cd(localRepoROI);
  [status,result] = mysystem(sprintf('hg pull'));
  [status,result] = mysystem(sprintf('hg update'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepoROI);
    return
  else
    disp(sprintf('(mlrAnatDBGetLocal) Successful update of %s',localRepoROI));
  end
else
  disp(sprintf('(mlrAnatDBGetLocal) This may take a few minutes...'));
  centralRepoROI = fullfile(centralRepo,sprintf('%s',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = mysystem(sprintf('hg -v clone %s %s',centralRepoROI,localRepoROI));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clone central Repo %s to local %s',centralRepoROI,localRepoROI));
    return    
  else
    disp(sprintf('(mlrAnatDBGetLocal) Successful clone of %s',centralRepoROI));
  end
end

if ROIonly
  tf = true;
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Now get Session repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetLocal) Getting local session repo for %s',subjectID));
localRepoSession = fullfile(localRepo,sprintf('.%s',subjectID));
if isdir(localRepoSession)
  % update it
  cd(localRepoSession);
  [status,result] = mysystem(sprintf('hg pull'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepoSession);
    return
  end
else
  centralRepoSession = fullfile(centralRepo,sprintf('%sd',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = mysystem(sprintf('hg clone %s %s',centralRepoSession,localRepoSession));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clone central Repo %s to local %s',centralRepoSession,localRepoSession));
    return    
  end
end

% success
tf = true;

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPlugin): %s',command));
[status,result] = system(command,'-echo');
