% mlrAnatDBGetRepo.m
%
%        $Id:$ 
%      usage: [localRepo] = mlrAnatDBGetRepo(subjectID)
%         or: [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID)
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
%             localRepo: string containing directory name where smaller structures like ROIs,
%                           sufaces and flat maps live (of general use). Returns empty on failure
%             localRepoLargeFiles: string containing direcotry name with larger structures like sessions
%                               and other raw data used to create ROIs/surfaces/flat maps (not necessarily of
%                               general use unless you want to check the original data or redraw rois)
%                               if you do not need these, then do not accept this argument (this is
%                               advised since getting this argument forces the large files repo to brought
%                               down from the server causing it to take up lots of disk space)
%
%             e.g.
%             localRepo = mlrAnatDBGetRepo('s0025');
%
%             or if you need to store large files like sessions and freesurfer directories
%
%             [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo('s0025');
%
function [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID)

% default return arguments
localRepo = [];
localRepoLargeFiles = [];

% current password
curpwd = pwd;

% get where the anatomy database lives
localRepoTop = mlrReplaceTilde(mrGetPref('mlrAnatDBLocalRepo'));
centralRepoTop = mlrReplaceTilde(mrGetPref('mlrAnatDBCentralRepo'));

% check existence of local repo
if isempty(localRepoTop)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) mlrAnatDBLocalRepo must be set to the location that you want the local repo to be in. You can change this in File/Anat DB/Anat DB Preferences'));
  return
end

% check existence of central repo
if isempty(centralRepoTop)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) mlrAnatDBCentralRepo must be set to the location (typically an https address or a shared drive) that contains the Anat DB central repository. You can change this in File/Anat DB/Anat DB Preferences'));
  return
end

% make the local repo directory if it does not exist
if ~isdir(localRepoTop)
  mkdir(localRepoTop);
end

% check again, if directory exists - in case the mkdir failed so that we can 
% report failure
if ~isdir(localRepoTop)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not make mlrAnatDB directory %s. Permission problem?',localRepoTop));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Now get repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetRepo) Getting local repo for %s',subjectID));
localRepo = fullfile(localRepoTop,sprintf('%s',subjectID));
if isdir(localRepo)
  % update it
  cd(localRepo);
  [status,result] = mysystem(sprintf('hg pull'));
  [status,result] = mysystem(sprintf('hg update'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepo);
    localRepo = [];
    return
  else
    disp(sprintf('(mlrAnatDBGetRepo) Successful update of %s',localRepo));
  end
else
  disp(sprintf('(mlrAnatDBGetRepo) This may take a few minutes...'));
  centralRepo = fullfile(centralRepoTop,sprintf('%s',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = mysystem(sprintf('hg -v clone %s %s',centralRepo,localRepo));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clone central Repo %s to local %s',centralRepo,localRepo));
    localRepo = [];
    return    
  else
    disp(sprintf('(mlrAnatDBGetRepo) Successful clone of %s',centralRepo));
  end
end

% if not getting largefiles (session repo then stop here
if nargout < 2
  tf = true;
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Now get Session repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetRepo) Getting local session repo for %s',subjectID));
localRepoLargeFiles = fullfile(localRepoTop,sprintf('.%s',subjectID));
if isdir(localRepoLargeFiles)
  % update it
  cd(localRepoLargeFiles);
  [status,result] = mysystem(sprintf('hg pull'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepoLargeFiles);
    localRepoLargeFiles = [];
    return
  end
else
  centralRepoLargeFiles = fullfile(centralRepoTop,sprintf('%sd',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = mysystem(sprintf('hg clone %s %s',centralRepoLargeFiles,localRepoLargeFiles));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clone central Repo %s to local %s',centralRepoLargeFiles,localRepoLargeFiles));
    localRepoLargeFiles = [];
    return    
  end
end

% now make links in local repo
curpwd = pwd;
cd(localRepo);
linkList = {'3D','localizers'};
for iLink = 1:length(linkList)
  linkFrom = fullfile('..',getLastDir(localRepoLargeFiles),linkList{iLink});
  system(sprintf('ln -sfh %s %s',linkFrom,linkList{iLink}));
end
cd(curpwd);

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPlugin): %s',command));
[status,result] = system(command,'-echo');
