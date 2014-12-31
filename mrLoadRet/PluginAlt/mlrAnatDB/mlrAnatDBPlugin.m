% mlrAnatDBPlugin
%
%        $Id:$ 
%      usage: mlrAnatDBPlugin(action,<v>)
%         by: justin gardner
%       date: 12/28/2014
%    purpose: Plugin function for git based anatomy database
%
function retval = mlrAnatDBPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrAnatDBPlugin) Need a valid view to install plugin'));
  else
    % add the Export for mlrAnatDB menu
    mlrAdjustGUI(v,'add','menu','Anat DB','/File/Export/','Callback',@mlrAnatDBExport);
    mlrAdjustGUI(v,'add','menu','Session','/File/Export/Anat DB/','Callback',@mlrAnatDBSession);
    mlrAdjustGUI(v,'add','menu','ROIs','/File/Export/Anat DB/Session','Callback',@mlrAnatDBROIs);

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This plugin support exporting sessions and ROIs to a git managed repository';
 otherwise
   disp(sprintf('(mlrAnatDBPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBExport    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBExport(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% see if we have any rois loaded, and gray out Export/AnatDB/ROIs menu accordingly
if viewGet(v,'nROIs')
  mlrAdjustGUI(v,'set','ROIs','Enable','on');
else
  mlrAdjustGUI(v,'set','ROIs','Enable','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBSession    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBSession(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get version control system (should be either mercurial -preferred or git)
vcs = mrGetPref('mlrAnatDBVCS');

% set up the commands used to run the version control system
vcsInit = 'init';
vcsClone = 'clone';
vcsAdd = 'add';
vcsVerbose ='-v';
vcsCommit = 'commit';
vcsPush = 'push';
vcsPull = 'pull';
if strcmp(lower(vcs),'mercurial')
  vcs = 'hg';
  vcsLargeFiles = '--large';
  % DEBUG: Need to enable largefiles extension
  % and make sure user has put their username correctly
%  system('echo ''[extensions]'' > .hg/hgrc');
%  system('echo ''largefiles ='' >> .hg/hgrc');
elseif strcmp(lower(vcs),'git')
  vcs = 'git';
  % git does not have any large file facility
  vcsLargeFiles = '';
  vcsAdd = 'add';
else
  mrWarnDlg('(mlrAnatDBPlugin) Your mlrAnatDBVCS variable (Edit/Preferences) SHould be set to one of mercurial (preferred) or git. Aborting');
  return
end
  
% check for VCS
[status,result] = system(sprintf('which %s',vcs));
if status ~= 0
  if strcmp(vcs,'git')
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) You do not have git installed. You will need to install git - typically on Mac OS by installing XCode with the Command Line Tools'));
    return
  else
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) You do not have mercurial installed. You will need to install mercurial. Typicaly by going to the website: http://mercurial.selenic.com and following download instructions.'));
    return
  end
end

% get where the anatomy database lives
mlrAnatDir = mrGetPref('volumeDirectory');
if isempty(mlrAnatDir)
  % set the default location for the mlrAnatDB
  mlrAnatDir = '~/data/mlrAnatDB';
  mlrAnatDir = mlrReplaceTilde(mlrAnatDir);
  % and save it in prefs
  mrSetPref('volumeDirectory',mlrAnatDir);
else
  % swap tilde for fullly qualified path
  mlrAnatDir = mlrReplaceTilde(mlrAnatDir);
end

% make the directory if it does not exist
if ~isdir(mlrAnatDir)
  mkdir(mlrAnatDir);
end

% check again, if directory exists - in case the mkdir failed so that we can 
% report failure
if ~isdir(mlrAnatDir)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not make mlrAnatDB directory %s. Permission problem?',mlrAnatDir));
  return
end

% get the subject
subjectID = viewGet(v,'subject');

% ask if this is the correct subjectID
paramsInfo{1} = {'subjectID',subjectID,'The subject ID that this session will be filed under in the anatDB. Usually this is of the form sXXX. If you do not know it, you may be able to look it up using mglSetSID if you are usuing the mgl ID database'};
params = mrParamsDialog(paramsInfo,'Set subjectID');
if isempty(params),return,end

% get central repo
mlrAnatDBCentralRepo = mrGetPref('mlrAnatDBCentralRepo');
if isempty(mlrAnatDBCentralRepo)
  h = oneTimeWarning('mlrAntDBCentralRepoEmpty','You have not set mlrAnatDBCentralRepo in your Edit/Preferences. This should be set to the remote central repository (typically on a shared server) that you use to centrally store your anat database. These functions will still run and work fine, but your repo will not be backed up to the central repository and only exist locally on your computer',false);
  if ~isempty(h),uiwait(h),end
end

% we should have the anat repo and subject ID, check to see if we have an entry for the subject
mlrAnatDirSession = fullfile(mlrAnatDir,sprintf('.%s',subjectID));
mlrAnatDirSessionExists = false;
if isdir(mlrAnatDirSession)
  % ok, it already exists a directory
  mlrAnatDirSessionExists = true;
else
  % see if we have a remote repo name
  if ~isempty(mlrAnatDBCentralRepo)
    % try to retrieve from git remote repo by cloning
    [status,result] = system(sprintf('%s %s %s %s',vcs,vcsClone,mlrAnatDBCentralRepo,mlrAnatDirSession));
    % if successful, then we have it
    if status==0
      mlrAnatDirSessionExists = true;
    end
  end
end

% remember what directory we started in
curpwd = pwd;

% if we don't have the session directory for this subject, then create it
if ~mlrAnatDirSessionExists
  % make the directory
  mkdir(mlrAnatDirSession)
  if ~isdir(mlrAnatDirSession)
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not make mlrAnatDB directory for subject: %s. Permission problem?',mlrAnatDirSession));
    return
  end
  % init local repo
  cd(mlrAnatDirSession);
  [status,result] = system(sprintf('%s %s',vcs,vcsInit));
  if status~=0
    disp(sprintf('(mlrAnatDBPlugin) %s init has failed on directory %s',vcs,mlrAnatDirSession));
    cd(curpwd);
    return
  end
  % add a single bogus file just so that we can start branching correctly
  % this is for git land
  if strcmp(vcs,'git');
    system('touch .mlrAnatDBInit');
    system('git add .mlrAnatDBInit');
    system('git commit -m ''Init repo''');
    % change the first branch to be called v0000
    system('git branch -m master v0000');
  else
    % simpler in mercurial
    system('hg branch v0000');
  end
end

% we should now have a directory with an initialized local repo. 
% so, now we update the branch number
if strcmp(vcs,'git')
  [status,result] = system(sprintf('git branch | grep ''*'''));
else
  % again, simpler in mercurial (just shows you current branch)
  [status,result] = system(sprintf('hg branch'));
end
branchNameLoc = regexp(result,'v\d');
if ~isempty(branchNameLoc)
  branchNum = str2num(result(branchNameLoc+1:end));
else
  mrWarnDlg(sprintf('(mlrAnatDbPlugin) Could not figure out version number. This should be the current branch of the git repository and should be in the format vXXXX where XXX is a number (e.g. v0001). Not able to commit changes. Aborting. you can fix by going to repo %s and assiging a valid version number as the branch name'));
  cd(curpwd);
  return
end
% update branch number
branchName = sprintf('v%04i',branchNum);
if strcmp(vcs,'git')
  [status,result] = system(sprintf('git checkout -b %s',branchName));
else
  [status,result] = system(sprintf('hg branch %s',branchName));
end
disp(sprintf('!!!! DEBUG: What happens here if you have uncommitted changes in the current branch? !!!!'));

% Check here to make sure that this session does not already live
% here (as would happen if you are running from that location
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));
if ~strcmp(homeDir,mlrAnatDirSession)
  % we are not, so we need to move data into that directory
  if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Will now copy (using hard links) your current session into directory: %s (which is part of the mlrAnatDB). To do so, will need to temporarily close the current session and then reopen in the mlrAnatDB session. Your current work will be saved as usual through the mrLastView mechanism which stores all your current settings. Also, this will not take any more hard disk space, since the files will be copied as hard links. Click OK to continue, or cancel to cancel this operation. If you hit cancel, you will be able to run File/Export/AnatDB/Session at a later time, as only a stub directory will have been created in the mlrAnatDB and none of your data will have yet been exported there.',mlrAnatDirSession),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
    cd(curpwd);
    return
  end
  % ok, user said we could close, so do it
  mrQuit;
  disppercent(-inf,sprintf('(mlrAnatDBPlugin) Copying %s to %s using hard links.',homeDir,mlrAnatDirSession));
  % tag on the name of the session
  mlrAnatDirSession = fullfile(mlrAnatDirSession,getLastDir(homeDir));
  mkdir(mlrAnatDirSession);
  % copy the data from this session over
  [status,result] = system(sprintf('rsync -a --link-dest=%s %s/ %s',homeDir,homeDir,mlrAnatDirSession));
  disppercent(inf);
  % check if everything worked ok.
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) rsync seems to have failed to copy data from %s to %s. Switching back to current location',homeDir,mlrAnatDirSession);
    cd(homeDir);
    mrLoadRet;
    cd(curpwd);
    return
  end
  % everything went ok, switch directories and start up over there
  cd(mlrAnatDirSession);
  curpwd = mlrAnatDirSession;
  mrLoadRet;
end

% add this directory to the local vcs repo
cd('..');
disppercent(-inf,sprintf('(mlrAnatDBPlugin) Adding files to local %s repo',vcs));
[status,result] = system(sprintf('%s %s %s %s %s',vcs,vcsAdd,getLastDir(mlrAnatDirSession),vcsLargeFiles,vcsVerbose));
disppercent(inf);
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not add files to local repo. Have you setup your config file for %s?',vcs));
  cd(curpwd);
  return
end

% commit files
disppercent(-inf,sprintf('(mlrAnatDBPlugin) Committing files to local %s repo. This may also take some time.',vcs));
[status,result] = system(sprintf('%s %s -m ''Saving snapshot labeled %s from MLR''',vcs,vcsCommit,branchName));
disp(sprintf('(mlrAnatDBPlugin) Done committing files to repo'));
disppercent(inf);
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not commit files to local repo. Have you setup your config file for %s?',vcs));
  cd(curpwd);
  return
end

% DEBUG: now need to upload to database
keyboard

cd(curpwd);
    

  




    

  


%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBROIs    %
%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

