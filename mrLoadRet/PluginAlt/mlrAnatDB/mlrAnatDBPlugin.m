% mlrAnatDBPlugin
%
%        $Id:$ 
%      usage: mlrAnatDBPlugin(action,<v>)
%         by: justin gardner
%       date: 12/28/2014
%    purpose: Plugin function for mercurial based anatomy database
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
    % add the Add for mlrAnatDB menu
    mlrAdjustGUI(v,'add','menu','Anat DB','/File/ROI','Callback',@mlrAnatDB,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Anat DB Preferences','/File/Anat DB/','Callback',@mlrAnatDBPreferences);
    mlrAdjustGUI(v,'add','menu','Load ROIs from Anat DB','/File/Anat DB/Anat DB Preferences','Callback',@mlrAnatDBLoadROIs,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Load Surfaces from Anat DB','/File/Anat DB/Load ROIs from Anat DB','Callback',@mlrAnatDBLoadSurfaces);
    mlrAdjustGUI(v,'add','menu','Add Session to Anat DB','/File/Anat DB/Load Surfaces from Anat DB','Callback',@mlrAnatDBAddSession,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Add ROIs to Anat DB','/File/Anat DB/Add Session to Anat DB','Callback',@mlrAnatDBAddROIs);
    mlrAdjustGUI(v,'add','menu','Add Surfaces to Anat DB','/File/Anat DB/Add ROIs to Anat DB','Callback',@mlrAnatDBAddSurfaces);
    mlrAdjustGUI(v,'add','menu','Examine ROI in Anat DB','/File/Anat DB/Add Surfaces to Anat DB','Callback',@mlrAnatDBEditROI,'Separator','on');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This plugin support exporting sessions and ROIs to a git managed repository';
 otherwise
   disp(sprintf('(mlrAnatDBPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%
%    mlrAnatDB    %
%%%%%%%%%%%%%%%%%%%
function mlrAnatDB(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get repo locations
centralRepo = mrGetPref('mlrAnatDBCentralRepo');
localRepo = mrGetPref('mlrAnatDBLocalRepo');

% see if the preference is set
if isempty(centralRepo) || isempty(localRepo) 
  % do not enable any thing, because we don't have correct
  % preferences set
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Load Surfaces from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Surfaces to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  return
end

% see if we are in an Anat DB session, get full path to local repo and 
% the current sessions home directory
localRepo = mlrReplaceTilde(localRepo);
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));
% if they are not the same, then offer add session as a menu item,
% but nothing else.
if ~strncmp(localRepo,homeDir,length(localRepo))
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Surfaces from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Surfaces to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
else
  % otherwise don't offer add seesion, but add everything else
  % contingent on whether there are ROIs loaded and so forth
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Surfaces from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Surfaces to Anat DB','Enable','on');

  % see if we have any rois loaded, and gray out Add/AnatDB/ROIs menu accordingly
  if viewGet(v,'nROIs')
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','on');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBPreferences    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBPreferences(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get repo locations
centralRepo = mrGetPref('mlrAnatDBCentralRepo');
localRepo = mrGetPref('mlrAnatDBLocalRepo');

% set defaults
if isempty(centralRepo),centralRepo = '';end
if isempty(localRepo),localRepo = '~/data/mlrAnatDB';end

% setup params info for mrParamsDialog
paramsInfo = {...
    {'mlrAnatDBCentralRepo',centralRepo,'Location of central repo, Typically on a shared server with an https address, but could be on a shared drive in the file structure.'}...
    {'mlrAnatDBLocalRepo',localRepo,'Location of local repo which is typically under a data directory - this will have local copies of ROIs and other data but can be removed the file system as copies will be stored in the central repo'}...
};

% and display the dialog
params = mrParamsDialog(paramsInfo);

% save params, if user did not hit cancel
if ~isempty(params)
  mrSetPref('mlrAnatDBCentralRepo',params.mlrAnatDBCentralRepo,false);
  mrSetPref('mlrAnatDBLocalRepo',params.mlrAnatDBLocalRepo,false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddSession    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddSession(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subjectID (should be sxxxx)
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% remember what directory we started in
curpwd = pwd;

% get the local repos
[tf localRepoROI localRepoSession] = mlrAnatDBGetLocal(subjectID);
if ~tf,return,end

% we should now have a directory with an initialized local repo. 
% so, now we get the branch number
branchNum = mlrAnatDBGetBranchNum(localRepoSession);
if isempty(branchNum),return,end

% Check here to make sure that this session does not already live
% here (as would happen if you are running from that location
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));
if ~strcmp(homeDir,localRepoSession)
  % we are not, so we need to move data into that directory
  if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Will now copy (using hard links) your current session into directory: %s (which is part of the mlrAnatDB). To do so, will need to temporarily close the current session and then reopen in the mlrAnatDB session. Your current work will be saved as usual through the mrLastView mechanism which stores all your current settings. Also, this will not take any more hard disk space, since the files will be copied as hard links. Click OK to continue, or cancel to cancel this operation. If you hit cancel, you will be able to run File/Anat DB/Add Session at a later time, as only a stub directory will have been created in the mlrAnatDB and none of your data will have yet been exported there.',localRepoSession),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
    cd(curpwd);
    return
  end
  % ok, user said we could close, so do it
  mrQuit;
  disppercent(-inf,sprintf('(mlrAnatDBPlugin) Copying %s to %s using hard links.',homeDir,localRepoSession));
  % tag on the name of the session
  localRepoSession = fullfile(localRepoSession,getLastDir(homeDir));
  mkdir(localRepoSession);
  % copy the data from this session over
  [status,result] = system(sprintf('rsync -a --link-dest=%s %s/ %s',homeDir,homeDir,localRepoSession));
  disppercent(inf);
  % check if everything worked ok.
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) rsync seems to have failed to copy data from %s to %s. Switching back to current location',homeDir,localRepoSession);
    cd(homeDir);
    mrLoadRet;
    cd(curpwd);
    return
  end
  % everything went ok, switch directories and start up over there
  cd(localRepoSession);
  curpwd = localRepoSession;
  mrLoadRet;
end

% add link to localizers in ROI repo
cd(fullfile(localRepoROI,'Localizers'));
[status,result] = system(sprintf('ln -s ../../%s %s',getLastDir(localRepoSession,2),getLastDir(localRepoSession)));
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not make link for %s in ROI database',getLastDir(localRepoSession)));
else
  % set branch number
  if ~mlrAnatDBSetBranchNum(localRepoROI,branchNum+1), return,end
  % add/commit/push the repo
  if ~mlrAnatDBAddCommitPush(localRepoROI,'Localizers',sprintf('Adding link to %s',getLastDir(localRepoSession)));
    return
  end
end


% set branch number
if ~mlrAnatDBSetBranchNum(fileparts(localRepoSession),branchNum+1), return,end

% add this directory to the local repo
% debug - FIX, FIX, FIX get comments from user to add as commit
if ~mlrAnatDBAddCommitPush(fileparts(localRepoSession),getLastDir(localRepoSession),'Saving initial session snapshot',true,true)
  cd(curpwd);
  return
end

cd(curpwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject number from directory
subjectID = getLastDir(fileparts(viewGet(v,'homeDir')));
subjectID = subjectID(2:end);

% get the local repos
[tf localRepoROI localRepoSession] = mlrAnatDBGetLocal(subjectID,true);
if ~tf,return,end

% get list of ROIs to save
roiList = selectInList(v,'rois','Select ROI(s) to save');

% where to save them to
localROIDir = fullfile(localRepoROI,'ROIs');

% get current branch number
branchNum = mlrAnatDBGetBranchNum(localRepoSession);
if isempty(branchNum),return,end

% increment branch num
if ~mlrAnatDBSetBranchNum(localRepoROI,branchNum+1), return,end
if ~mlrAnatDBSetBranchNum(localRepoSession,branchNum+1), return,end

% and save them
for iRoi = roiList
  % FIx, FIx, FIx, set new parameters for ROI here (like username
  % etc)
  saveROI(v,iRoi,false,localROIDir);
end

% Fix, FIx, FIx, 
% get user comments
comments = 'Saving ROIS';

% add/commit/push roi repo
if ~mlrAnatDBAddCommitPush(localRepoROI,'ROIs',comments)
  return
end

% update a bogus file with version number to make sure a version gets saved
bogusFile = fullfile(localRepoSession,'.mlrAnatDB');
[status,results] = system('echo s%04i >> %s',branchNum+1,bogusFile);

% add/commit/push
if ~mlrAnatDBAddCommitPush(localRepoSession,bogusFile,comments,true,true)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddSurfaces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddSurfaces(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBEditROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBEditROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadSurfaces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadSurfaces(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');



%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBCheckHg    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBSubjectID    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subjectID = mlrAnatDBSubjectID(v)

% get the subject
subjectID = viewGet(v,'subject');
if (length(subjectID) > 1) && isequal(subjectID(1),'s')
  % format subject ID
  subjectID = sprintf('s%04i',str2num(subjectID(2:end)));
end

% ask if this is the correct subjectID
paramsInfo{1} = {'subjectID',subjectID,'The subject ID that this session will be filed under in the anatDB. Usually this is of the form sXXXX. If you do not know it, you may be able to look it up using mglSetSID if you are usuing the mgl ID database'};
params = mrParamsDialog(paramsInfo,'Set subjectID');
if isempty(params)
  subjectID = '';
  return
end
subjectID = params.subjectID;
% check input format
if length(subjectID) < 5
  subjectID = sprintf('s%04i',str2num(subjectID(2:end)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBGetLocal    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Now get ROI repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetLocal) Getting local ROI repo for %s',subjectID));
localRepoROI = fullfile(localRepo,sprintf('%s',subjectID));
if isdir(localRepoROI)
  % update it
  cd(localRepoROI);
  [status,result] = system(sprintf('hg pull'));
  [status,result] = system(sprintf('hg update'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepoROI);
    return
  end
else
  centralRepoROI = fullfile(centralRepo,sprintf('%s',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = system(sprintf('hg clone %s %s',centralRepoROI,localRepoROI));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clonal central Repo %s to local %s',centralRepoROI,localRepoROI));
    return    
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Now get Session repo
%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('(mlrAnatDBGetLocal) Getting local session repo for %s',subjectID));
localRepoSession = fullfile(localRepo,sprintf('.%s',subjectID));
if isdir(localRepoSession)
  % update it
  cd(localRepoSession);
  [status,result] = system(sprintf('hg pull'));
  cd(curpwd);
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) Unable to update local Repo %s',localRepoSession);
    return
  end
else
  centralRepoSession = fullfile(centralRepo,sprintf('%sd',subjectID));
  % try to retrieve from remote repo by cloning
  [status,result] = system(sprintf('hg clone %s %s',centralRepoSession,localRepoSession));
  % if successful, then we have it
  if status~=0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Unable to clonal central Repo %s to local %s',centralRepoSession,localRepoSession));
    return    
  end
end

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddCommitPush    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrAnatDBAddCommitPush(localRepo,filename,comments,bigfiles,verbose)

if nargin < 4,bigfiles = false;end
if nargin < 5,verbose = false;end

% change paths
tf = false;
curpwd = pwd;
cd(localRepo);

disp(sprintf('(mlrAnatDBAddCommitPush) Adding files to repo %s',localRepo));

% add file to repo
if bigfiles
  [status,result] = system(sprintf('hg add %s',filename));
else
  [status,result] = system(sprintf('hg add %s --large',filename));
end
  
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not add files to local repo. Have you setup your config file for hg?'));
  cd(curpwd);
  return
end

% commit to repo
if bigfiles
  disppercent(-inf,sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s. This may take a minute or two...',localRepo));
else
  disp(sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s',localRepo));
end
  
[status,result] = system(sprintf('hg commit -m ''%s''',comments));
if bigfiles,disppercent(inf);,end
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not commit files to local repo. Have you setup your config file for hg?'));
  cd(curpwd);
  return
end

% and push
if ~verbose || isequal(questdlg(sprintf('Do you want to push to the central repo: %s? This can take several minutes depending on your connection. You will be able to work while this occurs (it will push in the background), but you should not shut off your matlab/computer. If you choose no now, you will need to push manually later.',mrGetPref('mlrAnatDBCentralRepo')),'Do push?','Yes','No','Yes'),'Yes')
  disp(sprintf('(mlrAnatDBAddCommitPush) Pushing repo %s in the background',localRepo));
  system(sprintf('hg push --new-branch &'));
end

tf = true;
cd(curpwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBGetBranchNum    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function branchNum = mlrAnatDBGetBranchNum(localRepo)

branchNum = [];

% set path
curpwd = pwd;
cd(localRepo);

% get the branch name
[status,result] = system(sprintf('hg branch'));
branchNameLoc = regexp(result,'v\d');
if ~isempty(branchNameLoc)
  branchNum = str2num(result(branchNameLoc+1:end));
else
  mrWarnDlg(sprintf('(mlrAnatDbPlugin) Could not figure out version number. This should be the current branch of the git repository and should be in the format vXXXX where XXX is a number (e.g. v0001). Not able to commit changes. Aborting. you can fix by going to repo %s and assiging a valid version number as the branch name'));
  cd(curpwd);
  return
end

cd(curpwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBSetBranchNum    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrAnatDBSetBranchNum(localRepo,branchNum)

tf = false;

% set path
curpwd = pwd;
cd(localRepo);

% update branch number
branchName = sprintf('v%04i',branchNum);
[status,result] = system(sprintf('hg branch %s',branchName));
if status == 0
  tf = true;
end

cd(curpwd);