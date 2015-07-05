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
    mlrAdjustGUI(v,'add','menu','Load Base Anatomies from Anat DB','/File/Anat DB/Load ROIs from Anat DB','Callback',@mlrAnatDBLoadBaseAnatomies);
    mlrAdjustGUI(v,'add','menu','Add Session to Anat DB','/File/Anat DB/Load Base Anatomies from Anat DB','Callback',@mlrAnatDBAddSession,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Add ROIs to Anat DB','/File/Anat DB/Add Session to Anat DB','Callback',@mlrAnatDBAddROIs);
    mlrAdjustGUI(v,'add','menu','Add Base Anatomies to Anat DB','/File/Anat DB/Add ROIs to Anat DB','Callback',@mlrAnatDBAddBaseAnatomies);
    mlrAdjustGUI(v,'add','menu','Examine ROI in Anat DB','/File/Anat DB/Add Base Anatomies to Anat DB','Callback',@mlrAnatDBExamineROI,'Separator','on');

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
localRepoTop = mrGetPref('mlrAnatDBLocalRepo');

% see if the preference is set
if isempty(centralRepo) || isempty(localRepoTop) 
  % do not enable any thing, because we don't have correct
  % preferences set
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  return
end

% see if we are in an Anat DB session
if ~mlrAnatDBInLocalRepo(v)
  % if not, then only allow add session and examine ROI
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','off');
  if viewGet(v,'nROIs')
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  end
else
  % otherwise don't offer add seesion, but add everything else
  % contingent on whether there are ROIs loaded and so forth
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');

  % see if any bases are loaded
  if viewGet(v,'numBase')
    mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','off');
  end    
  
  % see if we have any rois loaded, and gray out Add/AnatDB/ROIs menu accordingly
  if viewGet(v,'nROIs')
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','on');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  end
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
[localRepoSubject localRepoSubjectLargeFiles] = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject) || isempty(localRepoSubjectLargeFiles),return,end

% we should now have a directory with an initialized local repo. 
% so, now we get the branch number
branchNum = mlrAnatDBGetBranchNum(localRepoSubjectLargeFiles);
if isempty(branchNum),return,end

% Check here to make sure that this session does not already live
% here (as would happen if you are running from that location
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));
if ~strcmp(homeDir,localRepoSubjectLargeFiles)
  % we are not, so we need to move data into that directory
  if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Will now copy (using hard links) your current session into directory: %s (which is part of the mlrAnatDB). To do so, will need to temporarily close the current session and then reopen in the mlrAnatDB session. Your current work will be saved as usual through the mrLastView mechanism which stores all your current settings. Also, this will not take any more hard disk space, since the files will be copied as hard links. Click OK to continue, or cancel to cancel this operation. If you hit cancel, you will be able to run File/Anat DB/Add Session at a later time, as only a stub directory will have been created in the mlrAnatDB and none of your data will have yet been exported there.',localRepoSubjectLargeFiles),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
    cd(curpwd);
    return
  end
  % ok, user said we could close, so do it
  mrQuit;
  localRepoSubjectLargeFiles = fullfile(localRepoSubjectLargeFiles,'localizers',getLastDir(homeDir));
  disppercent(-inf,sprintf('(mlrAnatDBPlugin) Copying %s to %s using hard links.',homeDir,localRepoSubjectLargeFiles));
  % make the directory
  mkdir(localRepoSubjectLargeFiles);
  % copy the data from this session over
  [status,result] = system(sprintf('rsync -a --link-dest=%s %s/ %s',homeDir,homeDir,localRepoSubjectLargeFiles));
  disppercent(inf);
  % check if everything worked ok.
  if status ~= 0
    mrWarnDlg('(mlrAnatDBPlugin) rsync seems to have failed to copy data from %s to %s. Switching back to current location',homeDir,localRepoSubjectLargeFiles);
    cd(homeDir);
    mrLoadRet;
    cd(curpwd);
    return
  end
  % everything went ok, switch directories and start up over there
  cd(localRepoSubjectLargeFiles);
  curpwd = localRepoSubjectLargeFiles;
  mrLoadRet;
end

% set branch number
if ~mlrAnatDBSetBranchNum(localRepoSubject,branchNum+1), return,end

% set branch number
if ~mlrAnatDBSetBranchNum(fileparts(fileparts(localRepoSubjectLargeFiles)),branchNum+1), return,end

% add this directory to the local repo
% debug - FIX, FIX, FIX get comments from user to add as commit
if ~mlrAnatDBAddCommitPush(fileparts(localRepoSubjectLargeFiles),getLastDir(localRepoSubjectLargeFiles),'Saving initial session snapshot',true,true)
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
subjectID = mlrAnatDBSubjectID(v);

% get the local repos
[localRepoSubject localRepoSubjectLargeFiles] = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject) || isempty(localRepoSubjectLargeFiles),return,end

% get list of ROIs to save
roiList = selectInList(v,'rois','Select ROI(s) to save');

% where to save them to
localROIDir = fullfile(localRepoSubject,'mlrROIs');

% sessin name
sessionName = getLastDir(viewGet(v,'homeDir'));

% and save them
for iRoi = roiList
  % FIx, FIx, FIx, set new parameters for ROI here (like username
  % etc)
  v = viewSet(v,'roiCreatedFromSession',sessionName,iRoi);
  saveROI(v,iRoi,false,localROIDir);
end

% Fix, FIx, FIx, 
% get user comments
comments = 'Saving ROIS';

% get current branch number
branchNum = mlrAnatDBGetBranchNum(localRepoSubject);
if isempty(branchNum),return,end

% increment branch num
if ~mlrAnatDBSetBranchNum(localRepoSubject,branchNum+1), return,end
if ~mlrAnatDBSetBranchNum(localRepoSubjectLargeFiles,branchNum+1), return,end

% add/commit/push roi repo
if ~mlrAnatDBAddCommitPush(localRepoSubject,'mlrROIs',comments)
  return
end

% update a bogus file with version number to make sure a version gets saved
bogusFile = fullfile(localRepoSubjectLargeFiles,'.mlrAnatDB');
[status,results] = system(sprintf('echo s%04i >> %s',branchNum+1,bogusFile));

% add/commit/push
if ~mlrAnatDBAddCommitPush(localRepoSubjectLargeFiles,bogusFile,comments,true,true)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddBaseAnatomies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddBaseAnatomies(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject number from directory
subjectID = mlrAnatDBSubjectID(v);

% get the local repos
[localRepoSubject localRepoSubjectLargeFiles] = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject) || isempty(localRepoSubjectLargeFiles),return,end

% get list of ROIs to save
baseList = selectInList(v,'bases','Select MLR Base Anatomies to save');

% where to save them to
localBaseAnatomyDir = fullfile(localRepoSubject,'mlrBaseAnatomies');

% and save them
for iBase = baseList
  saveAnat(v,iBase,false,false,localBaseAnatomyDir);
end

% Fix, FIx, FIx, 
% get user comments
comments = 'Saving mlrBaseAnatomise';

% get current branch number
branchNum = mlrAnatDBGetBranchNum(localRepoSubject);
if isempty(branchNum),return,end

% increment branch num
if ~mlrAnatDBSetBranchNum(localRepoSubject,branchNum+1), return,end
if ~mlrAnatDBSetBranchNum(localRepoSubjectLargeFiles,branchNum+1), return,end

% add/commit/push roi repo
if ~mlrAnatDBAddCommitPush(localRepoSubject,'mlrROIs',comments)
  return
end

% update a bogus file with version number to make sure a version gets saved
bogusFile = fullfile(localRepoSubjectLargeFiles,'.mlrAnatDB');
[status,results] = system(sprintf('echo s%04i >> %s',branchNum+1,bogusFile));

% add/commit/push
if ~mlrAnatDBAddCommitPush(localRepoSubjectLargeFiles,bogusFile,comments,true,true)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBExamineROI    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBExamineROI(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% select an ROI to examine
roiNames = viewGet(v,'roiNames');
curROI = viewGet(v,'curROI');
if (curROI >1) && (length(roiNames)>=curROI)
  roiNames = putOnTopOfList(roiNames{curROI},roiNames);
end
paramsInfo = {{'roiToExamine',roiNames,'Choose which ROI to examine in its original localizer session'}};
params = mrParamsDialog(paramsInfo,'Choose ROI');
if isempty(params),return,end

% get the original session that the ROI was defined in and make
% sure we have it in the repo
roi = viewGet(v,'roi',params.roiToExamine);

% check for createdFromSession
if isempty(roi.createdFromSession)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) The ROI %s does not have the field createdFromSession set. Not sure which session it was created from',roi.name));
  return
end

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
[localRepoSubject localRepoSubjectLargeFiles] = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject) || isempty(localRepoSubjectLargeFiles),return,end

% check to see if session exists in repo
createdFromSession = fullfile(localRepoSubjectLargeFiles,roi.createdFromSession);
if ~isdir(createdFromSession)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) Could not find session %s in mlrAnatDB',createdFromSession));
  return
end
  
% now confirm that this is what the user really wants to do
if strcmp(questdlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) Will now close this current session (saving work as always) and will load up the localizer session where ROI: %s was defined',params.roiToExamine),'Switch to Localizer session','Ok','Cancel','Cancel'),'Cancel')
  return
end

% ok, user said we could close, so do it
mrQuit;
cd(createdFromSession);
mrLoadRet
v = getMLRView;

% check to see if ROI is loaded
if ~any(strcmp(viewGet(v,'roiNames'),params.roiToExamine))
  % then load it
  v = loadROI(v,params.roiToExamine,0,fullfile(localRepoSubject,'mlrROIs'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
localRepoSubject = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject),return,end

% load the rois
v = loadROI(v,[],[],fullfile(localRepoSubject,'mlrROIs'));

% and refresh
refreshMLRDisplay(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadBaseAnatomies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadBaseAnatomies(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
localRepoSubject = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject),return,end

% load the rois
v = loadAnat(v,[],fullfile(localRepoSubject,'mlrBaseAnatomies'));

% and refresh
refreshMLRDisplay(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddCommitPush    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrAnatDBAddCommitPush(localRepoTop,filename,comments,bigfiles,verbose)

if nargin < 4,bigfiles = false;end
if nargin < 5,verbose = false;end

% change paths
tf = false;
curpwd = pwd;
cd(localRepoTop);

disp(sprintf('(mlrAnatDBAddCommitPush) Adding files to repo %s',localRepoTop));

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
  disppercent(-inf,sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s. This may take a minute or two...',localRepoTop));
else
  disp(sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s',localRepoTop));
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
  disp(sprintf('(mlrAnatDBAddCommitPush) Pushing repo %s in the background',localRepoTop));
  system(sprintf('hg push --new-branch &'));
end

tf = true;
cd(curpwd);

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPlugin): %s',command));
[status,result] = system(command,'-echo');

