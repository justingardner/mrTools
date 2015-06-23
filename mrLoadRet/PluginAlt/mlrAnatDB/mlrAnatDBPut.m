% mlrAnatDBPut.m
%
%        $Id:$ 
%      usage: mlrAnatDBPut(subjectID,filePath,<fileType=[]>,<toPath=[]>,<bigfiles=0>,<verbose=0>,<comments=[]>)
%         by: justin gardner
%       date: 06/22/15
%    purpose: 
%
function tf = mlrAnatDBPut(subjectID,filePath,varargin)

tf = false;
% check arguments
if nargin < 2
  help mlrAnatDBPut;
  return;
end

% and get arguments
getArgs(varargin,{'fileType=[]','fileDir=[]','bigfiles=0','verbose=0','comments=[]'});

% check file
if ~isfile(filePath) && ~isdir(filePath)
  disp(sprintf('(mlrAnatDBPut) Could not find: %s',filePath));
  return
end

% first get local repo
[success localRepo localBigRepo] = mlrAnatDBGetLocal(subjectID,~bigfiles);
if ~success,return,end

% figure out what path to put it under
if ~isempty(fileDir)
  switch lower(fileType)
   case {'roi','rois'}
    fileDir = 'mlrROIs';
   case {'3d','freesurfer','canonical'}
    fileDir = '3D';
   case {'baseanat','mlrbaseAnat','mlrbaseanatomy','mlrbaseanatomies'}
    fileDir = 'mlrBaseAnatomies';
   case {'localizer','localizers'}
    fileDir = 'localizers';
  end    
end
if isempty(fileDir)
  disp(sprintf('(mlrAnatDBPut) Must specify a vaild fileType or fileDir'));
  return
end

% get location to copy to 
if bigfiles
  destDir = fullfile(localBigRepo,fileDir);
else
  destDir = fullfile(localRepo,fileDir);
end

% make directory if necessary
if ~isdir(destDir), mkdir(destDir); end

% copy using force
disp(sprintf('(mlrAnatDBPut) Copying %s to %s',filePath,destDir));
success = copyfile(filePath,destDir,'-f');
if ~success,return,end

curpwd = pwd;
cd(localRepo);

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

