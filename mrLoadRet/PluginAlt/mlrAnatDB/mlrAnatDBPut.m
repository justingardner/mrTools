% mlrAnatDBPut.m
%
%        $Id:$ 
%      usage: mlrAnatDBPut(subjectID,filePath,fileType,<toPath=[]>,<largefiles=[]>,<verbose=0>,<comments=[]>)
%         by: justin gardner
%       date: 06/22/15
%    purpose: Puts files into repo. The filePath should point to a directory outside of mlrAnatDB and
%             that file or directory will be copied into an appropriate place into mlrAnatDBm added to the
%             repo, committed and pushed (if the user accepts)
%
%             Typically you specify what kind of file you have (roi, 3D, baseAnat, localizer) and
%             this will put it into the correct directory
% 
%             e.g. To put an roi
%             mlrAnatDBPut(25,'~/data/mtroi.mat','roi');
%
%             e.g. To put a directory of rois
%             mlrAnatDBPut(25,'~/data/rois','roi');
%
%             e.g. To put a freesurfer directory
%             mlrAnatDBPut(25,'~/data/freesurfer','freesurfer');
%
%             Valid types can be found in the code below: roi, 3D, baseAnat, freesurfer, canonical, localizer
%
function tf = mlrAnatDBPut(subjectID,filePath,fileType,varargin)

tf = false;
% check arguments
if nargin < 2
  help mlrAnatDBPut;
  return;
end

% format the subject id
subjectID = formatSubjectID(subjectID);

% and get arguments
getArgs(varargin,{'fileDir=[]','largefiles=[]','verbose=0','comments=[]'});

% check file
if ~isfile(filePath) && ~isdir(filePath)
  disp(sprintf('(mlrAnatDBPut) Could not find: %s',filePath));
  return
end

% figure out what path to put it under and whether it is a large file or not
if isempty(fileDir)
  switch lower(fileType)
   case {'roi','rois'}
    fileDir = 'mlrROIs';
    if isempty(largefiles),largefiles = false;end
    if isdir(filePath),filePath = getFilenames(filePath,'*.mat');end
   case 'freesurfer'
    fileDir = '3D';
    if isempty(largefiles),largefiles = true;end
   case {'3d','canonical'}
    fileDir = '3D';
    if isempty(largefiles),largefiles = false;end
   case {'baseanat','mlrbaseAnat','mlrbaseanatomy','mlrbaseanatomies'}
    fileDir = 'mlrBaseAnatomies';
    if isempty(largefiles),largefiles = false;end
    if isdir(filePath),filePath = getFilenames(filePath,'*.mat');end
   case {'surfaces','surface'}
    fileDir = 'surfaces';
    if isempty(largefiles),largefiles = false;end
    if isdir(filePath),filePath = getFilenames(filePath,{'*.off','*.vff'});end
   case {'localizer','localizers'}
    fileDir = 'localizers';
    if isempty(largefiles),largefiles = true;end
  end    
end
if isempty(fileDir)
  disp(sprintf('(mlrAnatDBPut) Must specify a vaild fileType or fileDir'));
  return
end

% get local repo
if ~largefiles
  [localRepo] = mlrAnatDBGetRepo(subjectID);
else
  [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID);
end
if isempty(localRepo),return,end

% get location to copy to 
if largefiles
  destDir = fullfile(localRepoLargeFiles,fileDir);
else
  destDir = fullfile(localRepo,fileDir);
end

% make directory if necessary
if ~isdir(destDir), mkdir(destDir); end

% make sure filePath is a cell array (i.e. one cell for each file)
filePath = cellArray(filePath);

% remember current directory
curpwd = pwd;

% copy each file using force
for iFile = 1:length(filePath)
  disp(sprintf('(mlrAnatDBPut) Copying %s to %s',filePath{iFile},destDir));
  success = copyfile(filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile})),'f');
  if ~success,return,end
  % get where it was copied to
  toPath{iFile} = fullfile(fileDir,getLastDir(filePath{iFile}));
end

% now do the commit part
if largefiles
  % add/commit/push to large files repo
  cd(localRepoLargeFiles);
  comments = addCommitPush(toPath,largefiles,comments,verbose);
else
  % add/commit/push to the repo
  cd(localRepo);
  addCommitPush(toPath,false,comments,verbose);
end

cd(curpwd);

%%%%%%%%%%%%%%%%%%%%%%%
%    addCommitPush    %
%%%%%%%%%%%%%%%%%%%%%%%
function comments = addCommitPush(toPath,largefiles,comments,verbose)

% commit to repo
if largefiles
  disppercent(-inf,sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s. This may take a minute or two...',pwd));
else
  disp(sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s',pwd));
end

% add file to repo
for iFile = 1:length(toPath)
  if largefiles
    [status,result] = system(sprintf('hg add %s',toPath{iFile}));
  else
    [status,result] = system(sprintf('hg add %s --large',toPath{iFile}));
  end

  if status ~= 0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not add files to local repo. Have you setup your config file for hg?'));
    return
  end
end

% if there are no comments then ask user for comments
if isempty(comments)
  options.Resize='on';
  comments = inputdlg('Enter commit comments','Commit comments',[1 100],{''},options);
  if ~isempty(comments)
    comments = comments{1};
  end
end

[status,result] = system(sprintf('hg commit -m ''%s''',comments));
if largefiles,disppercent(inf);,end

% and push
if ~verbose || isequal(questdlg(sprintf('Do you want to push to the central repo: %s? This can take several minutes depending on your connection. You will be able to work while this occurs (it will push in the background), but you should not shut off your matlab/computer. If you choose no now, you will need to push manually later.',mrGetPref('mlrAnatDBCentralRepo')),'Do push?','Yes','No','Yes'),'Yes')
  disp(sprintf('(mlrAnatDBAddCommitPush) Pushing repo %s in the background',pwd));
  system(sprintf('hg push --new-branch &'));
end


%%%%%%%%%%%%%%%%%%%%%%
%    getFilenames    %
%%%%%%%%%%%%%%%%%%%%%%
function filePathOut = getFilenames(filePath,matchStr)

%getfilenames that match string
matchStr = cellArray(matchStr);

% go through and find matching files
filePathOut = {};
for iMatch = 1:length(matchStr)
  dirPath = dir(fullfile(filePath,matchStr{iMatch}));
  for iFile = 1:length(dirPath)
    filePathOut{end+1} = fullfile(filePath,dirPath(iFile).name);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    formatSubjectID    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function subjectID = formatSubjectID(subjectID)

if isstr(subjectID)
  if (length(subjectID) > 1) 
    [idLocStart idLocEnd] = regexp(subjectID,'s\d+');
    id = 0;
    if ~isempty(idLocStart)
      id = str2num(subjectID(idLocStart+1:idLocEnd));
    else
      [idLocStart idLocEnd] = regexp(subjectID,'\d+');
      if ~isempty(idLocStart)
	id = str2num(subjectID(idLocStart:idLocEnd));
      end
    end
    % format subject ID
    subjectID = sprintf('s%04i',id);
  end
else
  subjectID = sprintf('s%04i',subjectID);
end

