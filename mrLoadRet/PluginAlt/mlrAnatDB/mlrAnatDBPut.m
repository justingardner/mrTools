% mlrAnatDBPut.m
%
%        $Id:$ 
%      usage: mlrAnatDBPut(subjectID,filePath,fileType,<verbose=0>,<comments=[]>)
%         by: justin gardner
%       date: 06/22/15
%    purpose: Puts files into repo. The filePath should point to a directory outside of mlrAnatDB and
%             that file or directory will be copied into an appropriate place into mlrAnatDBm added to the
%             repo, committed and pushed (if the user accepts)
%
%             Typically you specify what kind of file you have (roi, freesurfer, baseAnat, localizer) and
%             this will put it into the correct directory
% 
%             e.g. To put mlr surfaces from a view
%             mlrAnatDBPut(25,v,'mlrBaseAnat');
%
%             e.g. To put a directory of rois
%             mlrAnatDBPut(25,'~/data/rois','roi');
%
%             e.g. To put a freesurfer directory
%             mlrAnatDBPut(25,'~/data/freesurfer','freesurfer');
%
%             Valid types can be found in the code below: roi, freesurfer, baseAnat, freesurfer, canonical, localizer
%
function tf = mlrAnatDBPut(subjectID,filePath,fileType,varargin)

tf = false;
% check arguments
if nargin < 2
  help mlrAnatDBPut;
  return;
end

% format the subject id
subjectID = mlrAnatDBSubjectID(subjectID);

% and get arguments
getArgs(varargin,{'fileDir=[]','largefiles=[]','verbose=0','comments=[]','freesurfer=[]','useHardLinks=1'});

% check file
if ~isview(filePath) && (~isfile(filePath) && ~isdir(filePath))
  disp(sprintf('(mlrAnatDBPut) Could not find: %s',filePath));
  return
end

% figure out what path to put it under and whether it is a large file or not
if isempty(fileDir)
  switch lower(fileType)
   case {'roi','rois'}
    fileDir = 'mlrROIs';
    if isempty(largefiles),largefiles = false;end
    if isview(filePath)
      % if a view then have user select which bases to get
      filePath = getROIFilenames(filePath,subjectID);
    elseif isdir(filePath)
      filePath = getFilenames(filePath,'*.mat');
    end
   case 'freesurfer'
    fileDir = 'anatomy';
    if isempty(largefiles),largefiles = true;end
   case {'3d','canonical','anatomy','anat'}
    fileDir = 'anatomy';
    if isempty(largefiles),largefiles = false;end
   case {'baseanat','mlrbaseanat','mlrbaseanatomy','mlrbaseanatomies'}
    fileDir = 'mlrBaseAnatomies';
    if isempty(largefiles),largefiles = false;end
    if isview(filePath)
      % if a view then have user select which bases to get
      filePath = getBaseFilenames(filePath);
    elseif isdir(filePath)
      filePath = getFilenames(filePath,'*.mat');
    end
   case {'surfaces','surface'}
    fileDir = 'surfaces';
    if isempty(largefiles),largefiles = false;end
    if isdir(filePath),filePath = getFilenames(filePath,{'*.off','*.vff','*.hdr','*.img','*.nii'});end
   case {'localizer','localizers'}
    fileDir = 'localizers';
    if isempty(largefiles),largefiles = true;end
  end    
end
if isempty(filePath),return,end
if isempty(fileDir)
  disp(sprintf('(mlrAnatDBPut) Must specify a vaild fileType or fileDir'));
  return
end

% if there are no comments then ask user for comments
if isempty(comments)
  options.Resize='on';
  comments = inputdlg('Enter commit comments','Commit comments',[1 100],{''},options);
  if ~isempty(comments)
    comments = comments{1};
  else
    comments = '';
  end
end

% get local repo
if ~largefiles
  [localRepo] = mlrAnatDBGetRepo(subjectID);
  if isempty(localRepo),return,end
else
  [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID);
  if isempty(localRepo)||isempty(localRepoLargeFiles),return,end
end

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

% hard link or copy each file using force
for iFile = 1:length(filePath)
  if useHardLinks
    disp(sprintf('(mlrAnatDBPut) Hard linking %s to %s',filePath{iFile},destDir));
    if isdir(filePath{iFile})
      % rsync if a directory
      [status,result] = system(sprintf('rsync -a --link-dest=%s %s/ %s',filePath{iFile},filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile}))));
    else
      % link if a file
      [status,result] = system(sprintf('ln -f %s %s',filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile}))));
    end
    % check success
    if status
      disp(sprintf('(mlrAnatDBPut) Could not link/copy file to repo: %s',result));
      return
    end
  else
    disp(sprintf('(mlrAnatDBPut) Copying %s to %s',filePath{iFile},destDir));
    success = copyfile(filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile})),'f');
    % check success
    if ~success,return,end
  end
  % get where it was copied to
  toPath{iFile} = fullfile(fileDir,getLastDir(filePath{iFile}));
end

% make links to canonical for surfaces
if strcmp(lower(fileType),'surfaces')
  % get list of nifti files in directory - these should be canonical
  canonicals = getFilenames(destDir,{'*.hdr','*.img','*.nii'});
  cd(localRepo);
  % make links to canonical
  for iFile = 1:length(canonicals)
    linkFrom = fullfile('surfaces',getLastDir(canonicals{iFile}));
    toPath{end+1} = setext(subjectID,getext(canonicals{iFile}));
    mysystem(sprintf('ln -sfh %s %s',linkFrom,toPath{end}));
  end
  % make a link to freesurfer dir(this will not be tracked by 
  % hg because thereis an .hgignore file already there to hide it)
  cd('surfaces');
  linkFrom = fullfile('..','anatomy',freesurfer);
  mysystem(sprintf('ln -sfh %s freesurfer',linkFrom));
  % and also make a text file that contains the freesurfer link
  if isfile('.freesurfer')
    mysystem(sprintf('rm -f .freesurfer'));
  end
  % store what the file links to in the .freesurfer file
  % this is so the mlrAnatDBGetRepo can remake the link as needed
  % since hg will not store a link to outside the repo
  mysystem(sprintf('echo %s > .freesurfer',fullfile('anatomy',freesurfer)));
  toPath{end+1} = fullfile('surfaces','.freesurfer');
end

% if this is a freesurfer directory, then go look for surfRelax
% so that we can remove it and then add it back later to the surfaces dir
if strcmp(lower(fileType),'freesurfer')
  surfRelaxDir = fullfile(localRepoLargeFiles,toPath{1},'surfRelax');
  if isdir(surfRelaxDir)
    % remove it (we are going to add it back later to surfaces)
    rmdir(surfRelaxDir,'s');
    surfRelaxDir = fullfile(filePath{1},'surfRelax');
  else
    surfRelaxDir = '';
  end
  % update the branch numbers for both the small and large files
  branchNum = mlrAnatDBGetBranchNum(localRepo);
  mlrAnatDBSetBranchNum(localRepo,branchNum+1);
  mlrAnatDBSetBranchNum(localRepoLargeFiles,branchNum+1);
end

% now do the commit part
if largefiles
  % add and commit to large files repo
  cd(localRepoLargeFiles);
  comments = addCommit(toPath,largefiles,comments,verbose);
else
  % add and commit to the repo
  cd(localRepo);
  addCommit(toPath,false,comments,verbose);
end

cd(curpwd);

% ask if the user wants to push
if ~verbose || isequal(questdlg(sprintf('Do you want to push to the central repo: %s? This can take several minutes depending on your connection. If you choose no now, you will need to push using mlrAnatDBPush later.',mrGetPref('mlrAnatDBCentralRepo')),'Do push?','Yes','No','Yes'),'Yes')
  mlrAnatDBPush(subjectID);
end

% for freesurfer, go and add the surfRelax directory as surfaces
if strcmp(lower(fileType),'freesurfer') && ~isempty(surfRelaxDir)
  tf = mlrAnatDBPut(subjectID,surfRelaxDir,'surfaces','comments',comments,'freesurfer',filePath{1});
  if ~tf,return,end
end

% success.
tf = true;

%%%%%%%%%%%%%%%%%%%
%    addCommit    %
%%%%%%%%%%%%%%%%%%%
function comments = addCommit(toPath,largefiles,comments,verbose)

% commit to repo
if largefiles
  disppercent(-inf,sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s. This may take a minute or two...',pwd));
else
  disp(sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s',pwd));
end

% add file to repo
for iFile = 1:length(toPath)
  if largefiles
    [status,result] = mysystem(sprintf('hg add %s --large',toPath{iFile}));
  else
    [status,result] = mysystem(sprintf('hg add %s',toPath{iFile}));
  end

  if status ~= 0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not add files to local repo. Have you setup your config file for hg?'));
    return
  end
end

% commit
[status,result] = mysystem(sprintf('hg commit -m ''%s''',comments));
if largefiles,disppercent(inf);,end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getBaseFilenames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function filePath = getBaseFilenames(v)

% get list of ROIs to save
baseList = selectInList(v,'bases','Select MLR Base Anatomies to save');

% list of extensions that might be generated
extList = {'hdr','img','nii','mat'};

% default to no files
filePath = {};

% save each anat locally, getting the filenames as the ones to put into the repo
for iBase = baseList
  thisAnat = saveAnat(v,iBase,false,false);
  [thisAnatPath thisAnatFilename] = fileparts(thisAnat);
  % get all the files associated with thisAnat
  for iExt = 1:length(extList)
    filename = fullfile(thisAnatPath,setext(thisAnatFilename,extList{iExt}));
    if isfile(filename)
      filePath{end+1} = filename;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getROIFilenames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function filePath = getROIFilenames(v,subjectID)

% default to no files
filePath = {};

% get list of ROIs to save
roiList = selectInList(v,'rois','Select ROI(s) to save');
if isempty(roiList),return,end

% get the branch number
branchNum = mlrAnatDBGetBranchNum(mlrAnatDBGetRepo(subjectID));

% session name
sessionName = getLastDir(viewGet(v,'homeDir'));
baseName = viewGet(v,'baseName');

% get username
[status,userName] = system('hg config ui.username');
userName = strtrim(userName);

% and save them
for iRoi = roiList
  % get the roi
  roi = viewGet(v,'roi',iRoi);
  % set user who is saving this
  v = viewSet(v,'roiCreatedBy',userName,iRoi);
  % set subjectID
  v = viewSet(v,'roiSubjectID',subjectID,iRoi);
  % set session if not already
  if isempty(roi.createdFromSession)
    % set that the ROI is created from this session
    v = viewSet(v,'roiCreatedFromSession',sessionName,iRoi);
  end
  % if createdOnBase is not set
  if isempty(roi.createdOnBase)
    % set to current base name
    v = viewSet(v,'roiCreatedOnBase',baseName,iRoi);
  end
  % if displayOnBase is not set
  if isempty(roi.displayOnBase)
    % set to current base name
    v = viewSet(v,'roiDisplayOnBase',baseName,iRoi);
  end
  % set the branch num of repo
  if isempty(roi.branchNum)
    % set that the ROI is created from this session
    v = viewSet(v,'roiBranchNum',branchNum+1,iRoi);
  end
  % save it  
  filePath{end+1} = saveROI(v,iRoi,false);
  % set the .mat extension
  filePath{end} = fullfile(fileparts(filePath{end}),setext(getLastDir(filePath{end}),'mat'));
end

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPut): %s',command));
[status,result] = system(command,'-echo');

