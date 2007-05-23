function [val val2] = viewGet(view,param,varargin)
%
%   val = viewGet(view,param,varargin)
%
% Read the parameters of a view (inplane, gray, or flat).
% Access to these structures should go through this routine and through
% viewSet.
%
% view can be either a view structure or a viewNum which is interpreted as
% MLR.views{viewNum}. For some params (e.g., 'subject', 'groupNum'), view
% can be [].
%
% For some params, you must pass additional arguments in varargin. For
% some params, the additional arguments are optional. An example that
% illustrates all of these possibilities is:
%       n = viewGet(view,'nScans',[groupNum]);
% where view can be either view structure or []. If view is a view
% structure, then and groupNum is optional (defaults to the current group).
%
% Examples:
%
% tseriesdir = viewGet(view,'tseriesdir');
% n = viewGet(view,'numberofGroups');
% groupNames = viewGet([],'groupNames');
% groupNum = viewGet(view,'currentGroup');
% n = viewGet(view,'nScans',groupNum);
% n = viewGet([],'nScans',groupNum);
%
% See comments in the source code for the full set of parameters and for
% optional arguments.
%
% 6/2004 djh
% 11/2006 jlg added help and
% originalScanNum,originalFileName,originalGroupName, scanNum,
% dicomName, dicom, stimFileName, stimFile, concatInfo, transforms, TR
mrGlobals

% if ieNotDefined('view'), error('No view defined.'); end
if ieNotDefined('param')
  dispViewGetHelp;
  return
end

% If 'view' is a viewNum, find the corresponding view structure
if isnumeric(view) & ~isempty(view)
  view = MLR.views{view};
end

% Initialize return value
val = [];val2 = [];

switch lower(param)

  case {'view'}
    % view = viewGet(view,'view')
    val = view;
  case {'viewtype','type'}
    % viewType = viewGet(view,'viewType')
    val = view.viewType;
  case {'viewnums'}
    % n = viewGet([],'viewnums')
    for v=1:length(MLR.views)
      if isview(MLR.views{v})
        val = [val,v];
      end
    end
  case {'viewnum'}
    % viewNum = viewGet(view,'viewNum')
    val = view.viewNum;
  case {'subject'}
    % subject = viewGet(view,'subject')
    val = MLR.session.subject;

    % subdirectories
  case {'homedir','homedirectory','sessiondirectory'}
    % homeDir = viewGet(view,'homeDir')
    val = MLR.homeDir;
  case {'anatomydir'}
    % anatomyDir = viewGet(view,'anatomydir')
    homeDir = viewGet(view,'homedir');
    val = fullfile(homeDir,'Anatomy');
    if ~exist(val,'dir')
      mkdir(homeDir,'Anatomy');
    end
  case {'datadir'}
    % dataDir = viewGet(view,'datadir',[groupNum])
    % dataDir = viewGet(view,'datadir',[])
    % dataDir = viewGet(view,'datadir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      viewDir = viewGet(view,'homedir');
      groupname = viewGet(view,'groupname',groupNum);
      val = fullfile(viewDir,groupname);
      if ~exist(val,'dir')
        mkdir(viewDir,groupname);
      end
    end
  case {'roidir'}
    % roiDir = viewGet(view,'roidir')
    viewDir = viewGet(view,'homedir');
    val = fullfile(viewDir,'ROIs');
    if ~exist(val,'dir')
      mkdir(viewDir,'ROIs');
    end
  case {'etcdir'}
    % etcDir = viewGet(view,'etcdir')
    etcDir = fullfile(fileparts(viewGet(view,'datadir',1)),'Etc');
    if isdir(etcDir)
      val = etcDir;
    end
  case {'tseriesdir'}
    % tseriesdir = viewGet(view,'tseriesdir',[groupNum])
    % tseriesdir = viewGet(view,'tseriesdir',[])
    % tseriesdir = viewGet(view,'tseriesdir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      dataDir = viewGet(view,'datadir',groupNum);
      val = fullfile(dataDir,'TSeries');
      if ~exist(val,'dir')
        mkdir(dataDir,'TSeries');
      end
    end
  case {'analysisdir'}
    % analysisdir = viewGet(view,'analysisdir',[groupNum],[analysisNum])
    % analysisdir = viewGet(view,'analysisdir',[groupNum],[])
    % analysisdir = viewGet(view,'analysisdir',[],[analysisNum])
    % analysisdir = viewGet(view,'analysisdir',[],[])
    % analysisdir = viewGet(view,'analysisdir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
      analysisNum = viewGet(view,'currentAnalysis');
    end
    switch length(varargin)
      case 1
        groupNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        groupNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = view.analyses{analysisNum};
    if ~isempty(groupNum) & ~isempty(analysis)
      dataDir = viewGet(view,'datadir',groupNum);
      val = fullfile(dataDir,analysis.type);
      if ~exist(val,'dir')
        mkdir(dataDir,analysis.type);
      end
    end
  case {'overlaydir'}
    % overlaydir = viewGet(view,'overlaydir',[groupNum],[analysisNum])
    val = viewGet(view,'analysisdir',varargin{:});

    % group
  case{'numberofgroups','numGroups','nGroups'}
    % n = viewGet(view,'numberofGroups')
    val = length(MLR.groups);
  case {'groupnames'}
    % groupNames = viewGet(view,'groupNames')
    val = {MLR.groups(:).name};
  case{'currentgroup','curgroup','selectedgroup'}
    % groupNum = viewGet(view,'currentGroup')
    val = view.curGroup;
  case{'group'}
    % groupStructure = viewGet(view,'group',[groupNum])
    % groupStructure = viewGet([],'group',groupNum)
    % groupStructure = viewGet(view,'group',[])
    % groupStructure = viewGet(view,'group')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    val = MLR.groups(groupNum);
  case{'groupname'}
    % groupName = viewGet(view,'groupName',[groupNum])
    % groupName = viewGet([],'groupName',groupNum)
    % groupName = viewGet(view,'groupName',[])
    % groupName = viewGet(view,'groupName')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    val = MLR.groups(groupNum).name;
  case{'groupnum'}
    % groupnum = viewGet(view,'groupnum',groupName)
    groupName = varargin{1};
    groupNames = {MLR.groups(:).name};
    val = find(strcmp(groupName,groupNames));

    % scan
  case{'nscans','numberofscans','numscans'}
    % n = viewGet(view,'nScans',[groupNum])
    % n = viewGet([],'nScans',groupNum)
    % n = viewGet(view,'nScans',[])
    % n = viewGet(view,'nScans')
    if ieNotDefined('varargin')
      g = viewGet(view,'currentGroup');
    else
      g = varargin{1};
    end
    if isempty(g)
      g = viewGet(view,'currentGroup');
    end
    val = length(MLR.groups(g).scanParams);
  case{'scanparams'}
    % n = viewGet(view,'scanParams',scanNum,[groupNum])
    % n = viewGet([],'scanParams',scanNum,groupNum)
    % n = viewGet(view,'scanParams',scanNum,[])
    % n = viewGet(view,'scanParams',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s);
    end
  case{'auxparams'}
    % n = viewGet(view,'auxParams',scanNum,[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).auxParams(s);
    end
  case{'totalframes'}
    % n = viewGet(view,'totalFrames',scanNum,[groupNum])
    % n = viewGet([],'totalFrames',scanNum,groupNum)
    % n = viewGet(view,'totalFrames',scanNum,[])
    % n = viewGet(view,'totalFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).totalFrames;
    end
  case{'junkframestotal'}
    % gets all junk frames from this scan and the scans this was made from
    % n = viewGet(view,'junkFramesTotal',scanNum,[groupNum])
    % n = viewGet([],'junkFramesTotal',scanNum,groupNum)
    % n = viewGet(view,'junkFramesTotal',scanNum,[])
    % n = viewGet(view,'junkFramesTotal',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    % get the junk frame from this scan
    junkFramesThisScan = viewGet(view,'junkFrames',s,g);
    if isempty(junkFramesThisScan),junkFramesThisScan=0;,end
    % find the original scan num
    [os og] = viewGet(view,'originalScanNum',s,g);
    if isempty(os)
      junkFramesTotal = junkFramesThisScan;
    else
      for osNum = 1:length(os)
        % get the junk frames from the original scan
        junkFramesOriginal = viewGet(view,'junkFramesTotal',os(osNum),og(osNum));
        % and add that to the total
        if ~isempty(junkFramesOriginal)
          junkFramesTotal(osNum) = junkFramesOriginal+junkFramesThisScan;
        else
          junkFramesTotal(osNum) = junkFramesThisScan;
        end
      end
    end
    val = junkFramesTotal;
  case{'junkframes'}
    % n = viewGet(view,'junkFrames',scanNum,[groupNum])
    % n = viewGet([],'junkFrames',scanNum,groupNum)
    % n = viewGet(view,'junkFrames',scanNum,[])
    % n = viewGet(view,'junkFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).junkFrames;
    end
  case{'nframes'}
    % n = viewGet(view,'nFrames',scanNum,[groupNum])
    % n = viewGet([],'nFrames',scanNum,groupNum)
    % n = viewGet(view,'nFrames',scanNum,[])
    % n = viewGet(view,'nFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).nFrames;
    end
  case{'frameperiod'}
    % framePeriod = viewGet(view,'framePeriod',scanNum,[groupNum])
    % framePeriod = viewGet([],'framePeriod',scanNum,groupNum)
    % framePeriod = viewGet(view,'framePeriod',scanNum,[])
    % framePeriod = viewGet(view,'framePeriod',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    if isempty(g)
      g = viewGet(view,'currentGroup');
    end
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).framePeriod;
    end
  case{'tseriespathstr','tseriespath'}
    % tseriesPath = viewGet(view,'tseriesPath',scanNum,[groupNum])
    % tseriesPath = viewGet([],'tseriesPath',scanNum,groupNum)
    % tseriesPath = viewGet(view,'tseriesPath',scanNum,[])
    % tseriesPath = viewGet(view,'tseriesPath',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = fullfile(viewGet(view,'tseriesDir',g),...
        MLR.groups(g).scanParams(s).fileName);
    end
  case{'tseriesfile'}
    % filename = viewGet(view,'tseriesFile',scanNum,groupNum)
    % filename = viewGet([],'tseriesFile',scanNum,groupNum)
    [s g] = getScanAndGroup(view,varargin,param);
    if s <= length(MLR.groups(g).scanParams)
      val = MLR.groups(g).scanParams(s).fileName;
    else
      val = [];
    end
  case {'scannum'}
    % get scanNum (and groupNum) of a scan by its file name
    % [scanNum,groupNum] = viewGet(view,'scanNum',tseriesFileName);
    % returns empty if the filename does not exist
    %
    % if you pass in a groupNum as well, then search for
    % tfilename is only done w/in the group.
    % [scanNum,groupNum] = viewGet(view,'scanNum',tseriesFileName,groupNum);
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet scanNum: Must specify tseriesFileName.');
    end
    [tseriesFilePath,tseriesFileName] = fileparts(varargin{1});
    if length(varargin) > 1
      groupNums = varargin{2};
    else
      groupNums = 1:length(MLR.groups);
    end
    scanNumMatch = [];groupNumMatch = [];
    for groupNum = groupNums
      for scanNum = 1:length(MLR.groups(groupNum).scanParams)
        [scanFilePath,scanFileName] = fileparts(MLR.groups(groupNum).scanParams(scanNum).fileName);
        if strcmp(scanFileName,tseriesFileName)
          scanNumMatch(end+1) = scanNum;
          groupNumMatch(end+1) = groupNum;
        end
      end
    end
    val = scanNumMatch;
    val2 = groupNumMatch;
  case {'concatinfo'}
    % concatInfo = viewGet(view,'stimfile',scanNum,[groupNum]);
    if ieNotDefined('varargin')
      mrErrorDlg(sprintf('viewGet %s: must specify scan.',param));
    end
    s = varargin{1};
    if length(varargin) > 1,g = varargin{2};else,g = viewGet(view,'currentGroup');,end
    if isempty(g),g = viewGet(view,'currentGroup');end
    [tseriesPath,tseriesFile] = fileparts(viewGet(view,'tseriesPath',s,g));
    % check for mat file
    matFileName = fullfile(tseriesPath,sprintf('%s.mat',tseriesFile));
    if isfile(matFileName)
      load(matFileName);
      if exist('concatInfo')
        val = concatInfo;
      end
    end
  case {'stimfile'}
    % stimfile = viewGet(view,'stimfile',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    stimFileName = viewGet(view,'stimFileName',s,g);
    val = {};
    % cycle over all stimfiles (concatenations and averages,
    % may have several stimfiles).
    if ~isempty(stimFileName)
      for j = 1:length(stimFileName)
        % load this stimfile
        if ~isfile(stimFileName{j})
          mrErrorDlg(sprintf('viewGet %s: Could not find stimfile %s',param,stimFileName{j}));
        else
          val{j}=load(stimFileName{j});
        end
        % check to see what type it is, and set the field appropriately
        if isfield(val{j},'mylog')
          val{j}.filetype = 'eventtimes';
        elseif isfield(val{j},'myscreen')
          val{j}.filetype = 'mgl';
        end
      end
    end
  case {'stimfilename'}
    % stimFileName = viewGet(view,'stimFileName',scanNum,[groupNum]);
    % stimFileName is returned as a cell array of all stimfiles
    % associated with the scan
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    stimFileName{1} = '';
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),'stimFileName')
        if ~isempty(MLR.groups(g).auxParams(s).stimFileName)
          stimFileName{1} = fullfile(viewGet(view,'etcDir'),MLR.groups(g).auxParams(s).stimFileName);
        end
      end
    end
    % if the file does not exist, then check original
    if isempty(stimFileName{1})
      stimFileName = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          newStimFileNames = viewGet(view,'stimFileName',os(osnum),og(osnum));
          for j = 1:length(newStimFileNames)
            stimFileName{end+1} = newStimFileNames{j};
          end
        end
      end
    end
    val = stimFileName;
  case {'tr'}
    % tr = viewGet(view,'tr',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    % first get dicom
    dicom = viewGet(view,'dicom',s,g);
    val = [];
    if isempty(dicom)
      return
    end
    for i = 1:length(dicom)
      if isfield(dicom{i},'ACQ') && isfield(dicom{i}.ACQ,'Repetition_Time')
        thistr = dicom{i}.ACQ.Repetition_Time/1000;
        if (isempty(val))
          val = thistr;
        elseif (val ~= thistr)
          disp(sprintf('(viewGet) Data collected with different TR: %0.2f vs %0.2f',val,thistr));
        end
      end
    end
  case {'dicom'}
    % dicom = viewGet(view,'dicom',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    dicomName = viewGet(view,'dicomName',s,g);
    val = {};
    if ~isempty(dicomName)
      for j = 1:length(dicomName)
        val{j} = readdicomheader(dicomName{j});
      end
    end
  case {'dicomname'}
    % dicomName = viewGet(view,'dicomName',scanNum,[groupNum]);
    % dicomName is returned as a cell array of all dicom headers
    % associated with the scan
    [s g] = getScanAndGroup(view,varargin,param);
    % create name of file, by putting -header.txt on to the filename
    [tseriesPath,tseriesFile] = fileparts(viewGet(view,'tseriesPath',s,g));
    dicomName{1} = fullfile(tseriesPath,sprintf('%s-header.txt',tseriesFile));
    % if the file does not exist, then check original
    if ~isfile(dicomName{1})
      dicomName = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          newDicomNames = viewGet(view,'dicomName',os(osnum),og(osnum));
          for j = 1:length(newDicomNames)
            dicomName{end+1} = newDicomNames{j};
          end
        end
      end
    end
    val = dicomName;
  case {'originalscannum'}
    % [scanNum,groupNum] = viewGet(view,'originalScanNum',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    originalFileNames = viewGet(view,'originalFileName',s,g);
    originalGroupNames = viewGet(view,'originalGroupName',s,g);
    originalScanNum = [];originalGroupNum = [];
    for i = 1:length(originalFileNames)
      % make sure there is a valid scan to return
      if ~isempty(viewGet(view,'scanNum',originalFileNames{i},viewGet(view,'groupNum',originalGroupNames{i})))
        [originalScanNum(end+1),originalGroupNum(end+1)] = viewGet(view,'scanNum',originalFileNames{i},viewGet(view,'groupNum',originalGroupNames{i}));
      end
    end
    val = originalScanNum;val2 = originalGroupNum;
  case {'originalfilename'}
    % originalFileName = viewGet(view,'originalFileName',scanNum,[groupNum]);
    % if the original is the same as the filename, then this returns []
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isfield(MLR.groups(g).scanParams(s),'originalFileName')
        val = MLR.groups(g).scanParams(s).originalFileName;
        if strcmp(val,MLR.groups(g).scanParams(s).fileName)
          val = [];
        end
      else
        val = [];
      end
    end
  case {'originalgroupname'}
    % originalGroupName = viewGet(view,'originalGroupName',scanNum,[groupNum]);
    % if the original is the same as the filename, then this returns []
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isfield(MLR.groups(g).scanParams(s),'originalGroupName')
        val = MLR.groups(g).scanParams(s).originalFileName;
        if strcmp(val,MLR.groups(g).scanParams(s).fileName)
          val = [];
        else
          val = MLR.groups(g).scanParams(s).originalGroupName;
        end
      else
        val = [];
      end
    end
  case{'transforms'}
    % transforms = viewGet(view,'transforms',scanNum,[groupNum]);
    % returns motion correction transformation matrices if they exists
    [s g] = getScanAndGroup(view,varargin,param);
    % get the tseries name
    [tSeriesPath tSeriesName] = fileparts(viewGet(view,'tSeriesPathStr',s,g));
    if isfile(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)))
      load(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)));
      if exist('transforms') == 1
        val = transforms;
      end
    end
  case{'params'}
    % transforms = viewGet(view,'transforms',scanNum,[groupNum]);
    % gets the .mat params file associated with this scan if it exists
    [s g] = getScanAndGroup(view,varargin,param);
    % get the tseries name
    [tSeriesPath tSeriesName] = fileparts(viewGet(view,'tSeriesPathStr',s,g));
    if isfile(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)))
      val = load(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)));
    end
  case{'niftihdr'}
    % hdr = viewGet(view,'niftiHdr',scanNum,[groupNum])
    % hdr = viewGet([],'niftiHdr',scanNum,groupNum)
    % hdr = viewGet(view,'niftiHdr',scanNum,[])
    % hdr = viewGet(view,'niftiHdr',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % make sure it is a nifti headr (not an analyze)
      val = MLR.groups(g).scanParams(s).niftiHdr;
    end
  case{'sformcode'}
    % xform = viewGet(view,'sformcode',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'sform_code')
        val = [];
      else
        val = MLR.groups(g).scanParams(s).niftiHdr.sform_code;
      end
    end
  case{'scanxform'}
    % xform = viewGet(view,'scanXform',scanNum,[groupNum])
    % xform = viewGet([],'scanXform',scanNum,groupNum)
    % xform = viewGet(view,'scanXform',scanNum,[])
    % xform = viewGet(view,'scanXform',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % make sure it is a nifti headr (not an analyze)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'sform_code')
        return
      end
      sformCode = MLR.groups(g).scanParams(s).niftiHdr.sform_code;
      if sformCode
        % if sform has been set, then use it.
        val = MLR.groups(g).scanParams(s).niftiHdr.sform44;
      else
        % If has not been set, then use the transform that
        % transforms this image directly on to the current anatomy
        % using the qform matrices. Note: There used to be code
        % here that reset val if it was the identity but that was a
        % bug (DJH 1/17/07).
        if strcmp(mrGetPref('verbose'),'Yes')
          disp('(viewGet) sform is not set. Using qform to align to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this');
        end
        baseqform = viewGet(view,'baseqform');
        val = pinv(baseqform)*MLR.groups(g).scanParams(s).niftiHdr.qform44;
      end
    end
  case{'scanvoxelsize'}
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum,[groupNum])
    % voxelSize = viewGet([],'scanVoxelSize',scanNum,groupNum)
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum,[])
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).voxelSize;
    end
  case{'scandims'}
    % dims = viewGet(view,'scanDims',scanNum,[groupNum])
    % dims = viewGet([],'scanDims',scanNum,groupNum)
    % dims = viewGet(view,'scanDims',scanNum,[])
    % dims = viewGet(view,'scanDims',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).dataSize;
    end
  case{'description'}
    % string = viewGet(view,'description',scanNum,[groupNum])
    % string = viewGet([],'description',scanNum,groupNum)
    % string = viewGet(view,'description',scanNum,[])
    % string = viewGet(view,'description',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).description;
    end

    % baseVolume (anatomy)
  case{'numberofbasevolumes'}
    % n = viewGet(view,'numberofBaseVolumes')
    val = length(view.baseVolumes);
  case {'curbase','currentbase'}
    % baseNum = viewGet(view,'currentBase')
    val = view.curBase;
  case{'basenum'}
    % baseNum = viewGet(view,'baseNum',baseName)
    baseName = varargin{1};
    baseNames = {view.baseVolumes(:).name};
    val = find(strcmp(baseName,baseNames));
  case{'basevolume','baseanatomy'}
    % basevolume = viewGet(view,'baseVolume',[baseNum])
    % basevolume = viewGet(view,'baseVolume',[])
    % basevolume = viewGet(view,'baseVolume')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b);
    end
  case {'basenames'}
    % baseNames = viewGet(view,'baseNames')
    if ieNotDefined('varargin')
      baseNum = viewGet(view,'curBase');
    else
      baseNum = varargin{1};
    end
    if isempty(baseNum)
      baseNum = viewGet(view,'curBase');
    end
    if ~isempty(view.baseVolumes)
      val = {view.baseVolumes(:).name};
    end
  case {'basename'}
    % basename = viewGet(view,'basename',[baseNum])
    % basename = viewGet(view,'basename',[])
    % basename = viewGet(view,'basename')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).name;
    end
  case {'basedata'}
    % basedata = viewGet(view,'basedata',[baseNum])
    % basedata = viewGet(view,'basedata',[])
    % basedata = viewGet(view,'basedata')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).data;
    end
  case {'baseclip'}
    % baseclip = viewGet(view,'baseclip',[baseNum])
    % baseclip = viewGet(view,'baseclip',[])
    % baseclip = viewGet(view,'baseclip')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).clip;
    end
  case {'baserange'}
    % baserange = viewGet(view,'baserange',[baseNum])
    % baserange = viewGet(view,'baserange',[])
    % baserange = viewGet(view,'baserange')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).range;
    end
  case {'basedims'}
    % basedims = viewGet(view,'basedims',[baseNum])
    % basedims = viewGet(view,'basedims',[])
    % basedims = viewGet(view,'basedims')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    baseVolume = viewGet(view,'baseVolume',b);
    if ~isempty(baseVolume)
      val = size(baseVolume.data);
    end
  case {'basexform'}
    % basexform = viewGet(view,'basexform',[baseNum])
    % basexform = viewGet(view,'basexform',[])
    % basexform = viewGet(view,'basexform')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    baseVolume = viewGet(view,'baseVolume',b);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.sform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'baseqform'}
    % basexform = viewGet(view,'baseqform',[baseNum])
    % basexform = viewGet(view,'baseqform',[])
    % basexform = viewGet(view,'baseqform')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    baseVolume = viewGet(view,'baseVolume',b);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.qform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'basevolpermutation'}
    % basevolpermutation = viewGet(view,'basevolpermutation',[baseNum])
    % basevolpermutation = viewGet(view,'basevolpermutation',[])
    % basevolpermutation = viewGet(view,'basevolpermutation')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    baseVolume = viewGet(view,'baseVolume',b);
    if ~isempty(baseVolume)
      val = baseVolume.permutationMatrix;
    else
      % placeholder when base images not loaded
      val = eye(3);
    end
  case {'basesliceindex'}
    % basesliceindex = viewGet(view,'basesliceindex',[baseNum])
    % basesliceindex = viewGet(view,'basesliceindex',[])
    % basesliceindex = viewGet(view,'basesliceindex')
    %
    % This admittedly arcane logic is also used in mrAlignGUI. If you
    % change this code, please make corresponding changes in that
    % function. Permutation matrix is set by loadAnat using even more
    % arcane logic.
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    sliceOrientation = viewGet(view,'sliceOrientation');
    permutation = viewGet(view,'baseVolPermutation',b);
    switch sliceOrientation
      case 3   % Sagittal
        [m,index] = max(permutation * [1 0 0]');
      case 2   % Coronal
        [m,index] = max(permutation * [0 1 0]');
      case 1   % Axial
        [m,index] = max(permutation * [0 0 1]');
    end
    val = index;
  case {'basevoxelsize'}
    % basevoxelsize = viewGet(view,'basevoxelsize',[baseNum])
    % basevoxelsize = viewGet(view,'basevoxelsize',[])
    % basevoxelsize = viewGet(view,'basevoxelsize')
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    baseVolume = viewGet(view,'baseVolume',b);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.pixdim([2,3,4])';;
    end

    % ROI
  case {'showrois'}
    % show = viewGet(view,'showROIs')
    val = view.showROIs;
  case{'numberofrois'}
    % n = viewGet(view,'numberofROIs')
    val = length(view.ROIs);
  case{'currentroi','currentroinum'}
    % roiNum = viewGet(view,'currentROI')
    val = view.curROI;
  case{'roinum'}
    % roiNum = viewGet(view,'roiNum',roiName)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet roiNum: must specify ROI name');
    end
    ROIname = varargin{1};
    ROInames = {view.ROIs(:).name};
    val = find(strcmp(ROIname,ROInames));
  case{'roi'}
    % roi = viewGet(view,'roi',[roiNum])
    % roi = viewGet(view,'roi',[])
    % roi = viewGet(view,'roi')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r);
    end
  case {'roinames'}
    % roiNames = viewGet(view,'roiNames')
    if ieNotDefined('varargin')
      roiNum = viewGet(view,'currentROI');
    else
      roiNum = varargin{1};
    end
    if isempty(roiNum)
      roiNum = viewGet(view,'currentROI');
    end
    if ~isempty(view.ROIs)
      val = {view.ROIs(:).name};
    end
  case{'roiname'}
    % roiName = viewGet(view,'roiName',[roiNum])
    % roiName = viewGet(view,'roiName',[])
    % roiName = viewGet(view,'roiName')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).name;
    end
  case{'roicoords'}
    % roiCoords = viewGet(view,'roiCoords',[roiNum])
    % roiCoords = viewGet(view,'roiCoords',[])
    % roiCoords = viewGet(view,'roiCoords')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).coords;
    end
  case{'prevroicoords'}
    % prevRoiCoords = viewGet(view,'prevRoiCoords')
    val = view.prevROIcoords;
  case{'roicolor'}
    % roicolor = viewGet(view,'roicolor',[roiNum])
    % roicolor = viewGet(view,'roicolor',[])
    % roicolor = viewGet(view,'roicolor')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).color;
    end
  case{'roixform'}
    % roixform = viewGet(view,'roixform',[roiNum])
    % roixform = viewGet(view,'roixform',[])
    % roixform = viewGet(view,'roixform')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).xform;
    else
      val = eye(3);
    end
  case{'roivoxelsize'}
    % roivoxelsize = viewGet(view,'roivoxelsize',[roiNum])
    % roivoxelsize = viewGet(view,'roivoxelsize',[])
    % roivoxelsize = viewGet(view,'roivoxelsize')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).voxelSize;
    end
  case('roivolume')
    % roivolume = viewGet(view,'roivolume',[roiNum])
    % roivolume = viewGet(view,'roivolume',[])
    % roivolume = viewGet(view,'roivolume')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      roiVoxelSize = view.ROIs(r).voxelSize;
      roiCoords = view.ROIs(r).coords;
      if ~isempty(roiVoxelSize) & ~isempty(roiCoords)
        val = prod(roiVoxelSize) * size(roiCoords,2);
      end
    end
  case{'roidate'}
    % roidate = viewGet(view,'roidate',[roiNum])
    % roidate = viewGet(view,'roidate',[])
    % roidate = viewGet(view,'roidate')
    if ieNotDefined('varargin')
      r = viewGet(view,'currentROI');
    else
      r = varargin{1};
    end
    if isempty(r)
      r = viewGet(view,'currentROI');
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).date;
    end

    % analysis
  case{'numberofanalyses'}
    % n = viewGet(view,'numberofAnalyses')
    val = length(view.analyses);
  case{'currentanalysis','curanalysis'}
    % n = viewGet(view,'currentAnalysis')
    val = view.curAnalysis;
  case{'analysisnum'}
    % n = viewGet(view,'analysisNum',analysisName)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet analysisNum: must specify analysisName');
    end
    analysisName = varargin{1};
    if ~isempty(view.analyses)
      analysisNames = viewGet(view,'analysisnames');
      val = find(strcmp(analysisName,analysisNames));
    end
  case {'analysis'}
    % analysis = viewGet(view,'analysis',[analysisNum])
    % analysis = viewGet(view,'analysis',[])
    % analysis = viewGet(view,'analysis')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum};
    end
  case {'analysisnames'}
    % analysisNames = viewGet(view,'analysisNames')
    if ~isempty(view.analyses)
      numAnalyses = viewGet(view,'numberofAnalyses');
      names = cell(1,numAnalyses);
      for a = 1:numAnalyses
        names{a} = view.analyses{a}.name;
      end
      val = names;
    end
  case {'analysisname'}
    % analysisName = viewGet(view,'analysisName',[analysisNum])
    % analysisName = viewGet(view,'analysisName',[])
    % analysisName = viewGet(view,'analysisName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.name;
    end
  case {'analysistypes'}
    % analysisTypes = viewGet(view,'analysisTypes')
    if ~isempty(view.analyses)
      numAnalyses = viewGet(view,'numberofAnalyses');
      types = cell(1,numAnalyses);
      for a = 1:numAnalyses
        types{a} = view.analyses{a}.type;
      end
      val = types;
    end
  case {'analysistype'}
    % analysisType = viewGet(view,'analysisType',[analysisNum])
    % analysisType = viewGet(view,'analysisType',[])
    % analysisType = viewGet(view,'analysisType')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.type;
    end
  case {'analysisgroupname'}
    % analysisGroup = viewGet(view,'analysisGroupName',[analysisNum])
    % analysisGroup = viewGet(view,'analysisGroupName',[])
    % analysisGroup = viewGet(view,'analysisGroupName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.groupName;
    end
  case {'analysisfunction'}
    % analysisFunction = viewGet(view,'analysisFunction',[analysisNum]);
    % analysisFunction = viewGet(view,'analysisFunction',[]);
    % analysisFunction = viewGet(view,'analysisFunction');
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.function;
    end
  case {'analysisreconcilefunction','reconcilefunction'}
    % reconcileFunction = viewGet(view,'reconcilefunction',[analysisNum])
    % reconcileFunction = viewGet(view,'reconcilefunction',[])
    % reconcileFunction = viewGet(view,'reconcilefunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.reconcileFunction;
    end
  case {'analysismergefunction','mergefunction'}
    % mergeFunction = viewGet(view,'mergefunction',[analysisNum])
    % mergeFunction = viewGet(view,'mergefunction',[])
    % mergeFunction = viewGet(view,'mergefunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.mergeFunction;
    end
  case {'analysisguifunction'}
    % analysisGui = viewGet(view,'analysisGuiFunction',[analysisNum])
    % analysisGui = viewGet(view,'analysisGuiFunction',[])
    % analysisGui = viewGet(view,'analysisGuiFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.guiFunction;
    end
  case {'analysisparams'}
    % analysisParams = viewGet(view,'analysisParams',[analysisNum])
    % analysisParams = viewGet(view,'analysisParams',[])
    % analysisParams = viewGet(view,'analysisParams')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.params;
    end
  case {'analysisdata'}
    % analysisData = viewGet(view,'analysisData',[analysisNum])
    % analysisData = viewGet(view,'analysisData',[])
    % analysisData = viewGet(view,'analysisData')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.d;
    end
  case {'analysisdate'}
    % analysisdate = viewGet(view,'analysisDate',[analysisNum])
    % analysisdate = viewGet(view,'analysisDate',[])
    % analysisdate = viewGet(view,'analysisDate')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.date;
    end
  case {'overlayinterpfunction'}
    % fnctn = viewGet(view,'overlayInterpFunction',[analysisNum])
    % fnctn = viewGet(view,'overlayInterpFunction',[])
    % fnctn = viewGet(view,'overlayInterpFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if isfield(view.analyses{analysisNum},'overlayInterpFunction');
        val = view.analyses{analysisNum}.overlayInterpFunction;
      else
        val = [];
      end
    end

    % overlay
  case{'overlays'}
    % overlays = viewGet(view,'overlays',[analysisNum])
    % overlays = viewGet(view,'overlays',[])
    % overlays = viewGet(view,'overlays')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = analysis.overlays;
    end
  case{'numberofoverlays'}
    % n = viewGet(view,'numberofOverlays',[analysisNum])
    % n = viewGet(view,'numberofOverlays',[])
    % n = viewGet(view,'numberofOverlays')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = length(analysis.overlays);
    end
  case{'currentoverlay'}
    % n = viewGet(view,'currentOverlay',[analysisNum])
    % n = viewGet(view,'currentOverlay',[])
    % n = viewGet(view,'currentOverlay')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = analysis.curOverlay;
    end
  case{'overlaynum'}
    % n = viewGet(view,'overlayNum',overlayName,[analysisNum])
    % n = viewGet(view,'overlayNum',overlayName,[])
    % n = viewGet(view,'overlayNum',overlayName)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayNum: must specify overlayName.');
    end
    overlayName = varargin{1};
    if (length(varargin) < 2)
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      overlayNames = {analysis.overlays(:).name};
      val = find(strcmp(overlayName,overlayNames));
    end
  case {'overlay'}
    % overlay = viewGet(view,'overlay',[overlayNum],[analysisNum])
    % overlay = viewGet(view,'overlay',overlayNum,[])
    % overlay = viewGet(view,'overlay',[],analysisNum)
    % overlay = viewGet(view,'overlay',[],[])
    % overlay = viewGet(view,'overlay',overlayNum)
    % overlay = viewGet(view,'overlay')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum);
      end
    end
  case {'overlaynames'}
    % overlayNames = viewGet(view,'overlayNames',[analysisNum])
    % overlayNames = viewGet(view,'overlayNames',[])
    % overlayNames = viewGet(view,'overlayNames')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(view.analyses) & ~isempty(analysis.overlays)
      val = {view.analyses{analysisNum}.overlays(:).name};
    end
  case {'overlayname'}
    % overlayName = viewGet(view,'overlayName',[overlayNum],[analysisNum])
    % overlayName = viewGet(view,'overlayName',overlayNum,[])
    % overlayName = viewGet(view,'overlayName',[],analysisNum)
    % overlayName = viewGet(view,'overlayName',[],[])
    % overlayName = viewGet(view,'overlayName',overlayNum)
    % overlayName = viewGet(view,'overlayName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).name;
      end
    end
  case {'overlaytype'}
    % overlayType = viewGet(view,'overlayType',[overlayNum],[analysisNum])
    % overlayType = viewGet(view,'overlayType',overlayNum,[])
    % overlayType = viewGet(view,'overlayType',[],analysisNum)
    % overlayType = viewGet(view,'overlayType',[],[])
    % overlayType = viewGet(view,'overlayType',overlayNum)
    % overlayType = viewGet(view,'overlayType')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).type;
      end
    end
  case{'overlayxform'}
    % xform = viewGet(view,'overlayXform',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayXform: must specify scan.');
    end
    scan = varargin{1};
    val = viewGet(view,'scanxform',scan);
  case {'overlaydata'}
    % overlaydata = viewGet(view,'overlaydata',scanNum,[overlayNum],[analysisNum])
    % overlaydata = viewGet(view,'overlaydata',scanNum,overlayNum,[])
    % overlaydata = viewGet(view,'overlaydata',scanNum,[],analysisNum)
    % overlaydata = viewGet(view,'overlaydata',scanNum,[],[])
    % overlaydata = viewGet(view,'overlaydata',scanNum,overlayNum)
    % overlaydata = viewGet(view,'overlaydata',scanNum)
    % scanNum can be a number or []. If [] then it returns the entire cell
    % array of data.
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayData: must specify scan.');
    end
    scan = varargin{1};
    switch (length(varargin))
      case 1
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay',analysisNum);
      case 2
        overlayNum = varargin{2};
        analysisNum = viewGet(view,'currentAnalysis');
      case 3
        overlayNum = varargin{2};
        analysisNum = varargin{3};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        if isempty(scan)
          val = analysis.overlays(overlayNum).data;
        else
          val = analysis.overlays(overlayNum).data{scan};
        end
      end
    end
  case {'overlaydims'}
    % overlaydims = viewGet(view,'overlaydims',scanNum,[overlayNum],[analysisNum])
    % overlaydims = viewGet(view,'overlaydims',scanNum,overlayNum,[])
    % overlaydims = viewGet(view,'overlaydims',scanNum,[],analysisNum)
    % overlaydims = viewGet(view,'overlaydims',scanNum,[],[])
    % overlaydims = viewGet(view,'overlaydims',scanNum,overlayNum)
    % overlaydims = viewGet(view,'overlaydims',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayDims: must specify scan.');
    end
    scan = varargin{1};
    switch (length(varargin))
      case 1
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay',analysisNum);
      case 2
        overlayNum = varargin{2};
        analysisNum = viewGet(view,'currentAnalysis');
      case 3
        overlayNum = varargin{2};
        analysisNum = varargin{3};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = size(analysis.overlays(overlayNum).data{scan});
      end
    end
  case {'overlayclip'}
    % overlayclip = viewGet(view,'overlayclip',[overlayNum],[analysisNum])
    % overlayclip = viewGet(view,'overlayclip',overlayNum,[])
    % overlayclip = viewGet(view,'overlayclip',[],analysisNum)
    % overlayclip = viewGet(view,'overlayclip',[],[])
    % overlayclip = viewGet(view,'overlayclip',overlayNum)
    % overlayclip = viewGet(view,'overlayclip')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).clip;
      end
    end
  case {'overlaycmap'}
    % overlaycmap = viewGet(view,'overlaycmap',[overlayNum],[analysisNum])
    % overlaycmap = viewGet(view,'overlaycmap',overlayNum,[])
    % overlaycmap = viewGet(view,'overlaycmap',[],analysisNum)
    % overlaycmap = viewGet(view,'overlaycmap',[],[])
    % overlaycmap = viewGet(view,'overlaycmap',overlayNum)
    % overlaycmap = viewGet(view,'overlaycmap')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).colormap;
      end
    end
  case {'overlayctype'}
    % overlayctype = viewGet(view,'overlayctype',[overlayNum],[analysisNum])
    % overlayctype = viewGet(view,'overlayctype',overlayNum,[])
    % overlayctype = viewGet(view,'overlayctype',[],analysisNum)
    % overlayctype = viewGet(view,'overlayctype',[],[])
    % overlayctype = viewGet(view,'overlayctype',overlayNum)
    % overlayctype = viewGet(view,'overlayctype')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        if isfield(analysis.overlays(overlayNum),'colormapType')
          val = analysis.overlays(overlayNum).colormapType;
        else
          val = 'normal';
        end
      end
    end
  case {'overlayrange'}
    % overlayrange = viewGet(view,'overlayrange',[overlayNum],[analysisNum])
    % overlayrange = viewGet(view,'overlayrange',overlayNum,[])
    % overlayrange = viewGet(view,'overlayrange',[],analysisNum)
    % overlayrange = viewGet(view,'overlayrange',[],[])
    % overlayrange = viewGet(view,'overlayrange',overlayNum)
    % overlayrange = viewGet(view,'overlayrange')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).range;
      end
    end
  case {'overlayalpha'}
    % overlayalpha = viewGet(view,'overlayalpha',[overlayNum],[analysisNum])
    % overlayalpha = viewGet(view,'overlayalpha',overlayNum,[])
    % overlayalpha = viewGet(view,'overlayalpha',[],analysisNum)
    % overlayalpha = viewGet(view,'overlayalpha',[],[])
    % overlayalpha = viewGet(view,'overlayalpha',overlayNum)
    % overlayalpha = viewGet(view,'overlayalpha')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).alpha;
      end
    end
  case {'overlaydate'}
    % overlaydate = viewGet(view,'overlaydate',[overlayNum],[analysisNum])
    % overlaydate = viewGet(view,'overlaydate',overlayNum,[])
    % overlaydate = viewGet(view,'overlaydate',[],analysisNum)
    % overlaydate = viewGet(view,'overlaydate',[],[])
    % overlaydate = viewGet(view,'overlaydate',overlayNum)
    % overlaydate = viewGet(view,'overlaydate')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).date;
      end
    end
  case {'interrogator'}
    % interrogator = viewGet(view,'interrogator',[overlayNum],[analysisNum])
    % interrogator = viewGet(view,'interrogator',overlayNum,[])
    % interrogator = viewGet(view,'interrogator',[],analysisNum)
    % interrogator = viewGet(view,'interrogator',[],[])
    % interrogator = viewGet(view,'interrogator',overlayNum)
    % interrogator = viewGet(view,'interrogator')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).interrogator;
      end
    end
    if isempty(val)
      val = 'mrDefaultInterrogator';
    end
  case {'overlaygroupname'}
    % groupName = viewGet(view,'overlayGroupName',[overlayNum],[analysisNum])
    % groupName = viewGet(view,'overlayGroupName',overlayNum,[])
    % groupName = viewGet(view,'overlayGroupName',[],analysisNum)
    % groupName = viewGet(view,'overlayGroupName',[],[])
    % groupName = viewGet(view,'overlayGroupName',overlayNum)
    % groupName = viewGet(view,'overlayGroupName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).groupName;
      end
    end
  case {'overlayparams'}
    % params = viewGet(view,'overlayParams',[overlayNum],[analysisNum])
    % params = viewGet(view,'overlayParams',overlayNum,[])
    % params = viewGet(view,'overlayParams',[],analysisNum)
    % params = viewGet(view,'overlayParams',[],[])
    % params = viewGet(view,'overlayParams',overlayNum)
    % params = viewGet(view,'overlayParams')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).params;
      end
    end
  case {'overlayreconcilefunction'}
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[overlayNum],[analysisNum])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',overlayNum,[])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[],analysisNum)
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[],[])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',overlayNum)
    % reconcileFunction = viewGet(view,'overlayReconcileFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).reconcileFunction;
      end
    end
  case {'overlaymergefunction'}
    % mergeFunction = viewGet(view,'overlayMergeFunction',[overlayNum],[analysisNum])
    % mergeFunction = viewGet(view,'overlayMergeFunction',overlayNum,[])
    % mergeFunction = viewGet(view,'overlayMergeFunction',[],analysisNum)
    % mergeFunction = viewGet(view,'overlayMergeFunction',[],[])
    % mergeFunction = viewGet(view,'overlayMergeFunction',overlayNum)
    % mergeFunction = viewGet(view,'overlayMergeFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).mergeFunction;
      end
    end

    % corAnal
  case {'coranal'}
    % corAnal = viewGet(view,'corAnal')
    % Find all corAnals
    analysisTypes = viewGet(view,'analysisTypes');
    corAnalNums = find(strcmp('corAnal',analysisTypes));
    if ~isempty(corAnalNums)
      % Check if current analysis is corAnal. Otherwise, use the
      % first corAnal on the analyses list.
      curAnalysisNum = viewGet(view,'currentAnalysis');
      if find(curAnalysisNum == corAnalNums)
        n = curAnalysisNum;
      else
        n = corAnalNums(1);
      end
      val = viewGet(view,'analysis',n);
    end
  case {'coranalparams'}
    % corAnalParams = viewGet(view,'corAnalParams')
    corAnal = viewGet(view,'corAnal');
    if ~isempty(corAnal)
      val = corAnal.params;
    end
  case {'co'}
    % co = viewGet(view,'co')
    m = viewGet(view,'analysisNum','corAnal');
    n = viewGet(view,'overlayNum','co',m);
    val = viewGet(view,'overlay',n,m);
  case {'amp'}
    % amp = viewGet(view,'amp')
    m = viewGet(view,'analysisNum','corAnal');
    n = viewGet(view,'overlayNum','amp',m);
    val = viewGet(view,'overlay',n,m);
  case {'ph'}
    % ph = viewGet(view,'ph')
    m = viewGet(view,'analysisNum','corAnal');
    n = viewGet(view,'overlayNum','ph',m);
    val = viewGet(view,'overlay',n,m);
  case {'ncycles'}
    % ncycles = viewGet(view,'ncycles',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet ncycles: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      val = params.ncycles(scanNum);
    end
  case {'detrend'}
    % detrend = viewGet(view,'detrend',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet detrend: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      val = params.detrend{scanNum};
    end
  case {'spatialnorm'}
    % spatialnorm = viewGet(view,'spatialnorm',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet spatialnorm: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      val = params.spatialnorm{scanNum};
    end

    % GUI
  case {'fignum','figurenumber'}
    % figureHandle = viewGet(view,'figureNumber');
    val = view.figure;
  case {'curscan','currentscan'}
    % scan = viewGet(view,'currentScan');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = round(get(handles.scanSlider,'Value'));
  case {'curslice','currentslice'}
    % slice = viewGet(view,'currentSlice');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = round(get(handles.sliceSlider,'Value'));
  case {'curcoords','currentcoordinates'}
    % coords = viewGet(view,'currentCoordinates');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = handles.coords;
  case {'alpha'}
    % alpha = viewGet(view,'alpha');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = get(handles.alphaSlider,'Value');
  case {'overlaymin'}
    % overlayMin = viewGet(view,'overlayMin');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = get(handles.overlayMinSlider,'Value');
  case {'overlaymax'}
    % overlayMax = viewGet(view,'overlayMax');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = get(handles.overlayMaxSlider,'Value');
  case {'rotate'}
    % rotate = viewGet(view,'rotate');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = get(handles.rotateSlider,'Value');
  case {'sliceorientation'}
    % sliceorientation = viewGet(view,'sliceorientation');
    fig = viewGet(view,'fignum');
    handles = guidata(fig);
    val = handles.sliceOrientation;
  case {'cursliceoverlaycoords'}
    % overlayCoords = viewGet(view,'cursliceoverlaycoords');
    val = view.curslice.overlayCoords;
  case {'curslicebasecoords'}
    % baseCoords = viewGet(view,'curslicebasecoords');
    val = view.curslice.baseCoords;

  otherwise
    % viewType dependent
    switch(lower(view.viewType))
      case 'volume'
        val = volumeGet(view,param,varargin{:});
      case 'surface'
        val = surfaceGet(view,param,varargin{:});
      case 'flat'
        val = flatGet(view,param,varargin{:});
      otherwise
        error('Unknown type of View.');
    end
end

return;


% N.B.  For all the routines below, varargin{1} is the varargin passed in by the
% main routine. Hence, the first entry varargin{1}, is the entire varargin
% in the calling routine.
%------------------------------

%------------------------------
function val  = volumeGet(view,param,varargin)

global MLR

val = [];
switch lower(param)

  case {'tseriessize'}
    % size of entire tseries: [ysize xsize zsize nframes]
    %
    % tseriessize = viewGet(view,'tseriessize',scanNum,[groupNum])
    % tseriessize = viewGet(view,'tseriessize',scanNum,[])
    % tseriessize = viewGet(view,'tseriessize',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet tseriessize: must specify scan.');
    end
    scan = varargin{1};
    switch length(varargin)
      case 1
        groupNum = viewGet(view,'currentGroup');
      case 2
        groupNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    datasize = viewGet(view,'dataSize',scan,groupNum);
    nframes = viewGet(view,'nFrames',scan,groupNum);
    val = [datasize nframes];

  case {'datasize','dims'}
    % dims of single temporal frame of functional volume (same as size
    % of parameter map) for a given scan.
    %
    % datasize = viewGet(view,'datasize',scanNum,[groupNum])
    % datasize = viewGet(view,'datasize',scanNum,[])
    % datasize = viewGet(view,'datasize',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet datasize: must specify scan.');
    end
    scan = varargin{1};
    switch length(varargin)
      case 1
        groupNum = viewGet(view,'currentGroup');
      case 2
        groupNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    nscans = viewGet(view,'nscans',groupNum);
    if (nscans >= scan) & (scan > 0)
      val = MLR.groups(groupNum).scanParams(scan).dataSize;
    end

  case {'slicedims'}
    % dims of single slice and single temporal frame of functional
    % volume (same as dims of single slice of parameter map) for a
    % given scan.
    %
    % slicedims = viewGet(view,'slicedims',scanNum,[groupNum])
    % slicedims = viewGet(view,'slicedims',scanNum,[])
    % slicedims = viewGet(view,'slicedims',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet slicedims: must specify scan.');
    end
    scan = varargin{1};
    switch length(varargin)
      case 1
        groupNum = viewGet(view,'currentGroup');
      case 2
        groupNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    datasize = viewGet(view,'datasize',scan,groupNum);
    val = datasize(1:2);

  case {'nslices'}
    % nslices for a given scan
    %
    % nslices = viewGet(view,'nslices',scanNum,[groupNum])
    % nslices = viewGet(view,'nslices',scanNum,[])
    % nslices = viewGet(view,'nslices',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet nslices: must specify scan.');
    end
    scan = varargin{1};
    switch length(varargin)
      case 1
        groupNum = viewGet(view,'currentGroup');
      case 2
        groupNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    datasize = viewGet(view,'datasize',scan,groupNum);
    val = datasize(3);

  otherwise
    if strcmp(mrGetPref('verbose'),'Yes')
      dispViewGetHelp;
    end
    disp(['Invalid parameter for volume view: ',param]);
end
return;

%------------------------------
function val  = surfaceGet(view,param,varargin)

val = [];
switch lower(param)

  otherwise
    error(['Invalid parameter for gray view: ',param]);
end
return;


%------------------------------
function val  = flatGet(view,param,varargin)

val = [];
switch lower(param)

  otherwise
    error(['Invalid parameter for gray view: ',param]);
end
return;

%------------------------------
function dispViewGetHelp()

% open this file and scan the text
fid = fopen(which('viewGet'));
C = textscan(fid,'%s%s','delimiter','{}');
fclose(fid);

collectComments = 0;
commands = {};
for i = 1:length(C{1})
  % collect commens after we have found a case word
  if collectComments & (C{1}{i}(1) == '%')
    if collectComments == 1
      commandComments{length(commands)} = sprintf('   %s',C{1}{i}(3:end));
    else
      commandComments{length(commands)} = sprintf('%s\n   %s',commandComments{length(commands)},C{1}{i}(3:end));
    end
    collectComments = collectComments+1;
  else
    % not a comment, stop collecting comments
    collectComments = 0;
  end
  % get a command
  if (strncmp(C{1}{i},'case',4))
    collectComments = 1;
    commands{end+1} = C{2}{i};
  end
end

% sort everything
[commands sortedIndex] = sort(commands);

% now print out
for i = 1:length(commands)
  disp(commands{i});
  if length(commandComments)>=sortedIndex(i)
    disp(commandComments{sortedIndex(i)});
  end
end

disp(sprintf('\n'));
disp('------------------------- All possible parameters ---------------------');
for i = 1:length(commands)
  fprintf(1,commands{i});
end
fprintf(1,sprintf('\n'));
disp('-----------------------------------------------------------------------');

%------------------------------
function [s g] = getScanAndGroup(view,varg,param)

if ieNotDefined('varg')
  s = viewGet(view,'curScan');
else
  s = varg{1};
end
if length(varg) > 1
  g = varg{2};
else
  g = viewGet(view,'currentGroup');
end
if isempty(g)
  g = viewGet(view,'currentGroup');
end
