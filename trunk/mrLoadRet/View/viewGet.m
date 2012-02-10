function [val val2] = viewGet(view,param,varargin)
%
%   val = viewGet(view,param,varargin)
%
% Read the parameters of a view
% Access to these structures should go through this routine and through
% viewSet.
%
% ====================================================================
% Type viewGet w/out any arguments for the full set of parameters and for
% optional arguments. If you want help on a specific parameter, you
% can also do:
% viewGet([],'parameterName','help')
% ====================================================================
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
% 6/2004 djh
% 11/2006 jlg added help and
% originalScanNum,originalFileName,originalGroupName, scanNum,
% dicomName, dicom, stimFileName, stimFile, concatInfo, transforms, TR
%	$Id$

mrGlobals

% if ieNotDefined('view'), error('No view defined.'); end
if ieNotDefined('param')
  dispViewGetHelp;
  return
elseif ~ieNotDefined('varargin') && (length(varargin)==1) && isstr(varargin{1}) && (strcmp(lower(varargin{1}),'help') || strcmp(lower(varargin{1}),'?'))
  dispViewGetHelp(param);
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
    % view = viewGet([],'view',viewNum)
    if length(varargin) == 1
      if (varargin{1} > 0) & (varargin{1} <= length(MLR.views))
        val = MLR.views{varargin{1}};
      end
    else
      val = view;
    end
  case {'viewtype','type'}
    % viewType = viewGet(view,'viewType')
    val = view.viewType;
  case {'viewnums'}
    % Returns an array containing the viewNum of every active view.
    % Returns empty if there are no open views
    % viewNums = viewGet([],'viewnums')
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
  case {'sessiondescription'}
    % subject = viewGet(view,'sessionDescription')
    val = MLR.session.description;
  case {'operator'}
    % subject = viewGet(view,'operator')
    val = MLR.session.operator;
  case {'magnet'}
    % subject = viewGet(view,'magnet')
    val = MLR.session.magnet;
  case {'coil'}
    % subject = viewGet(view,'coil')
    val = MLR.session.coil;
  case {'protocol'}
    % subject = viewGet(view,'protocol')
    val = MLR.session.protocol;
    
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
  case {'loadedanalyses'}
    % loadedanalyses = viewGet(view,'loadedAnalyses',[groupNum])
    % this stores all the loaded analyses so that we can switch
    % between groups
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      if (groupNum >= 1) & (groupNum <= length(view.loadedAnalyses))
        val = view.loadedAnalyses{groupNum};
      end
    end
  case {'groupscannum'}
    % groupscannum = viewGet(view,'groupscannum',[groupNum])
    % this stores what scan number we were on when we switch
    % between groups
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      if (groupNum >= 1) && (groupNum <= length(view.groupScanNum))
        val = view.groupScanNum(groupNum);
      end
    end
    % default to scan 1
    if isempty(val),val = 1;end
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
  case{'numberofgroups','numgroups','ngroups'}
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
    elseif isstr(varargin{1})
      % if passed in a string, look for the string
      % as a group name
      groupNum = viewGet(view,'groupNum',varargin{1});
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    val = MLR.groups(groupNum).name;
  case{'groupnum'}
    % groupnum = viewGet(view,'groupnum',groupName)
    if isstr(varargin{1})
      groupName = varargin{1};
      groupNames = {MLR.groups(:).name};
      val = find(strcmp(groupName,groupNames));
      % if passed in a valid number just return that number
    elseif isnumeric(varargin{1}) && isequal(size(varargin{1}),[1 1])
      if (varargin{1} >= 1) && (varargin{1} <= viewGet(view,'nGroups'))
        val = varargin{1};
      end
    end
    % scan
  case{'nscans','numberofscans','numscans'}
    % n = viewGet(view,'nScans',[groupNum])
    % n = viewGet([],'nScans',groupNum)
    % n = viewGet(view,'nScans',[])
    % n = viewGet(view,'nScans')
    g = getGroup(view,varargin);
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
    % get the whole auxParams structure. Compare
    % with auxparam which gets one field of auxParams
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).auxParams(s);
    end
  case{'stimfile'}
    % stimfile = viewGet(view,'stimfile',[scanNum],[groupNum],[removeEmpties]);
    % returns a cell array with all the stimfiles associated with the
    % scan (including if it is a concatenation or an average, will show
    % all ths stimfiles from the original scans used to make the scan).
    % If you set removeEmpties (defaults to true) to false, then it
    % will return empties for any original scan that has a missing
    % stimfile
    [s g] = getScanAndGroup(view,varargin,param);
    % get removeEp=mpties
    if length(varargin) >= 3,removeEmpties = varargin{3};else removeEmpties = true;end
    % get the stimFileNames from auxParam, calling viewGetLoadStimFileile
    % to prepend etc directory and load
    val = viewGet(view,'auxParam','stimFileName',s,g,@viewGetLoadStimFile,removeEmpties);
    if ~isempty(val),val = cellArray(val);end
  case{'stimfilename'}
    % stimFileName = viewGet(view,'stimFileName',[scanNum],[groupNum],[removeEmpties]);
    % returns the names of all the stimfiles associated with the scan
    % this will prepend the ETc directory where they should live so
    % gives a fully qualified path. If it is a concatenation or an average,
    % will dispaly all stimfile names from the original scans used to make the scan.
    % If you set removeEmpties (defaults to true) to false, then it
    % will return empties for any original scan that has a missing
    % stimfile
    [s g] = getScanAndGroup(view,varargin,param);
    % get the stimfilenames from auxParam, calling viewGetPrependEtc
    % to prepend etc directory
    val = viewGet(view,'auxParam','stimFileName',s,g,@viewGetPrependEtc,1);
    if ~isempty(val),val = cellArray(val);else val = {};end
  case{'fidinfo'}
    % fidinfo = viewGet(view,'fidInfo',[scanNum],[groupNum],[removeEmpties]);
    % returns a cell array with the fidInfo associated with each scan
    % this is only for varian systems in which there is a fid file
    % linked with the auxParam 'fidFilename'
    % (if it is a concatenation or an average, will show
    % all ths fidInfo from the original scans used to make the scan).
    % If you set removeEmpties (defaults to false) to true, then it
    % will not return empties for any original scan that has a missing
    % fidFilename or missing fid
    [s g] = getScanAndGroup(view,varargin,param);
    % get removeEp=mpties
    if length(varargin) >= 3,removeEmpties = varargin{3};else removeEmpties = false;end
    % get the fidFilename from auxParam, calling viewGetFidInfo
    % to prepend Pre directory and get fidInfo
    val = viewGet(view,'auxParam','fidFilename',s,g,@viewGetFidInfo,removeEmpties);
    if ~isempty(val),val = cellArray(val);end
  case{'fidfilename'}
    % fidFileName = viewGet(view,'fidFileName',[scanNum],[groupNum],[removeEmpties]);
    % returns the names of all the fidfiles associated with the scan
    % this will prepend the Pre directory where they should live so
    % gives a fully qualified path. If it is a concatenation or an average,
    % will dispaly all fidfile names from the original scans used to make the scan.
    % If you set removeEmpties (defaults to false) to true, then it
    % will not return empties for any original scan that has a missing
    % fidFileName
    [s g] = getScanAndGroup(view,varargin,param);
    % get the fidFilenames from auxParam, calling viewGetPrependPre
    % to prepend Pre directory
    val = viewGet(view,'auxParam','fidFilename',s,g,@viewGetPrependPre,1);
    if ~isempty(val),val = cellArray(val);else val = {};end
  case{'auxparam'}
    % n = viewGet(view,'auxParam','paramName',[scanNum],[groupNum],[functionptr],<removeEmpties>)
    % get the named parameter from auxParams. Compare
    % with auxParams which gets the whole structure. Note that
    % this will return a single value if paramName exists for that
    % scan. If it doesn't and exist in the original scan (i.e.
    % for a concat or average), then it will return a cell array
    % of the parameter - one for each of the original scans. <functionptr> is
    % an argument used by other viewGet commands - it calls the function
    % on the auxParam field before returning it.
    if length(varargin) < 1
      disp(sprintf('(viewGet) Need to specify paramName for auxParam'));
      return
    end
    [s g] = getScanAndGroup(view,{varargin{2:end}});
    if isempty(g),g = viewGet(view,'curGroup');end
    nscans = viewGet(view,'nscans',g);
    val = [];
    paramName = fixBadChars(varargin{1});
    % handle auxillary arguments
    if length(varargin) >= 4,fhandle = varargin{4};else fhandle = [];end
    if length(varargin) >= 5,removeEmpties = varargin{5};else removeEmpties = false;end
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),paramName)
        if ~isempty(MLR.groups(g).auxParams(s).(paramName))
          % get the stored stimFileNames
          val = MLR.groups(g).auxParams(s).(paramName);
	  % call function on value (if function handle is specified)
	  if ~isempty(fhandle) val = feval(fhandle,view,val); end
	end
      end
    end
    % if the field does not exist, then check original
    if isempty(val)
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          val{osnum} = viewGet(view,'auxParam',paramName,os(osnum),og(osnum),fhandle,removeEmpties);
	  % don't let the cell arrays get too nested
	  if iscell(val{osnum}) && (length(val{osnum}) == 1)
	    val{osnum} = val{osnum}{1};
	  end
        end
	% remove empties
	if removeEmpties
	  newval = {};
	  for j = 1:length(val)
	    if ~isempty(val{j})
	      newval{end+1} = val{j};
	    end
	  end
	  val = newval;
	end
      end
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
  case{'junkframestotal','totaljunkedframes'}
    % gets all junk frames from this scan and the scans this was made from
    % that is the total of all frames that have been junked. This
    % does not included the number of frames in the junkFrames parameter
    % n = viewGet(view,'totalJunkedFrames',scanNum,[groupNum])
    % n = viewGet([],'totalJunkedFrames',scanNum,groupNum)
    % n = viewGet(view,'totalJunkedFrames',scanNum,[])
    % n = viewGet(view,'totalJunkedFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    % get the frames that _have already been_ junked from this scan
    val = [];
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0) & isfield(MLR.groups(g).scanParams(s),'totalJunkedFrames')
      val = MLR.groups(g).scanParams(s).totalJunkedFrames;
    end
    if isempty(val)
      % figure out how many time series
      numTimeSeries = max(1,length(viewGet(view,'originalFilename',s,g)));
      % add just set the totalJunkedFrames to 0 for each one
      val = zeros(1,numTimeSeries);
    end
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
  case{'nframes','numframes','nvolumes','numvolumes'}
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
    % Frame period is the time from one volume to the next.
    % It is initially extracted from the nifti header. It
    % is *not* necessarily the TR.
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
    ngroups = viewGet(view,'ngroups');
    if (g > 0) & (g <= ngroups)
      nscans = viewGet(view,'nscans',g);
      if (nscans >= s) & (s > 0)
	val = fullfile(viewGet(view,'tseriesDir',g), MLR.groups(g).scanParams(s).fileName);
      end
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
      mrErrorDlg('(viewGet) scanNum: Must specify tseriesFileName.');
    end
    [tseriesFilePath,tseriesFileName] = fileparts(varargin{1});
    if length(varargin) > 1
      groupNums = viewGet(view,'groupNum',varargin{2});
      if isempty(groupNums)
        disp(sprintf('(viewGet) scanNum: Unknown group'));
        return
      end
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
    % check for more than one scan in the same group that
    % has the same filename. This should never happen since
    % filenames are meant to be unique identifiers of the scan
    uniqueGroupNumMatch = unique(groupNumMatch);
    if length(uniqueGroupNumMatch) ~= length(groupNumMatch)
      % init new arrays
      val = [];
      val2 = [];
      % go through all unique group nums
      for i = 1:length(uniqueGroupNumMatch)
	% find which ones are duplicated so we can display
	whichMatch = find(uniqueGroupNumMatch(i) == groupNumMatch);
	if length(whichMatch) > 1
	  % display warning
	  mrWarnDlg(sprintf('(viewGet:scanNum) !!! Scans %s in group %s have the same filename: %s. Returning only scan %i. !!!\nNote this should not happen because tseries names for scans are meant to be unique identifiers. You should probably delete one of the scans.',num2str(scanNumMatch(whichMatch),'%i '),viewGet(view,'groupName',groupNumMatch(whichMatch(1))),viewGet(view,'tSeriesFile',scanNumMatch(whichMatch(1)),groupNumMatch(whichMatch(1))),scanNumMatch(whichMatch(1))));
	end
	% now keep only the first match
	val(end+1) = scanNumMatch(whichMatch(1));
	val2(end+1) = groupNumMatch(whichMatch(1));
      end
    end
    

  case {'scannumfromdescription'}
    % scanNum = viewGet(view,'scanNumFromDescription',description,<groupNum>,<searchType>);
    % get scanNum(s) in current group that have a matching description
    % returns empty if there are no matches.
    % searchType can be any one of:
    %   'anywhere': case insensitive match where the description string can be found anywhere in the
    %      scans description. e.g. a scan with description 'CCW wedges' will match 'ccw'. Be aware
    %      that 'CCW wedges' will also match the string 'CW'.
    %   'anywhereCaseSensitive': like above, but case sensitive.
    %   'exact': case sensitive match where the description string has to exactly match
    %      the scans description. e.g. 'CCW wedges' will only match 'CCW wedges', not 'CCW'
    %   'exactCaseInsensitive': like exact, but not case sensitive.
    %   'regexp': Will use matlab's regular expression function (regexp) to do the match. For
    %      instance, if you want to match 'CW wedges', but not 'CCW wedges', you could
    %      search for '^CW '.
    % The default serachType is 'anywhere'
    if ieNotDefined('varargin')
      disp('(viewGet) scanNumFromDescription: Must specify description.');
      return
    end
    % get the search string
    searchString = varargin{1};
    % get the group num
    if length(varargin)>1
      groupNum = viewGet(view,'groupNum',varargin{2});
    else
      groupNum = viewGet(view,'curGroup');
    end
    % get the serach Type
    if length(varargin)>2
      searchType = varargin{3};
    else
      searchType = 'anywhere';
    end
    if ~any(strcmp(searchType,{'anywhere','anywhereCaseSensitive','exact','exactCaseInsensitive','regexp'}))
      disp(sprintf('(viewGet) scanNumFromDescription: Unknwon searchType %s',searchType));
      return
    end
    % now go do the search
    for scanNum = 1:viewGet(view,'nScans',groupNum)
      thisDescription = viewGet(view,'description',scanNum,groupNum);
      switch searchType
        case 'anywhere'
          if ~isempty(strfind(lower(thisDescription),lower(searchString)));
            val(end+1) = scanNum;
          end
        case 'anywhereCaseSensitive'
          if ~isempty(strfind(thisDescription,searchString))
            val(end+1) = scanNum;
          end
        case 'exact'
          if strcmp(thisDescription,searchString)
            val(end+1) = scanNum;
          end
        case 'exactCaseInsensitive'
          if strcmp(lower(thisDescription),lower(searchString))
            val(end+1) = scanNum;
          end
        case 'regexp'
          if ~isempty(regexp(thisDescription,searchString))
            val(end+1) = scanNum;
          end
      end
    end
  case {'concatinfo'}
    % concatInfo = viewGet(view,'concatInfo',[scanNum],[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    [tseriesPath,tseriesFile] = fileparts(viewGet(view,'tseriesPath',s,g));
    % check for mat file
    matFileName = fullfile(tseriesPath,sprintf('%s.mat',tseriesFile));
    if isfile(matFileName)
      load(matFileName);
      if exist('concatInfo')
        val = concatInfo;
      end
    end
  case {'stimfileold'}
    % stimfile = viewGet(view,'stimfile',scanNum,[groupNum]);
    % this is the old way of getting stimfiles. The new version
    % should be completely compatible with this old version,
    % but leaving this code in here, just in case, for now.
    % came be deleted later. jg 11/2011
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
        elseif isfield(val{j},'stimts')
          val{j}.filetype = 'afni';
        elseif isfield(val{j},'myscreen')
          val{j}.filetype = 'mgl';
        elseif isfield(val{j},'stimvol')
          val{j}.filetype = 'stimvol';
        else
          val{j}.filetype = 'unknown';
        end
      end
    end
  case {'stimfilenameold'}
    % stimFileName = viewGet(view,'stimFileName',scanNum,[groupNum]);
    % stimFileName is returned as a cell array of all stimfiles
    % associated with the scan
    % this is the old way of getting stimfiles. The new version
    % should be completely compatible with this old version,
    % but leaving this code in here, just in case, for now.
    % came be deleted later. jg 11/2011
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    stimFileName{1} = '';
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),'stimFileName')
        if ~isempty(MLR.groups(g).auxParams(s).stimFileName)
          % get the stored stimFileNames
          stimFileName = cellArray(MLR.groups(g).auxParams(s).stimFileName);
          % and prepend the etcDir path on to them
          for i = 1:length(stimFileName)
            stimFileName{i} = fullfile(viewGet(view,'etcDir'),stimFileName{i});
          end
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
  case{'mousedownscancoords'}
    % scanCoords = viewGet(view,'mouseDownScanCoords')
    viewNum = viewGet(view,'viewNum');
    if isfield(MLR,'interrogator')
      if length(MLR.interrogator) >= viewNum
        if isfield(MLR.interrogator{viewNum},'mouseDownScanCoords')
          val = MLR.interrogator{viewNum}.mouseDownScanCoords;
        end
      end
    end
  case{'mousedownbasecoords'}
    % baseCoords = viewGet(view,'mouseDownBaseCoords')
    viewNum = viewGet(view,'viewNum');
    if isfield(MLR,'interrogator')
      if length(MLR.interrogator) >= viewNum
        if isfield(MLR.interrogator{viewNum},'mouseDownBaseCoords')
          val = MLR.interrogator{viewNum}.mouseDownBaseCoords;
        end
      end
    end
  case {'spikeinfo'}
    % eyepos = viewGet(view,'spikeinfo',scanNum,[groupNum]);
    val = [];
    [s g] = getScanAndGroup(view,varargin,param);
    nScans = viewGet(view,'nscans',g);
    nGroups = viewGet(view,'nGroups');
    if (g > 0) &  (length(MLR.groups(g).auxParams) >= s) & (s > 0)
      if isfield(MLR.groups(g).auxParams(s),'spikeinfo')
        val = MLR.groups(g).auxParams(s).spikeinfo;
        % check tSeries name
        if isfield(val,'filename') & ~strcmp(val.filename,viewGet(view,'tSeriesFile',s,g))
          val = [];
        end
      end
    end
  case {'eyepos'}
    % eyepos = viewGet(view,'eyepos',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    [eyeposFileName eyeposNum] = viewGet(view,'eyeposFileName',s,g);
    val = {};
    % cycle over all eyepos files (concatenations and averages,
    if ~isempty(eyeposFileName)
      for j = 1:length(eyeposFileName)
        % load this eyepos file
        if ~isfile(sprintf('%s.mat',stripext(eyeposFileName{j})))
          mrErrorDlg(sprintf('viewGet %s: Could not find eyepos file %s',param,eyeposFileName{j}));
        else
          eyepos=load(eyeposFileName{j});
          f = fieldnames(eyepos);
          % here we get the proper scan.
          val{j} = eyepos.(f{1}).scan{eyeposNum{j}};
          val{j}.hdr = rmfield(eyepos.(f{1}),'scan');
        end
      end
    end
  case {'eyeposfilename'}
    % eyeposFileName = viewGet(view,'eyeposFileName',scanNum,[groupNum]);
    % eyeposFileName is returned as a cell array of all eyepos files
    % associated with the scan
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    eyeposFileName{1} = '';
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),'eyeposFileName')
        if ~isempty(MLR.groups(g).auxParams(s).eyeposFileName)
          eyeposFileName{1} = fullfile(viewGet(view,'etcDir'),MLR.groups(g).auxParams(s).eyeposFileName);
          eyeposNum{1} = MLR.groups(g).auxParams(s).eyeposNum;
        end
      end
    end
    % if the file does not exist, then check original
    if isempty(eyeposFileName{1})
      eyeposFileName = {};
      eyeposNum = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          [newEyeposFileNames newEyeposNum] = viewGet(view,'eyeposFileName',os(osnum),og(osnum));
          for j = 1:length(newEyeposFileNames)
            eyeposFileName{end+1} = newEyeposFileNames{j};
            eyeposNum{end+1} = newEyeposNum{j};
          end
        end
      end
    end
    val = eyeposFileName;
    val2 = eyeposNum;
  case {'3d'}
    % is3d = viewGet(view,'tr',scanNum,[groupNum]);
    % returns whether sequence is 3d or not
    [s g] = getScanAndGroup(view,varargin,param);
    % first get dicom
    dicom = viewGet(view,'dicom',s,g);
    val = [];
    if isempty(dicom)
      return
    end
    for i = 1:length(dicom)
      if isfield(dicom{i},'ACQ') && isfield(dicom{i}.ACQ,'MR_Acquisition_Type_')
        thisval = strcmp(dicom{i}.ACQ.MR_Acquisition_Type_,'3D');
      end
      if (isempty(val))
        val = thisval;
      elseif (val ~= thisval)
        disp(sprintf('(viewGet:3D) Data collected with both 3D and 2D sequence types'))
      end
    end
  case {'tr'}
    % TR is the TR of the pulse sequence. It is *not* necessarily
    % the time between volumes (framePeriod). The TR is extracted
    % from the DICOM header
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
        % if originalFilename is the same as the filename
        % then return empty
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
    % params = viewGet(view,'params',scanNum,[groupNum]);
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
  case{'sformcode','scansformcode'}
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
  case{'scanqformcode'}
    % xform = viewGet(view,'qformCode',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'qform_code')
        val = [];
      else
        val = MLR.groups(g).scanParams(s).niftiHdr.qform_code;
      end
    end
  case{'scanvol2tal'}
    % xform = viewGet(view,'scanVol2tal',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to talairach
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).vol2tal;
    end
  case{'scanvol2mag'}
    % xform = viewGet(view,'scanVol2mag',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to the magnet
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).vol2mag;
    end
  case{'scan2tal'}
    % xform = viewGet(view,'scan2tal',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from this scan to the volume
    % in talairach coordinates. Note that this also has
    % composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      sform_code = viewGet(view,'sformCode',s,g);
      scanSform = viewGet(view,'scanSform',s,g);
      if ~isempty(scanSform)
        if (sform_code == 3)
          val = scanSform * shiftOriginXform;
        elseif (sform_code == 1)
          vol2tal = viewGet(view,'scanVol2tal',s,g);
          vol2mag = viewGet(view,'scanVol2mag',s,g);
          if ~isempty(vol2tal) && ~isempty(vol2mag)
            val = vol2tal * inv(vol2mag) * scanSform * shiftOriginXform;
          end
        end
      end
    end
  case{'scan2mag','scanxform'}
    % xform = viewGet(view,'scan2mag',[scanNum],[groupNum])
    % The scan2mag xform specifies the xform from the scan
    % to the volume in magnet coordinates.
    % If the sform_code is set to 1, then this is the same
    % as scanSform *except* that the origin has been shifted
    % to start at 1,1,1. i.e. It is scanSform * shiftOriginXform
    %
    % Note that before adding the talairach xform code, scan2mag
    % was called scanxform and the code assumed that the sform
    % of the scan contained this xform. If you ask for scanXform
    % you will get the scan2mag xform, but it does *not* have
    % the shiftOriginXform composited so as to be compatible
    % with the old code
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      sform_code = viewGet(view,'sformCode',s,g);
      scanSform = viewGet(view,'scanSform',s,g);
      if ~isempty(scanSform)
        if (sform_code == 1)
          val = scanSform * shiftOriginXform;
        elseif (sform_code == 3)
          vol2tal = viewGet(view,'scanVol2tal',s,g);
          vol2mag = viewGet(view,'scanVol2mag',s,g);
          if ~isempty(vol2tal) && ~isempty(vol2mag)
            val = vol2mag * inv(vol2tal) * scanSform * shiftOriginXform;
          end
        elseif (sform_code == 0)
          % If sform has not been set, then use the transform that
          % transforms this image directly on to the current anatomy
          % using the qform matrices.
          if strcmp(mrGetPref('verbose'),'Yes')
            oneTimeWarning(sprintf('noScanSform_%i_%i',s,g),...
              ['(viewGet:scanXform) sform is not set. Using qform to align '...
              'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
          end
          val = MLR.groups(g).scanParams(s).niftiHdr.qform44 * shiftOriginXform;
        end
      end
      if strcmp(lower(param),'scanxform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    end
  case{'scan2roi'}
    % xform = viewGet(view,'scan2roi',[roiNum],[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from scan coordinates to ROI coordinates
    % It checks whether the ROI and the scan are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param,2);
    nscans = viewGet(view,'nscans',g);
    r = getRoiNum(view,varargin);
    nRois = viewGet(view,'numrois');
    if (nscans >= s) && (s > 0) && (nRois >= r) && (r > 0)
      roi2tal = viewGet(view,'roi2tal',r);
      if ~isempty(roi2tal) % The roi has a Tal xform
        scan2tal = viewGet(view,'scan2tal',s,g); % check scan
        if ~isempty(scan2tal) % -CASE 1-: both the roi and the scan have a Tal xform
          val = inv(roi2tal) * scan2tal; % use it
        else %  -CASE 2-: the roi has a Tal xform but the scan does not
          roi2mag = viewGet(view,'roi2mag',r); % check if the roi has a Mag xform
          if ~isempty(roi2mag) % if it does,
            scan2mag = viewGet(view,'scan2mag',s,g); % check if the scan has a mag xform too
            if ~isempty(scan2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiScanMismatch_%i_%i_%i',r,s,g),...
                ['WARNING: Ignoring the roi talairach transformation, '...
                'because roi has a Talairach transformation, but '...
                'the scan does not. Coordinates are being converted '...
                'through the magnet coordinate frame. If you wanted '...
                'the transformations to use Talairach coordinates '...
                'instead of magnet coordinates, you need to use '...
                'mrAlign to export the talairach transformation to the scan']);
              val = inv(roi2mag) * scan2mag;
            else % if the scan doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noScanXform_%i_%i_%i',r,s,g),...
                ['ERROR: Scan does not have a transform for '...
                'magnet or talairach space. '...
                'Using the identity matrix to transform from '...
                'scan to ROI. Run mrAlign to fix the scan.']);
              val = eye(4);
            end
          else % if the roi does not have a mag xform and the scan does not have a tal xform
            oneTimeWarning(sprintf('incompatibleRoiScan_%i_%i_%i',r,s,g),...
              ['Scan and ROI are not compatible: ROI is in Talairach space '...
              '(and not magnet space) but Scan is not. '...
              'Using the identity matrix to transform from scan to ROI. '...
              'Run mrAlign to get scan and ROI into the same space.']);
            val = eye(4);
          end
        end
      else % The ROI doesn't have a Tal xform...
        roi2mag = viewGet(view,'roi2mag',r);
        if ~isempty(roi2mag) % ... but the ROI does have a mag transform
          scan2mag = viewGet(view,'scan2mag',s,g); % check the scan:
          if ~isempty(scan2mag) % -CASE 3-: both scan and ROI have mag transform
            val = inv(roi2mag) * scan2mag; % use it;
            % but check if scan had a Tal xform so can warn user that it's being ignored
            scan2tal = viewGet(view,'scan2tal',s,g);
            if ~isempty(scan2tal)
              oneTimeWarning(sprintf('roiScanMismatch_%i_%i_%i',r,s,g),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because the scan has a talairach '...
                'transformation but the ROI does not. Coordinates '...
                'are being converted through the magnet '...
                'coordinate frame. If you want convert using the '...
                'talairach transformations, you need '...
                'to export a talairach transformation to the '...
                'ROI by running  mrAlign.']);
            end
          else % -CASE 4-: ROI has a mag xform but scan does not
            % (and ROI doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleRoiScan_%i_%i_%i',r,s,g),...
              ['ROI and Scan are not compatible: ROI is in magnet '...
              'space and Scan is not. Using the '...
              'identity matrix to transform from scan to ROI '...
              'Run mrAlign to get scan and ROI into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a magnet nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i',r,s,g),...
            ['ROI is neither in Magnet space nor in Talairach '...
            'Space. Using the identity matrix '...
            'to transform from base to ROI. Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case{'scan2scan'}
    % xform = viewGet(view,'scan2scan',[fromScanNum],[fromGroupNum],[toScanNum],[toGroupNum])
    % This will return the xform matrix that specifies the
    % transformation from coordinates of the 'From' scan to coordinates of the 'To' scan
    % (E.g., given x,y,s from the 'From' Scan, multiply by the xform calculated in this
    % call to convert to x',y',s' in the 'To' scan.)
    % It checks whether the scans are each in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [sFrom gFrom] = getScanAndGroup(view,varargin,param);
    [sTo gTo] = getScanAndGroup(view,varargin,param,3);
    nscansFrom = viewGet(view,'nscans',gFrom);
    nscansTo = viewGet(view,'nscans',gTo);
    if (nscansFrom >= sFrom) && (sFrom > 0) && (nscansTo >= sTo) && (sTo > 0)
      scan2talFrom = viewGet(view,'scan2tal',sFrom,gFrom); % check if the FROMscan
      % has a tal xform
      if ~isempty(scan2talFrom) % if the FROMscan has a Tal xform
        scan2talTo = viewGet(view,'scan2tal',sTo,gTo); % check TOscan
        if ~isempty(scan2talTo) % -CASE 1-: both scans have a Tal xform, then use it
          % first check if they are the same, then just
          % return identity matrix, % so there's no roundoff
          %  error when compute the inverse
          if isequal(scan2talTo, scan2talFrom)
            val = eye(4);
          else
            val = inv(scan2talTo) * scan2talFrom;
          end
        else %  -CASE 2-: the FROMscan has a Tal xform but the TOscan doesn't
          scan2magFrom = viewGet(view,'scan2mag',sFrom,gFrom); % check if the FROMscan
          % has a Mag xform
          if ~isempty(scan2magFrom) % if it does,
            scan2magTo = viewGet(view,'scan2mag',sTo,gTo); % check if the TOscan
            % has a Mag xform too
            if ~isempty(scan2magTo) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('scanMismatch_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because FROM scan has a Talairach transformation,'...
                ' but TO scan does not. Coordinates are being converted'...
                ' through the magnet coordinate frame.'...
                ' If you wanted the transformations to use Talairach'...
                ' coordinates instead of magnet coordinates, you'...
                ' need to use mrAlign to export the talairach '...
                ' transformation to the TO scan']);
              if isequal(scan2magTo,scan2magFrom) % if they're the same, avoid
                % the roundoff errors of calling inv
                val = eye(4);
              else
                val = inv(scan2magTo) * scan2magFrom;
              end
            else % if the TOscan doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noScanXform_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['ERROR: TO scan does not have a transform for magnet '...
                'or talairach space. Using the identity'...
                ' matrix to transform from scan to scan. '...
                'Run mrAlign to fix the TO scan.']);
              val = eye(4);
            end
          else % if the FROM scan does not have a mag xform
            % and the TO scan does not have a tal xform
            oneTimeWarning(sprintf('incompatibleScans_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
              ['Scans are not compatible: FROM scan is in '...
              'Talairach space but TO scan is not. '...
              'Using the identity matrix to transform from scan to scan. '...
              'Run mrAlign to get scans into the same space.']);
            val = eye(4);
          end
        end
      else % The FROM scan doesn't have a Tal xform...
        scan2magFrom = viewGet(view,'scan2mag',sFrom,gFrom);
        if ~isempty(scan2magFrom) % ... but the FROM scan does have a mag transform
          scan2magTo = viewGet(view,'scan2mag',sTo,gTo); % check the TO scan
          if ~isempty(scan2magTo) % -CASE 3-: both scans have mag transform
            if isequal(scan2magTo,scan2magFrom) % if they're the same, avoid
              % the roundoff errors of calling inv
              val = eye(4);
            else
              val = inv(scan2magTo) * scan2magFrom;
            end
            % but check if TO scan had a Tal xform so can warn user that it's being ignored
            scan2talTo = viewGet(view,'scan2tal',sTo,gTo);
            if ~isempty(scan2talTo)
              oneTimeWarning(sprintf('scanMismatch_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['WARNING: Ignoring the TO scans talairach transformation, '...
                'because the TO scan has a talairach transformation but '...
                'the FROM scan does not. Coordinates are being converted '...
                'through the magnet coordinate frame. '...
                'If you want convert using the talairach transformations, '...
                'you need to export a talairach transformation to '...
                'the TO scan by running  mrAlign.']);
            end
          else % -CASE 4-: FROM scan has a mag xform but TO scan does not
            % (and scan doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleScans_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
              ['Scans are not compatible: FROM scan is in magnet space but '...
              'TO scan is not. Using the identity matrix to transform from '...
              'scan to scan. Run mrAlign to get scans into the same space.']);
            val = eye(4);
          end
        else % error if scan has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
            ['FROM scan is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix to transform from scan to scan.'...
            ' Run mrAlign to fix the scan.']);
          val = eye(4);
        end
      end
    end
  case{'scanqform'}
    % sform = viewGet(view,'scanQform',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % get the field
      if isfield(MLR.groups(g).scanParams(s).niftiHdr,'qform44')
        val = MLR.groups(g).scanParams(s).niftiHdr.qform44;
      end
    end
  case{'scansform'}
    % sform = viewGet(view,'scanSform',[scanNum],[groupNum])
    % The scanSform is the sform set in the nifti header.
    % If the sform_code is set to 1, then this field has
    % been set by mrAlign to be the tansformation from this
    % scan to the volume in magnet coordinates. Note that
    % scanSform does not shift the origin to start at 1,1,1
    % You usually will need to composite this sform with
    % shiftOriginXform to get the scan2mag xform.
    % If the sform_code is set to 3, then this field has
    % been set by mrAlign to be the transformation from
    % this scan to the volume in talairach coordinates.
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
          oneTimeWarning(sprintf('noScanSform_%i_%i',s,g),...
            ['(viewGet:scanXform) sform is not set. Using qform to align '...
            'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
        end
        val = MLR.groups(g).scanParams(s).niftiHdr.qform44;
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
  case{'numberofbasevolumes','numbase'}
    % n = viewGet(view,'numberofBaseVolumes')
    val = length(view.baseVolumes);
  case {'curbase','currentbase'}
    % baseNum = viewGet(view,'currentBase')
    val = view.curBase;
  case{'basenum'}
    % baseNum = viewGet(view,'baseNum',baseName)
    baseName = varargin{1};
    % if numeric there is nothing to do, just return value
    if isnumeric(baseName)
      val = baseName;
    else
      % otherwise look up the baseNum
      baseNames = {view.baseVolumes(:).name};
      val = find(strcmp(baseName,baseNames));
    end
  case{'basevolume','baseanatomy'}
    % basevolume = viewGet(view,'baseVolume',[baseNum])
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
  case {'currentbasename','curbasename'}
    % baseName = viewGet(view,'currentBaseName')
    baseNames = viewGet(view,'baseNames');
    curBase = viewGet(view,'curBase');
    if isempty(curBase)
      val = 'NoBase';
    else
      val = baseNames{curBase};
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
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).name;
    end
  case {'basedata'}
    % basedata = viewGet(view,'basedata',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).data;
    end
  case {'basehdr'}
    % basedata = viewGet(view,'basehdr',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).hdr;
    end
  case {'base'}
    % basedata = viewGet(view,'base',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b);
    end
  case {'basecoordmappath'}
    % basedata = viewGet(view,'baseCoordMapPath',[baseNum],[corticalDepth])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    val = [];
    if b & (b > 0) & (b <= n)
      if isfield(view.baseVolumes(b),'coordMap')
        if ~isempty(view.baseVolumes(b).coordMap)
          val = view.baseVolumes(b).coordMap.path;
	  % check if the path exists
	  if ~isdir(val)
	    % guess where the path is in from the canonical volume directory
	    volumeDirectory = mrGetPref('volumeDirectory');
	    volumeDirectoryList = dir(volumeDirectory);
	    innerCoordsFilename = view.baseVolumes(b).coordMap.innerCoordsFileName;
	    subjectDir = '';
	    % tell user what we are doing
	    disp(sprintf('(viewGet:baseCoordMapPath) Surface directory %s for base %s does not exist, searching in volumeDirectory: %s',val,viewGet(view,'baseName',b),volumeDirectory));
	    for i = 1:length(volumeDirectoryList)
	      % for each volume directory in the list, see if the directory name
	      % matches the first part of the baseVolumes anatomy (this assumes
	      % that people use a convention like calling the directory s001 and
	      % calling the anatomy file s001anatomy or something like that.
	      matchName = strfind(view.baseVolumes(b).coordMap.anatFileName,volumeDirectoryList(i).name);
	      if ~isempty(matchName) && isequal(matchName(1),1)
		% we have a match, for the subject directory under the volume direcotry
		subjectDir = fullfile(volumeDirectory,volumeDirectoryList(i).name);
		break;
	      end
	    end
	    % not found, give up
	    if isempty(subjectDir)
	      disp(sprintf('(viewGet:baseCoordMapPath) Could not find a matchind subjectDir in %s',volumeDirectory));
	    else
	      % set to return subjectDir in case we don't find the surfRelax directory
	      val = subjectDir;
	      % Check for a directory right under this one called any of the following
	      possibleSurfDirNames = {'surfRelax','surfaces','surface','Surface','Surfaces'};
	      for i = 1:length(possibleSurfDirNames)
		if isdir(fullfile(subjectDir,possibleSurfDirNames{i}))
		  % then return that, if it contains the WM surface file
		  val = fullfile(subjectDir,possibleSurfDirNames{i});
		  if isfile(fullfile(val,innerCoordsFilename)),return,end
		end
	      end
	      % look for the FreeSurfer directory, which should either directly contain the innerCoordsFilename
	      % or have a sub directory called surfRelax which does
	      subjectDirListing = dir(subjectDir);
	      for i = 1:length(subjectDirListing)
		if subjectDirListing(i).isdir
		  % see if the file is directly here
		  if isfile(fullfile(subjectDir,subjectDirListing(i).name,innerCoordsFilename))
		    % found it, return the directory
		    val = fullfile(subjectDir,subjectDirListing(i).name);
		    return
		  end
		  % look for surfRelaxDir
		  surfRelaxDir = fullfile(subjectDir,subjectDirListing(i).name,'surfRelax');
		  if isdir(surfRelaxDir)
		    % found surfRelaxDir
		    val = surfRelaxDir;
		    % return this one if we find the innerCoords - it is likely correct, if not
		    % will continue searching to see if we find the file in some other directory
		    if isfile(fullfile(val,innerCoordsFilename))
		      return
		    end
		  end
		end
	      end
	    end
	  end
        end
      end
    end
  case {'basecoordmap'}
    % basedata = viewGet(view,'baseCoordMap',[baseNum],[corticalDepth])
    b = getBaseNum(view,varargin);
    % get cortical depth
    if ieNotDefined('varargin') || (length(varargin)<2)
      corticalDepth = viewGet(view,'corticalDepth');
    else
      corticalDepth = varargin{2};
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).coordMap;
      % see if the coordMap is calculated for the correct cortical depth
      if ~isempty(val) && (~isfield(val,'corticalDepth') || (val.corticalDepth ~= corticalDepth))
        if isfield(val,'innerCoords') && isfield(val,'outerCoords')
          % if not, then we have to do it
          %	  val.coords = (1-corticalDepth)*val.innerCoords + corticalDepth*val.outerCoords;
          val.coords = val.innerCoords + corticalDepth*(val.outerCoords-val.innerCoords);
          val.corticalDepth = corticalDepth;
        end
      end
    end
  case {'basesurface'}
    % basedata = viewGet(view,'baseSurface',[baseNum],[corticalDepth])
    b = getBaseNum(view,varargin);
    % get cortical depth
    if ieNotDefined('varargin') || (length(varargin)<2)
      corticalDepth = viewGet(view,'corticalDepth');
    else
      corticalDepth = varargin{2};
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      % get coordmap
      coordMap = view.baseVolumes(b).coordMap;
      if isfield(coordMap,'innerVtcs') && isfield(coordMap,'outerVtcs')
        % find intermideate values
        vtcs = coordMap.innerVtcs + corticalDepth*(coordMap.outerVtcs-coordMap.innerVtcs);
        val.tris = coordMap.tris;
      elseif isfield(coordMap,'innerVtcs')
        % if only inner is present then just return that
        vtcs = coordMap.innerVtcs;
        val.tris = coordMap.tris;
      else
        vtcs = [];
        val.tris = [];
      end
      % center surface
      if ~isempty(vtcs)
        val.vtcs(:,1) = vtcs(:,2);
        val.vtcs(:,2) = vtcs(:,1);
        val.vtcs(:,3) = vtcs(:,3);
      end
    end
  case {'baseclip'}
    % baseclip = viewGet(view,'baseclip',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).clip;
    end
  case {'baserotate'}
    % baserotate = viewGet(view,'baserotate',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      baseType = viewGet(view,'baseType',b);
      if baseType <= 1
        val = view.baseVolumes(b).rotate;
      else
        if isfield(view.baseVolumes(b).coordMap,'rotate')
          val = view.baseVolumes(b).coordMap.rotate;
        else
          val = 0;
        end
      end
    end
  case {'basetilt'}
    % baserotate = viewGet(view,'baseTilt',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).tilt;
    end
  case {'basecurslice','baseslice'}
    % baseslice = viewGet(view,'baseslice',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).curSlice;
    end
  case {'basesliceorientation'}
    % basesliceOrientation = viewGet(view,'baseSliceOrientation',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).sliceOrientation;
    end
  case {'baserange'}
    % baserange = viewGet(view,'baserange',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).range;
    end
  case {'basegamma'}
    % baseGamma = viewGet(view,'baseGamma',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).gamma;
    end
  case {'basetype'}
    % baserange = viewGet(view,'baserange',[baseNum])
    % 0 is for a regular volume, 1 is for a flat, 2 is for a surface
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).type;
    end
  case {'basesurfacedims'}
    % basedims = viewGet(view,'basedims',[baseNum])
    % basedims = viewGet(view,'basedims',[])
    % basedims = viewGet(view,'basedims')
    b = getBaseNum(view,varargin);
    baseType = viewGet(view,'baseType',b);
    if baseType == 2
      baseDims = view.baseVolumes(b).coordMap.dims;
      val(1) = baseDims(2);
      val(2) = baseDims(1);
      val(3) = baseDims(3);
    end
  case {'basedims'}
    % basedims = viewGet(view,'basedims',[baseNum])
    % basedims = viewGet(view,'basedims',[])
    % basedims = viewGet(view,'basedims')
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = size(baseVolume.data);
      % for an image with only one slice
      if length(val) == 2
        val(3) = 1;
      end
    end
  case{'base2tal'}
    % xform = viewGet(view,'base2tal',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from this base to the volume
    % in talairach coordinates. Note that this also has
    % composited the shiftOriginXform.
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numBase');
    if (b > 0) && (b <= n)
      sform_code = viewGet(view,'baseSformCode',b);
      baseSform = viewGet(view,'baseSform',b);
      if ~isempty(baseSform)
        if (sform_code == 3)
          val = baseSform * shiftOriginXform;
        elseif (sform_code == 1)
          vol2mag = viewGet(view,'baseVol2mag',b);
          vol2tal = viewGet(view,'baseVol2tal',b);
          if ~isempty(vol2mag) && ~isempty(vol2tal)
            val = vol2tal * inv(vol2mag) * baseSform * shiftOriginXform;
          end
        end
      end
    end
  case{'base2mag','basexform'}
    % xform = viewGet(view,'base2mag',[baseNum])
    % The base2mag xform specifies the xform from the base
    % to the volume in magnet coordinates.
    % If the sform_code is set to 1, then this is the same
    % as baseSform *except* that the origin has been shifted
    % to start at 1,1,1. i.e. It is baseSform * shiftOriginXform
    %
    % Note that previous to adding the talairach transformation
    % code, base2mag used to be called baseXform which used
    % to be assumed to be the same as baseSform (note that
    % if you ask for basexform then shiftOriginXform is *not*
    % composited with the transformation--since this is the
    % way the old code worked).
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numBase');
    if (b > 0) && (b <= n)
      sform_code = viewGet(view,'baseSformCode',b);
      baseSform = viewGet(view,'baseSform',b);
      if ~isempty(baseSform)
        if (sform_code == 1)
          val = baseSform * shiftOriginXform;
        elseif (sform_code == 3)
          vol2mag = viewGet(view,'baseVol2mag',b);
          vol2tal = viewGet(view,'baseVol2tal',b);
          if ~isempty(vol2mag) && ~isempty(vol2tal)
            val = vol2mag * inv(vol2tal) * baseSform * shiftOriginXform;
          end
        elseif (sform_code == 0)
          % If sform has not been set, then use the transform that
          % transforms this image directly on to the current anatomy
          % using the qform matrices.
          if strcmp(mrGetPref('verbose'),'Yes')
            oneTimeWarning(sprintf('noBaseSform_%i',b),...
              ['(viewGet:baseXform) sform is not set. Using qform to align '...
              'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
          end
          baseqform = viewGet(view,'baseqform');
          val = baseqform * shiftOriginXform;
        end
      end
      % if we are being asked for baseXform
      if strcmp(lower(param),'basexform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    end
  case{'base2base'}
    % xform = viewGet(view,'base2base',[baseNum1],[baseNum2])
    % This will return the xform matrix that specifies the
    % transformation from the current base coordinates to the specified
    % base (baseNum1)'s coordinates. If you specify baseNum2 it will
    % calculate the xform matrix from baseNum2 to baseNum1
    if length(varargin) < 1, disp(sprintf('(viewGet:base2base) Most specify baseNum1'));return,end
    % get the from base
    if length(varargin) == 1
      baseFromNum = viewGet(view,'curBase');
    else
      baseFromNum = varargin{2};
    end
    % make sure we have base numbers (and not names)
    baseToNum = viewGet(view,'baseNum',varargin{1});
    baseFromNum = viewGet(view,'baseNum',baseFromNum);
    % go through the magnet coordinates
    base2magTo = viewGet(view,'base2mag',baseToNum);
    base2magFrom = viewGet(view,'base2mag',baseFromNum);
    % if neither is empty, then return it
    if ~isempty(base2magTo) && ~isempty(base2magFrom)
      val = inv(base2magTo)*base2magFrom;
    else
      disp(sprintf('(viewGet:base2base) Could not compute transform of base %s to base %s',viewGet(view,'baseName',baseToNum),viewGet(view,'baseName',baseFromNum)));
    end
  case{'base2scan'}
    % xform = viewGet(view,'base2scan',[scanNum],[groupNum],[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from base coordinates to scan coordinates
    % It checks whether the scan and the base are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    b = getBaseNum(view,varargin,3);
    n = viewGet(view,'numBase');
    if (b > 0) & (b <= n) & (nscans >= s) & (s > 0)
      scan2tal = viewGet(view,'scan2tal',s,g);
      if ~isempty(scan2tal) % The scan has a Tal xform
        base2tal = viewGet(view,'base2tal',b); % check base
        if ~isempty(base2tal) % -CASE 1-: both the scan and the base have a Tal xform
          val = inv(scan2tal) * base2tal; % use it
        else %  -CASE 2-: the scan has a Tal xform but the base does not
          scan2mag = viewGet(view,'scan2mag',s,g); % check if the scan has a Mag xform
          if ~isempty(scan2mag) % if it does,
            base2mag = viewGet(view,'base2mag',b); % check if the base has a mag xform too
            if ~isempty(base2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('scanBaseMismatch_%i_%i_%i',s,g,b),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because scan '...
                'has a Talairach transformation, but the base does not. '...
                'Coordinates '...
                'are being converted through the magnet coordinate frame. '...
                'If you wanted '...
                'the transformations to use Talairach coordinates instead '...
                'of magnet '...
                'coordinates, you need to use mrAlign to export the talairach '...
                'transformation to the base']);
              val = inv(scan2mag) * base2mag;
            else % if the base doesn't have either xform, that's an error. Give a warning
              % and use the identity matrix
              oneTimeWarning(sprintf('noBaseXform_%i_%i_%i',s,g,b),...
                ['ERROR: Base does not have a transform for magnet or talairach '...
                'space. '...
                'Using the identity matrix to transform from base to scan. Run '...
                'mrAlign to fix the base.']);
              val = eye(4);
            end
          else % if the scan does not have a mag xform and the base does not have a tal xform
            oneTimeWarning(sprintf('incompatibleB1S3_%i_%i_%i',s,g,b),...
              ['Base and Scan are not compatible: Scan is in Talairach space '...
              '(and not magnet space) but Base is not. '...
              'Using the identity matrix to transform from base to scan. Run '...
              'mrAlign to get base and scan into the same space.']);
            val = eye(4);
          end
        end
      else % The scan doesn't have a Tal xform...
        scan2mag = viewGet(view,'scan2mag',s,g);
        if ~isempty(scan2mag) % ... but the scan does have a mag transform
          base2mag = viewGet(view,'base2mag',b); % check the base:
          if ~isempty(base2mag) % -CASE 3-: both base and scan have mag transform
            val = inv(scan2mag) * base2mag; % use it
            base2tal = viewGet(view,'base2tal',b); % but check if base had a Tal xform
            % so can warn user that it's being ignored
            if ~isempty(base2tal)
              oneTimeWarning(sprintf('scanBaseMismatch_%i_%i_%i',s,g,b),...
                ['WARNING: Ignoring the base talairach transformation, because '...
                'the base has a talairach '...
                'transformation but the scan does not. Coordinates are being '...
                'converted through the magnet '...
                'coordinate frame. If you want convert using the talairach '...
                'transformations, you need '...
                'to export a talairach transformation to the scan by running '...
                ' mrAlign.']);
            end
          else % -CASE 4-: Scan has a mag xform but base does not (and scan doesn't have
            % a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleB3S1_%i_%i_%i',s,g,b),...
              ['Base and Scan are not compatible: Scan is in magnet space '...
              'and Base is not. Using the '...
              'identity matrix to transform from base to scan. Run mrAlign '...
              'to get base and scan into the same space.']);
            val = eye(4);
          end
        else % error if scan has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i',s,g,b),...
            ['Scan is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix '...
            'to transform from base to scan. Run mrAlign to fix the scan.']);
          val = eye(4);
        end
      end
    end
  case{'base2roi'}
    % xform = viewGet(view,'base2roi',[roiNum],[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from base coordinates to ROI coordinates
    % It checks whether the ROI and the base are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    b = getBaseNum(view,varargin,2);
    n = viewGet(view,'numBase');
    r = getRoiNum(view,varargin);
    nRois = viewGet(view,'numrois');
    if (b > 0) & (b <= n) & (nRois >= r) & (r > 0)
      roi2tal = viewGet(view,'roi2tal',r);
      if ~isempty(roi2tal) % The roi has a Tal xform
        base2tal = viewGet(view,'base2tal',b); % check base
        if ~isempty(base2tal) % -CASE 1-: both the roi and the base have a Tal xform
          val = inv(roi2tal) * base2tal; % use it
        else %  -CASE 2-: the roi has a Tal xform but the base does not
          roi2mag = viewGet(view,'roi2mag',r); % check if the roi has a Mag xform
          if ~isempty(roi2mag) % if it does,
            base2mag = viewGet(view,'base2mag',b); % check if the base has a mag xform too
            if ~isempty(base2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiBaseMismatch_%i_%i',r,b),...
                ['WARNING: Ignoring the roi talairach transformation, '...
                'because roi has a Talairach transformation, but '...
                'the base does not. Coordinates are being converted '...
                'through the magnet coordinate frame. If you wanted '...
                'the transformations to use Talairach coordinates '...
                'instead of magnet coordinates, you need to use '...
                'mrAlign to export the talairach transformation to the base']);
              val = inv(roi2mag) * base2mag;
            else % if the base doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noBaseXform_%i_%i',r,b),...
                ['ERROR: Base does not have a transform for '...
                'magnet or talairach space. '...
                'Using the identity matrix to transform from '...
                'base to ROI. Run mrAlign to fix the base.']);
              val = eye(4);
            end
          else % if the roi does not have a mag xform and the base does not have a tal xform
            oneTimeWarning(sprintf('incompatibleRoiBase_%i_%i',r,b),...
              ['Base and ROI are not compatible: ROI is in Talairach space '...
              '(and not magnet space) but Base is not. '...
              'Using the identity matrix to transform from base to ROI. '...
              'Run mrAlign to get base and ROI into the same space.']);
            val = eye(4);
          end
        end
      else % The ROI doesn't have a Tal xform...
        roi2mag = viewGet(view,'roi2mag',r);
        if ~isempty(roi2mag) % ... but the ROI does have a mag transform
          base2mag = viewGet(view,'base2mag',b); % check the base:
          if ~isempty(base2mag) % -CASE 3-: both base and ROI have mag transform
            val = inv(roi2mag) * base2mag; % use it
            base2tal = viewGet(view,'base2tal',b); % but check if base had a Tal xform
            % so can warn user that it's being ignored
            if ~isempty(base2tal)
              oneTimeWarning(sprintf('roiBaseMismatch_%i_%i',r,b),...
                ['WARNING: Ignoring the base talairach transformation, '...
                'because the base has a talairach '...
                'transformation but the ROI does not. Coordinates '...
                'are being converted through the magnet '...
                'coordinate frame. If you want convert using the '...
                'talairach transformations, you need '...
                'to export a talairach transformation to the '...
                'ROI by running  mrAlign.']);
            end
          else % -CASE 4-: ROI has a mag xform but base does not (and ROI doesn't have
            % a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleRoiBase_%i_%i',r,b),...
              ['Base and ROI are not compatible: ROI is in magnet '...
              'space and Base is not. Using the '...
              'identity matrix to transform from base to ROI. '...
              'Run mrAlign to get base and ROI into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a magnet nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i',r,b),...
            ['ROI is neither in Magnet space nor in Talairach '...
            'Space. Using the identity matrix '...
            'to transform from base to ROI. Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case {'basevol2tal'}
    % basexform = viewGet(view,'baseVol2tal',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to talairach
    % coordinates of the base volume that this base was aligned to.
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.vol2tal;
    end
  case {'basevol2mag'}
    % basexform = viewGet(view,'baseVol2mag',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to magnet
    % coordinates of the base volume that this base was aligned to.
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.vol2mag;
    end
  case {'basesform'}
    % baseSform = viewGet(view,'baseSform',[baseNum])
    % This returns the sform in the nifit header
    % which if the sform_code is set to 1 by mrAlign
    % should be the xformation of the base to the
    % volume in magnet coordinates. Note that if
    % this base is the "canonical" volume then
    % this is the same as the qform.
    % if sform_code is set to 3 then this is the
    % trasnformation of the base to the volume in
    % talairach coordinates
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.sform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'baseqform'}
    % basexform = viewGet(view,'baseqform',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.qform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'basesformcode'}
    % baseSformCode = viewGet(view,'baseSformCode',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.sform_code;
    else
      % placeholder when base images not loaded
      val = [];
    end
  case {'basevolpermutation'}
    % basevolpermutation = viewGet(view,'basevolpermutation',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.permutationMatrix;
    else
      % placeholder when base images not loaded
      val = eye(3);
    end
  case {'basesliceindex'}
    % basesliceindex = viewGet(view,'basesliceindex',[baseNum])
    % This admittedly arcane logic is also used in mrAlignGUI. If you
    % change this code, please make corresponding changes in that
    % function. Permutation matrix is set by loadAnat using even more
    % arcane logic.
    b = getBaseNum(view,varargin);
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
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.pixdim([2,3,4])';;
    end
    
    % ROI
  case {'showrois'}
    % show = viewGet(view,'showROIs')
    val = view.showROIs;
  case {'labelrois'}
    % show = viewGet(view,'showROIs')
    % returns a boolean that is set to 1 if
    % ROI labels are turned on, or 0 if not.
    val = view.labelROIs;
  case{'numberofrois','numrois','nrois'}
    % n = viewGet(view,'numberofROIs')
    val = length(view.ROIs);
  case{'currentroi','currentroinum','curroi'}
    % roiNum = viewGet(view,'currentROI')
    val = view.curROI;
  case{'roigroup','currentroigroup'}
    % roiNum = viewGet(view,'roigroup')
    % gets the roiNums of the current roi group.
    roiGroupNames = view.roiGroup;
    val = [];
    for i = 1:viewGet(view,'nrois');
      roiName = viewGet(view,'roiName',i);
      if any(strcmp(roiName,roiGroupNames))
	val(end+1) = i;
      end
    end
  case{'roigroupnames','currentroigroupnames'}
    % roiNum = viewGet(view,'roigroup')
    % gets the roiNames of the current roi group.
    val = view.roiGroup;
  case{'roinum'}
    % roiNum = viewGet(view,'roiNum',roiName)
    % This will return the roiNum that correspondes
    % to the passed in name, if that roi is loaded into
    % the view. If roiName is a number, it will return
    % that number if roiNum is a valid roiNum.
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet roiNum: must specify ROI name');
    end
    ROIname = varargin{1};
    % if no ROIs are loaded than return empty
    nROIs = viewGet(view,'nROIs');
    if nROIs == 0
      return
    end
    % check for number
    if isscalar(ROIname)
      if (ROIname >= 1) && (ROIname <= viewGet(view,'nROIs'))
        val = round(ROIname);
      end
      return
    end
    ROInames = {view.ROIs(:).name};
    % note that it is possible to have more than one
    % ROI with the same name, so pick the last one in
    % the list, since we usually are getting the number
    % for the last ROI that was just created (This should not
    % happen anymore since viewSet checks for duplicates).
    val = find(strcmp(ROIname,ROInames));
    if length(val) > 1
      mrWarnDlg(sprintf('(viewGet) you have multiple ROIs with the name %s',ROIname));
      val = last(val);
    end
  case{'roi'}
    % roi = viewGet(view,'roi',[roiNum])
    r = getRoiNum(view,varargin);
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
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).name;
    end
  case{'roicoords'}
    % roiCoords = viewGet(view,'roiCoords',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).coords;
    end
  case{'prevroicoords'}
    % prevRoiCoords = viewGet(view,'prevRoiCoords')
    val = view.prevROIcoords;
  case{'roicolor'}
    % roicolor = viewGet(view,'roicolor',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).color;
    end
  case{'roicolorrgb'}
    % roicolor = viewGet(view,'roicolorRGB',[roiNum])
    % roicolor = viewGet(view,'roicolorRGB',[])
    % roicolor = viewGet(view,'roicolorRGB')
    if ieNotDefined('varargin')
      roicolor = viewGet(view,'roiColor');
    else
      roicolor = viewGet(view,'roiColor',varargin{1});
    end
    val = color2RGB(roicolor);
  case{'roivol2mag'}
    % roiVol2mag = viewGet(view,'roiVol2mag',[roiNum])
    % returns the xform of the volume coordinates to
    % the magnet coordinates for the canonical volume
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).vol2mag;
    end
  case{'roivol2tal'}
    % roiVol2tal = viewGet(view,'roiVol2tal',[roiNum])
    % returns the xform of the volume coordinates to
    % the talairach coordinates for the canonical volume
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).vol2tal;
    end
  case{'roisformcode'}
    % roiSformCode = viewGet(view,'roiSformCode',[roiNum])
    % returns the sFormCode of the transform saved in
    % view.ROIs(r).sformcode, which is the sformcode of the
    % base on which the ROI was originally defined.
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).sformCode;
    end
  case{'roi2mag','roixform'}
    % roi2mag = viewGet(view,'roi2mag',[roiNum])
    % roi2mag returns the transformation of the roi
    % to the volume in magnet coordinates. It has
    % composited on it the shiftOriginXform.
    % The only difference between this and roiXform
    % is the incluson of the shiftOriginXform
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      sform_code = viewGet(view,'roiSformCode',r);
      if (sform_code == 1)
        val = view.ROIs(r).xform * shiftOriginXform;
        % if the sform_code was set to 0, then either
        % it has a xform that came from the bases qform
        % or none at all (in which case return identity)
      elseif (sform_code == 0)
        if isempty(view.ROIs(r).xform)
          val = eye(4) * shiftOriginXform;
        else
          val = view.ROIs(r).xform * shiftOriginXform;
        end
      elseif (sform_code == 3)
        vol2tal = viewGet(view,'roiVol2tal',r);
        vol2mag = viewGet(view,'roiVol2mag',r);
        if ~isempty(vol2tal) && ~isempty(vol2mag)
          val = vol2mag * inv(vol2tal) * view.ROIs(r).xform * shiftOriginXform;
        end
      end
      if strcmp(lower(param),'roixform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    else
      val = eye(4);
    end
  case{'roi2tal'}
    % roi2tal = viewGet(view,'roi2tal',[roiNum])
    % roi2tal returns the transformation of the roi
    % to the volume in talairach coordinates. It has
    % composited on it the shiftOriginXform.
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      sform_code = viewGet(view,'roisformCode',r);
      if (sform_code == 3)
        val = view.ROIs(r).xform * shiftOriginXform;
      elseif (sform_code == 1)
        vol2tal = viewGet(view,'roiVol2tal',r);
        vol2mag = viewGet(view,'roiVol2mag',r);
        if ~isempty(vol2tal) && ~isempty(vol2mag)
          val = vol2tal * inv(vol2mag) * view.ROIs(r).xform * shiftOriginXform;
        end
      end
    else
      val = eye(4);
    end
  case{'roi2roi'}
    % xform = viewGet(view,'roi2roi',[fromROINum],[toROINum])
    % This will return the xform matrix that specifies the
    % transformation from coordinates of the 'From' ROI to coordinates of the 'To' ROI
    % (E.g., given x,y,s from the 'From' ROI, multiply by the xform calculated in this
    % call to convert to x',y',s' in the 'To' ROI.)
    % It checks whether the ROIs are each in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    rFrom = getRoiNum(view,varargin);
    rTo = getRoiNum(view,varargin,2);
    n = viewGet(view,'numberofROIs');
    if (n >= rFrom) && (rFrom > 0) && (n >= rTo) && (rTo > 0)
      roi2talFrom = viewGet(view,'roi2tal',rFrom); % check if FROMroi has talXform
      if ~isempty(roi2talFrom) % if the FROMroi has a Tal xform
        roi2talTo = viewGet(view,'roi2tal',rTo); % check TOroi
        if ~isempty(roi2talTo) % -CASE 1-: both rois have a Tal xform, then use it
          % first check if they are the same, then just
          % return identity matrix, so there's no roundoff
          % error when compute the inverse
          if isequal(roi2talTo, roi2talFrom)
            val = eye(4);
          else
            val = inv(roi2talTo) * roi2talFrom;
          end
        else %  -CASE 2-: the FROM roi has a Tal xform but the TO roi doesn't
          roi2magFrom = viewGet(view,'roi2mag',rFrom); % check if the FROMroi has magXform
          if ~isempty(roi2magFrom) % if it does,
            roi2magTo = viewGet(view,'roi2mag',rTo); % check if the TOroi has magXform
            if ~isempty(roi2magTo) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiMismatch_%i_%i',rFrom,rTo),...
                ['WARNING: Ignoring the ROI talairach transformation, '...
                'because FROM ROI has a Talairach transformation,'...
                ' but TO ROI does not. Coordinates are being converted'...
                ' through the magnet coordinate frame.'...
                ' If you wanted the transformations to use Talairach'...
                ' coordinates instead of magnet coordinates, you'...
                ' need to use mrAlign to export the talairach '...
                ' transformation to the TO ROI']);
              if isequal(roi2magTo,roi2magFrom) % if they're the same, avoid
                % the roundoff errors of calling inv
                val = eye(4);
              else
                val = inv(roi2magTo) * roi2magFrom;
              end
              
            else % if the TOroi doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noRoiXform_%i_%i',rFrom,rTo),...
                ['ERROR: TO ROI does not have a transform for magnet '...
                'or talairach space. Using the identity '...
                'matrix to transform from ROI to ROI. '...
                'Run mrAlign to fix the TO ROI.']);
              val = eye(4);
            end
          else % if the FROM ROI does not have a mag xform
            % and the TO ROI does not have a tal xform
            oneTimeWarning(sprintf('incompatibleROIs_%i_%i',rFrom,rTo),...
              ['ROIs are not compatible: FROM ROI is in '...
              'Talairach space but TO ROI is not. '...
              'Using the identity matrix to transform from ROI to ROI. '...
              'Run mrAlign to get ROIs into the same space.']);
            val = eye(4);
          end
        end
      else % The FROM ROI doesn't have a Tal xform...
        roi2magFrom = viewGet(view,'roi2mag',rFrom);
        if ~isempty(roi2magFrom) % ... but the FROM ROI does have a mag transform
          roi2magTo = viewGet(view,'roi2mag',rTo); % check the TO ROI
          if ~isempty(roi2magTo) % -CASE 3-: both ROIs have mag transform
            if isequal(roi2magTo,roi2magFrom) % if they're the same, avoid
              % the roundoff errors of calling inv
              val = eye(4);
            else
              val = inv(roi2magTo) * roi2magFrom;
            end
            % but check if TO scan had a Tal xform so can warn user that it's being ignored
            roi2talTo = viewGet(view,'roi2tal',rTo);
            if ~isempty(roi2talTo)
              oneTimeWarning(sprintf('roiMismatch_%i_%i',rFrom,rTo),...
                ['WARNING: Ignoring the TO ROIs talairach transformation, '...
                'because the TO ROI has a talairach transformation but '...
                'the FROM ROI does not. Coordinates are being converted '...
                'through the magnet coordinate frame. '...
                'If you want convert using the talairach transformations, '...
                'you need to export a talairach transformation to '...
                'the TO ROI by running  mrAlign.']);
            end
          else % -CASE 4-: FROM ROI has a mag xform but TO ROI does not
            % (and FROM ROI doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleROIs_%i_%i',rFrom,rTo),...
              ['ROIs are not compatible: FROM ROI is in magnet space but '...
              'TO ROI is not. Using the identity matrix to transform from '...
              'ROI to ROI. Run mrAlign to get ROIs into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i',rFrom,rTo),...
            ['FROM ROI is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix to transform from ROI to ROI.'...
            ' Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case{'roinotes'}
    % roinotesm = viewGet(view,'roinotes',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).notes;
    else
      val = '';
    end
  case{'roivoxelsize'}
    % roivoxelsize = viewGet(view,'roivoxelsize',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).voxelSize;
    end
  case('roivolume')
    % roivolume = viewGet(view,'roivolume',[roiNum])
    r = getRoiNum(view,varargin);
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
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).date;
    end
    
    % Cache
  case{'roicacheid'}
    % cacheID = viewGet(view,'ROICacheID')
    % cacheID = viewGet(view,'ROICacheID',roinum)
    if ieNotDefined('varargin')
      % if we are not passed in a specific ROI, then
      % we are doing all ROIs
      roiNums = 1:length(view.ROIs);
    else
      roiNums = varargin{1};
    end
    rotate = viewGet(view,'rotate');
    % go through all the bases and look for
    % the first one with a matching xform and
    % voxel size. We do this because flat maps
    % all have the same baseXform and voxle size
    % and so there is no need to recompute them
    % for each base map
    baseXform = viewGet(view,'baseXform');
    baseVoxelSize = viewGet(view,'baseVoxelSize');
    baseMatch = viewGet(view,'curbase');
    for bnum = 1:viewGet(view,'numBase')
      if (isequal(baseXform,viewGet(view,'baseXform',bnum)) && ...
          isequal(baseVoxelSize,viewGet(view,'baseVoxelSize',bnum)))
        baseMatch = bnum;
        break;
      end
    end
    baseName = viewGet(view,'baseName',baseMatch);
    val = sprintf('%s_%i',baseName,rotate);
    for i = roiNums
      val = sprintf('%s_%s_%i',val,view.ROIs(i).name,size(view.ROIs(i).coords,2));
    end
  case{'roicache'}
    % cacheVal = viewGet(view,'ROICache')
    % cacheVal = viewGet(view,'ROICache',roinum)
    if ieNotDefined('varargin')
      % if we are not passed in a specific ROI, then
      % we are doing all ROIs
      roiID = viewGet(view,'ROICacheID');
    else
      roiID = viewGet(view,'ROICacheID',varargin{1});
    end
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.roiCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.roiCache,roiID);
  case{'basecacheid'}
    % cacheID = viewGet(view,'baseCacheID')
    rotate = viewGet(view,'rotate');
    baseName = viewGet(view,'curBaseName');
    currentBase = viewGet(view,'currentBase');
    clip = viewGet(view,'baseClip',currentBase);
    gamma = viewGet(view,'baseGamma',currentBase);
    currentSlice = viewGet(view,'curSlice');
    sliceIndex = viewGet(view,'baseSliceIndex');
    % only use the corticalDepth if this is a flat
    if viewGet(view,'baseType')
      corticalDepth = viewGet(view,'corticalDepth');
    else
      corticalDepth = 0;
    end
    val = sprintf('%s_%i_%i_%i_%s_%0.2f_%0.2f',baseName,currentSlice,sliceIndex,rotate,num2str(clip),gamma,corticalDepth);
  case{'basecache'}
    % cacheVal = viewGet(view,'baseCache')
    baseID = viewGet(view,'baseCacheID');
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.baseCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.baseCache,baseID);
  case{'overlaycacheid'}
    % cacheID = viewGet(view,'overlayCacheID')
    %curSlice = viewGet(view,'curSlice');
    %analysisNum = viewGet(view,'currentAnalysis');
    %curOverlay = viewGet(view,'currentOverlay');
    %clip = viewGet(view,'overlayClip',curOverlay);
    %overlayType = viewGet(view,'overlayCtype',curOverlay);
    %overlayRange = viewGet(view,'overlayRange',curOverlay);
    % forgoe viewgets, and just grab stuff here explicitly
    % this saves about 100ms
    val = -1;
    analysisNum = view.curAnalysis;
    if ~isempty(analysisNum)
      curSlice = viewGet(view,'curSlice');
      curOverlay = view.analyses{analysisNum}.curOverlay;
      if ~isempty(curOverlay)
        % get all clips
        clip = [];
        for i = 1:length(view.analyses{analysisNum}.overlays)
          clip = [clip view.analyses{analysisNum}.overlays(i).clip];
        end
        overlayRange = view.analyses{analysisNum}.overlays(curOverlay).range;
        baseName = viewGet(view,'curBaseName');
        scanNum = viewGet(view,'curScan');
        rotate = viewGet(view,'rotate');
        alpha = viewGet(view,'alpha');
        sliceIndex = viewGet(view,'baseSliceIndex');
        % need to recalculate overlay if this is aflat
        % and the cortical depth has changed
        if viewGet(view,'baseType')
          corticalDepth = viewGet(view,'corticalDepth');
        else
          corticalDepth = 0;
        end
        % calculate string
        val = sprintf('%i_%s_%i_%i_%i_%i_%s_%s_%i_%i_%0.2f',scanNum,baseName,curSlice,sliceIndex,analysisNum,curOverlay,num2str(clip),num2str(overlayRange),rotate,alpha,corticalDepth);
      end
    end
    %    val = curSlice*analysisNum*curOverlay;
  case{'overlaycache'}
    % cacheVal = viewGet(view,'overlayCache')
    % get the overlay ID
    overlayID = viewGet(view,'overlayCacheID');
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.overlayCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.overlayCache,overlayID);
    
    
    % analysis
  case{'numberofanalyses','nanalyses','numanalyses'}
    % n = viewGet(view,'numberofAnalyses')
    % n = viewGet(view,'numberofAnalyses',[groupNum])
    if ieNotDefined('varargin')
      val = length(view.analyses);
    else
      % if this is the current group
      % then just get it
      if varargin{1} == viewGet(view,'curGroup')
        val = length(view.analyses);
      elseif (varargin{1} > 0) && (varargin{1} <= viewGet(view,'nGroups'))
        val = length(view.loadedAnalyses{varargin{1}});
      else
        val = [];
      end
    end
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
    % analysis = viewGet(view,'analysis',[analysisNum],[groupNum])
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    % if the user passed in group num, then that
    % means to check the "loadedAnalyses" for that group
    if length(varargin) <= 1
      groupNum = viewGet(view,'curGroup');
    else
      groupNum = varargin{2};
    end
    n = viewGet(view,'numberofAnalyses',groupNum);
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if groupNum == viewGet(view,'curGroup')
        val = view.analyses{analysisNum};
      else
        val = view.loadedAnalyses{groupNum}{analysisNum};
      end
    end
  case {'d'}
    % The d structure is a data structure that contains
    % outcomes of analyses that are not overlays. For
    % example, event-related analyses have a d strucutre
    % with the estimated hemodynamic responses and other
    % computed fields
    % d = viewGet(view,'d',[scanNum],[analysisNum])
    if ieNotDefined('varargin')
      scanNum = viewGet(view,'currentScan');
    else
      scanNum = varargin{1};
    end
    if isempty(scanNum)
      scanNum = viewGet(view,'currentScan');
    end
    % get analysis num
    if length(varargin) <= 1
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{2};
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if isfield(view.analyses{analysisNum},'d')
        if ~isempty(scanNum) && (scanNum > 0 ) && (scanNum <= length(view.analyses{analysisNum}.d))
          val = view.analyses{analysisNum}.d{scanNum};
        end
      end
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
  case{'numberofoverlays','numoverlays','noverlays'}
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
  case{'currentoverlay','curoverlay'}
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
    % handle if passed in a number
    if isscalar(overlayName)
      if (overlayName > 0) && (overlayName <= viewGet(view,'numOverlays'))
	val = overlayName;
      end
    elseif ~isempty(analysis) & ~isempty(analysis.overlays)
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
    % if overlayNum is a string, then convert to a number
    if isstr(overlayNum),overlayNum = viewGet(view,'overlayNum',overlayNum);end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum);
      end
    end
    if isempty(val),disp(sprintf('(viewGet:overlay) Could not find requested overlay.'));end
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
    if isstr(overlayNum),overlayNum = viewGet(view,'overlayNum',overlayNum);end
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
  case {'overlaydataval'}
    % overlaydata = viewGet(view,'overlayDataVal',x,y,s)
    % gets the value of the current overlay at the x,y,s point
    val = [];
    if ~isempty(view.curAnalysis)
      curOverlay = view.analyses{view.curAnalysis}.curOverlay;
      if ~isempty(curOverlay)
	overlay = view.analyses{view.curAnalysis}.overlays(curOverlay).data;
	if ~isempty(overlay{view.curScan})
	  val = overlay{view.curScan}(varargin{1},varargin{2},varargin{3});
	end
      end
    end
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
          if length(analysis.overlays(overlayNum).data) >= scan
            val = analysis.overlays(overlayNum).data{scan};
          end
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
  case {'alphaoverlay'}
    % overlayclip = viewGet(view,'alphaOverlay')
    analysisNum = viewGet(view,'currentAnalysis');
    overlayNum = viewGet(view,'currentOverlay',analysisNum);
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).alphaOverlay;
      end
    end
  case {'alphaoverlayexponent'}
    % overlayclip = viewGet(view,'alphaOverlayExponent')
    val = 1;
    analysisNum = viewGet(view,'currentAnalysis');
    overlayNum = viewGet(view,'currentOverlay',analysisNum);
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).alphaOverlayExponent;
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
  case {'overlayclip'}
    % overlayrange = viewGet(view,'overlayClip',[overlayNum],[analysisNum])
    % overlayrange = viewGet(view,'overlayClip',overlayNum,[])
    % overlayrange = viewGet(view,'overlayClip',[],analysisNum)
    % overlayrange = viewGet(view,'overlayClip',[],[])
    % overlayrange = viewGet(view,'overlayClip',overlayNum)
    % overlayrange = viewGet(view,'overlayClip')
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
      val = 'timecoursePlot';
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
    n = viewGet(view,'overlayNum','co');
    if n
      val = viewGet(view,'overlay',n);
    end
  case {'amp'}
    % amp = viewGet(view,'amp')
    n = viewGet(view,'overlayNum','amp');
    if n
      val = viewGet(view,'overlay',n);
    end
  case {'ph'}
    % ph = viewGet(view,'ph')
    n = viewGet(view,'overlayNum','ph');
    if n
      val = viewGet(view,'overlay',n);
    end
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
    val = view.curScan;
    if isempty(val),val = 1;end
  case {'curslice','currentslice'}
    % slice = viewGet(view,'currentSlice');
    if isfield(view.curslice,'sliceNum')
      val = view.curslice.sliceNum;
    else
      % never set in structure. Default to 1
      val = 1;
    end
  case {'curcoords','currentcoordinates'}
    % coords = viewGet(view,'currentCoordinates');
    fig = viewGet(view,'fignum');
    if ~isempty(fig)
      handles = guidata(fig);
      val = handles.coords;
    else
      val = [];
    end
  case {'alpha'}
    % alpha = viewGet(view,'alpha');
    fig = viewGet(view,'fignum');
    if ~isempty(fig)
      handles = guidata(fig);
      val = get(handles.alphaSlider,'Value');
    else
      % get alpha from analysis structure
      overlayNum = viewGet(view,'currentOverlay');
      analysisNum = viewGet(view,'currentAnalysis');
      if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
        val = view.analyses{analysisNum}.overlays(overlayNum).alpha;
      else
        val = 1;
      end
    end
  case {'overlaymin'}
    % overlayMin = viewGet(view,'overlayMin');
    % overlayMin = viewGet(view,'overlayMin',<overlayName>);
    if (length(varargin) > 0)
      overlayNum = viewGet(view,'overlayNum',varargin{1});
    else
      overlayNum = viewGet(view,'currentOverlay');
    end
    if isempty(overlayNum)
      disp(sprintf('(viewGet:overlayMin) Unknown overlay'));
    end
    % get overlayMin from analysis structure
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
      val = view.analyses{analysisNum}.overlays(overlayNum).clip(1);
    else
      val = 0;
    end
  case {'overlaymax'}
    % overlayMax = viewGet(view,'overlayMax');
    % overlayMax = viewGet(view,'overlayMax',<overlayName>);
    if (length(varargin) > 0)
      overlayNum = viewGet(view,'overlayNum',varargin{1});
    else
      overlayNum = viewGet(view,'currentOverlay');
    end
    if isempty(overlayNum)
      disp(sprintf('(viewGet:overlayMin) Unknown overlay'));
    end
    % get overlayMax from analysis structure
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
      val = view.analyses{analysisNum}.overlays(overlayNum).clip(2);
    else
      val = 1;
    end
  case {'rotate'}
    % rotate = viewGet(view,'rotate');
    baseType = viewGet(view,'baseType');
    % rotation is dependent on the base type
    % since surfaces aren't rotated in the
    % same way, they are rotated through
    % rotateSurface
    if baseType <= 1
      fig = viewGet(view,'fignum');
      if ~isempty(fig)
        handles = guidata(fig);
        val = get(handles.rotateSlider,'Value');
      else
        % otherwise gui is not running get from the structure
        curBase = viewGet(view,'curBase');
        if ~isempty(curBase)
          val = view.baseVolumes(curBase).rotate;
        else
          val = 0;
        end
      end
    else
      val = 0;
    end
  case {'rotatesurface'}
    % rotate = viewGet(view,'rotateSurface');
    baseType = viewGet(view,'baseType');
    if baseType == 2
      fig = viewGet(view,'fignum');
      if ~isempty(fig)
        handles = guidata(fig);
        % 360- is so that the anatomy rotates correct direction
        val = 360-get(handles.rotateSlider,'Value');
      else
        val = 0;
      end
    else
      val = 0;
    end
  case {'corticaldepth'}
    % rotate = viewGet(view,'rotate');
    fig = viewGet(view,'fignum');
    if ~isempty(fig)
      handles = guidata(fig);
      val = get(handles.corticalDepthSlider,'Value');
    else
      val = 0.5;
    end
  case {'sliceorientation'}
    % sliceorientation = viewGet(view,'sliceorientation');
    val = view.sliceOrientation;
  case {'cursliceoverlaycoords'}
    % overlayCoords = viewGet(view,'cursliceoverlaycoords');
    val = view.curslice.overlayCoords;
  case {'curslicebasecoords'}
    % baseCoords = viewGet(view,'curslicebasecoords');
    val = view.curslice.baseCoords;
  case {'defaultinterrogators'}
    % defaultInterrogators = viewGet(view,'defaultInterrogators');
    % see the viewSet for more info
    if isfield(MLR,'defaultInterrogators')
      val = MLR.defaultInterrogators;
    else
      val = [];
    end
  case {'colormaps'}
    % colormaps = viewGet(view,'colormaps');
    % see the viewSet for more info
    if isfield(MLR,'colormaps')
      val = MLR.colormaps;
    else
      val = [];
    end
 otherwise
    if isempty(view)
      disp(sprintf('(viewGet) No viewGet for %s',param));
      return
    else
      switch(lower(view.viewType))
        case 'volume'
          val = volumeGet(view,param,varargin{:});
        otherwise
          error('Unknown type of View.');
      end
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
    % tseriessize = viewGet(view,'tseriessize',[scanNum],[groupNum])
    % tseriessize = viewGet(view,'tseriessize',scanNum,[])
    % tseriessize = viewGet(view,'tseriessize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      nframes = viewGet(view,'nframes',s,g);
      val = [datasize nframes];
    end
    
  case {'datasize','dims'}
    % dims of single temporal frame of functional volume (same as size
    % of parameter map) for a given scan.
    %
    % datasize = viewGet(view,'datasize',[scanNum],[groupNum])
    % datasize = viewGet(view,'datasize',scanNum,[])
    % datasize = viewGet(view,'datasize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).dataSize;
    end
    
  case {'slicedims'}
    % dims of single slice and single temporal frame of functional
    % volume (same as dims of single slice of parameter map) for a
    % given scan.
    %
    % slicedims = viewGet(view,'slicedims',[scanNum],[groupNum])
    % slicedims = viewGet(view,'slicedims',scanNum,[])
    % slicedims = viewGet(view,'slicedims',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      val = datasize(1:2);
    end
    
  case {'nslices'}
    % nslices = viewGet(view,'nslices',[scanNum],[groupNum])
    % nslices = viewGet(view,'nslices',scanNum,[])
    % nslices = viewGet(view,'nslices',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      val = datasize(3);
    end
    
  case{'sliceorder'}
    % n = viewGet(view,'sliceOrder',[scanNum],[groupNum])
    % n = viewGet(view,'sliceOrder',scanNum,[])
    % n = viewGet(view,'sliceOrder',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      nslices = viewGet(view,'nslices',s,g);
      if strcmp(mrGetPref('site'),'NYU')
        if isodd(nslices)
          val = [1:2:nslices,2:2:nslices-1];
        else
          val = [2:2:nslices,1:2:nslices-1];
        end
      else
	% check for fidinfo
	fidInfo = viewGet(v,'fidInfo',s,g);
	if isempty(fidInfo) 
	  % DEFAULT, warn user and return slices in slice order
	  mrWarnDlg('(viewGet) Slice ordering is unknown for this site. Using default order: [1:nslices]. If this is incorrect, then edit viewGet sliceOrder to add the convention for your site.');
	  val = [1:nslices];
	else
	  % extract from fidInfo
	  if length(fidInfo) > 1
	    disp(sprintf('(viewGet:sliceOrder) There seems to be more than one associated FIDs with this scan, returning the sliceOrder for the first fid: %s',fidInfo{1}.fidname));
	  end
	  % get slice order
	  val = fidInfo{1}.sliceOrder;
	end
      end
    end
    
  case{'slicetimes'}
    % n = viewGet(view,'sliceTimes',[scanNum],[groupNum])
    % n = viewGet(view,'sliceTimes',scanNum,[])
    % n = viewGet(view,'sliceTimes',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    sliceOrder = viewGet(view,'sliceOrder',s,g);
    nslices = viewGet(view,'nslices',s,g);
    val = zeros(1,nslices);
    val(sliceOrder) = [0:nslices-1]/nslices;
    
  otherwise
    if strcmp(mrGetPref('verbose'),'Yes')
      dispViewGetHelp;
    end
    disp(['Invalid parameter for volume view: ',param]);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispViewGetHelp   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispViewGetHelp(param)

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

% if the user wants help on a particular command
if ~ieNotDefined('param')
  commandNum = [];
  % find the command
  for i = 1:length(commands)
    if ~isempty(findstr(lower(param),commands{i}))
      commandNum(end+1) = i;
    end
  end
  if isempty(commandNum)
    disp(sprintf('(viewGet) No viewGet for %s',param));
  else
    % display all matching commands
    for i = 1:length(commandNum)
      if (commandNum(i) > 0) && (commandNum(i) <= length(commandComments)) && ~isempty(commandComments{commandNum(i)})
        commandString = sprintf('%s %s %s',repmat('=',1,12),commands{commandNum(i)},repmat('=',1,12));
        disp(commandString);
        disp(commandComments{commandNum(i)});
      end
    end
  end
  return
end

% sort everything
[commands sortedIndex] = sort(commands);

maxlen = -inf;
% now print out
for i = 1:length(commands)
  disp(commands{i});
  lens(i) = length(commands{i});
  if length(commandComments)>=sortedIndex(i)
    disp(commandComments{sortedIndex(i)});
  end
end

maxlen = round(median(lens)+4);
nColumns = 6;
disp(sprintf('\n'));
disp('------------------------- All possible parameters ---------------------');
columnNum = 1;
currentChar = '.';
for i = 1:length(commands)
  % make a command with enough spaces
  if ~isempty(commands{i})
    % get this command, remove ',' and replace with /
    thisCommand = fixBadChars(commands{i}(2:end-1),{''',''','/'});
    % if we are switching the start letter, then put a new line
    if thisCommand(1) ~= currentChar
      currentChar = thisCommand(1);
      if columnNum ~= 1
        fprintf(1,sprintf('\n'));
      end
      fprintf(1,'(%s): ',upper(currentChar));
      columnNum = 1;
    elseif columnNum == 1
      fprintf(1,'     ');
    end
    % find out how many columns we need to print out this key word
    thisColumns = ceil((length(thisCommand)+1)/maxlen);
    % now fill up that many columns worth of space
    dispcommand = sprintf(' ');
    dispcommand = repmat(dispcommand,1,thisColumns*maxlen);
    % insert the command name
    dispcommand(1:(length(thisCommand))) = thisCommand;
    % and update what column we are on
    columnNum = columnNum + thisColumns;
    % print the command
    mrDisp(dispcommand);
    % print out new line if we have printed all columns
    if columnNum > nColumns
      fprintf(1,sprintf('\n'));
      columnNum = 1;
    end
  end
end
fprintf(1,sprintf('\n'));
disp('-----------------------------------------------------------------------');

%%%%%%%%%%%%%%%%%%
%%   getGroup   %%
%%%%%%%%%%%%%%%%%%
function g = getGroup(view,varg)

if ieNotDefined('varg')
  g = viewGet(view,'currentGroup');
else
  g = varg{1};
end
if isempty(g)
  g = viewGet(view,'currentGroup');
end
% if group is a string, then convert it to a number
if isstr(g)
  g = viewGet(view,'groupNum',g);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getScanAndGroup   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [s g] = getScanAndGroup(view,varg,param,argnum)

if ieNotDefined('argnum'), argnum = 1; end
if ieNotDefined('varg') || (length(varg) < argnum)
  s = viewGet(view,'curScan');
else
  s = varg{argnum};
end
if length(varg) > argnum
  g = varg{argnum+1};
else
  g = viewGet(view,'currentGroup');
end
if isempty(g)
  g = viewGet(view,'currentGroup');
end
if isempty(s)
  s = viewGet(view,'curScan');
end
% if group is a string, then convert it to a number
if isstr(g)
  gnum = viewGet(view,'groupNum',g);
  if isempty(gnum),disp(sprintf('(viewGet:getScanAndGroup) Unknown group: %s', g));end
  g = gnum;
end

%%%%%%%%%%%%%%%%%%%%
%%   getBaseNum   %%
%%%%%%%%%%%%%%%%%%%%
function [b baseVolume] = getBaseNum(view,varg,argnum)

if ieNotDefined('argnum'),argnum = 1;end
if ieNotDefined('varg') || (length(varg) < argnum)
  b = viewGet(view,'currentBase');
else
  b = varg{argnum};
end
if isempty(b)
  b = viewGet(view,'currentBase');
end
if nargout == 2
  baseVolume = viewGet(view,'baseVolume',b);
end

%%%%%%%%%%%%%%%%%%%
%%   getRoiNum   %%
%%%%%%%%%%%%%%%%%%%
function r= getRoiNum(view,varg,argnum)

if ieNotDefined('argnum'), argnum = 1; end
if ieNotDefined('varg') || (length(varg) < argnum)
  r = viewGet(view,'currentROI');
else
  r = varg{argnum};
  if isstr(r)
    r = viewGet(view,'roiNum',r);
  end
end
if isempty(r)
  r = viewGet(view,'currentROI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetPrependEtc    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileName = viewGetPrependEtc(view,fileName);

% prepend etc
fileName = fullfile(viewGet(view,'etcDir'),fileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetLoadStimFile    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimFile = viewGetLoadStimFile(view,stimFileName)

% prepend etc
stimFileName = fullfile(viewGet(view,'etcDir'),stimFileName);

% load this stimfile
if ~isfile(stimFileName)
  mrErrorDlg(sprintf('(viewGet:viewGetLoadStimfile): Could not find stimfile %s',stimFileName));
else
  stimFile = load(stimFileName);
end
% check to see what type it is, and set the field appropriately
if isfield(stimFile,'mylog')
  stimFile.filetype = 'eventtimes';
elseif isfield(stimFile,'stimts')
  stimFile.filetype = 'afni';
elseif isfield(stimFile,'myscreen')
  stimFile.filetype = 'mgl';
elseif isfield(stimFile,'stimvol')
  stimFile.filetype = 'stimvol';
else
  stimFile.filetype = 'unknown';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetPrependPre    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileName = viewGetPrependPre(view,fileName);

% prepend Pre
fileName = fullfile(viewGet(view,'homeDir'),'Pre',fileName);

%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetFidInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%
function fidInfo = viewGetFidInfo(view,fileName)

% prepend Pre
fileName = fullfile(viewGet(view,'homeDir'),'Pre',fileName);

% get fidInfo
[xform fidInfo] = fid2xform(fileName);
