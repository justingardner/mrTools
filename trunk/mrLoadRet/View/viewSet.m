function view = viewSet(view,param,val,varargin)
%
%   view = viewSet(view,param,val,varargin)
%
% view can be either a view structure or a viewNum which is interpreted as
% MLR.views{viewNum}. Modifies the global variable MLR.views{viewNum} as
% well as the local and returned copy.
%
% Examples:
%
% view = viewSet(view,'currentGroup',n);
%
% view = viewSet(view,'newbase',baseStructure);
% view = viewSet(view,'deletebase',baseNum);
% view = viewSet(view,'currentbase',baseNum);
% view = viewSet(view,'basemin',number,[baseNum]);
% view = viewSet(view,'basemax',number,[baseNum]);
%
% view = viewSet(view,'newanalysis',analysisStructure);
% view = viewSet(view,'deleteAnalysis',analysisNum);
% view = viewSet(view,'currentAnalysis',analysisNum);
%
% view = viewSet(view,'newoverlay',overlayStructure);
% view = viewSet(view,'currentOverlay',overlayNum);
% view = viewSet(view,'overlayname',nameString,[overlayNum]);
% view = viewSet(view,'overlaycmap',cmapName,[overlayNum]);
% view = viewSet(view,'overlaymin',number,[overlayNum]);
% view = viewSet(view,'overlaymax',number,[overlayNum]);
%
% view = viewSet(view,'newROI',roiStructure);
% view = viewSet(view,'deleteroi',roiNum);
% view = viewSet(view,'currentROI',roiNum);
% view = viewSet(view,'roiCoords',array,[roiNum]);
%
% view = viewSet(view,'newScan',fileName,description);
% view = viewSet(view,'deleteScan',scanNum);
%
% Unlike viewGet, very few of the viewSet arguments are optional although
% there are a few (some examples of which are listed above. Specifically,
% 1hen setting the subfields of baseVolumes, overlays, and ROIs, there is
% an optional argument for specifying the baseNumber, overlayNumber, or
% roiNumber. Default (if that argument is not passed) is to use the current
% base, overlay, or ROI. Please follow this convention if you add to this
% function.
%
% See comments in the source code for the full set of parameters and for
% optional arguments.
%
% 6/2004 djh

mrGlobals

if ieNotDefined('param')
  dispViewSetHelp;
  return
end

if ieNotDefined('view'), mrErrorDlg('No view specified.'); end
if ieNotDefined('param'), mrErrorDlg('No parameter specified'); end
if ieNotDefined('val'), val = []; end

% view can be either a view structure or a view number. Either way, make
% sure that the view passed is consistent with the global variable
% MRL.views{viewNum}
if isnumeric(view)
  viewNum = view;
  view = MLR.views{viewNum};
else
  viewNum = view.viewNum;
  MLR.views{viewNum} = view;
end
if (viewNum < 1) | (viewNum > length(MLR.views))
  mrErrorDlg('Invalid viewNum.');
end

switch lower(param)

  % Session
  case {'homedir','homedirectory','sessiondirectory'}
    % view = viewSet(view,'homedir',string);
    MLR.homeDir = val;

    % View
  case {'viewtype','type'}
    % view = viewSet(view,'viewtype',string);
    view.viewType = val;

    % -------------------------------------------
    % Group

  case{'currentgroup','curgroup'}
    % view = viewSet(view,'currentGroup',n);
    % view = viewSet(view,'currentGroup','groupName')
    if isstr(val)
      % convert to a groupNum if we got passed in a groupname
      groupNum = viewGet(view,'groupNum',val);
      if isempty(groupNum)
	disp(sprintf('(viewSet:curGroup) No group %s found',val));
      end
      val = groupNum;
    end
    if (view.curGroup ~= val)
      % save loaded analysis if there are any
      view = viewSet(view,'loadedAnalyses',view.analyses,view.curGroup);
      % save the current scan number
      view = viewSet(view,'groupScanNum',viewGet(view,'curScan'),view.curGroup);
      MLR.views{view.viewNum} = view;
      % set the current group
      view.curGroup = val;
      view.analyses = [];
      view.curAnalysis = [];
      % Update the gui
      mlrGuiSet(view,'group',val);
      nScans = viewGet(view,'nScans',val);
      mlrGuiSet(view,'nScans',nScans);
      scanNum = viewGet(view,'groupScanNum',view.curGroup);
      mlrGuiSet(view,'scan',min(nScans,scanNum));
      mlrGuiSet(view,'analysis',1);
      mlrGuiSet(view,'analysisPopup',{'none'});
      mlrGuiSet(view,'overlay',1);
      mlrGuiSet(view,'overlayPopup',{'none'});
      % update the interrogator
      if isfield(MLR,'interrogator') && (view.viewNum <=length(MLR.interrogator))
        mrInterrogator('updateInterrogator',view.viewNum,viewGet(view,'interrogator'));

      end
      % load up saved analysis if they exist
      loadedAnalyses = viewGet(view,'loadedAnalyses',val);
      for i = 1:length(loadedAnalyses)
	view = viewSet(view,'newAnalysis',loadedAnalyses{i});
      end
      % delete the analyses from the loaded cache
      view = viewSet(view,'loadedAnalyses',{},view.curGroup);
    end

  case {'loadedanalyses'}
    % view = viewSet(view,'loadedAnalyses',analyses,groupNum);
    % this is used to remember the analyses that were loaded
    % for a group for switching between groups
    groupNum = varargin{1};
    if (groupNum >=1 ) && (groupNum <= viewGet(view,'numGroups'))
      view.loadedAnalyses{groupNum} = val;
    end
  case {'groupscannum'}
    % view = viewSet(view,'groupScanNum',scanNum,groupNum);
    % this is used to remember which scan we were on when
    % we switch between groups
    groupNum = varargin{1};
    if (groupNum >=1 ) && (groupNum <= viewGet(view,'numGroups'))
      if ~isempty(val)
	view.groupScanNum(groupNum) = val;
      end
    end
  case{'groupname'}
    % view = viewSet(view,'currentGroup',string);
    n = viewGet(view,'groupnum',val);
    view = viewSet(view,'currentGroup',n);

  case {'newgroup'}
    % view = viewSet(view,'newGroup',groupName);
    newgroup.name = val;
    newgroup.scanParams = [];
    newgroup.auxParams = [];
    % Add it to MLR.groups
    pos = length(MLR.groups)+1;
    if (pos == 1)
      MLR.groups = newgroup;
    else
      MLR.groups(pos) = newgroup;
    end
    % Update the GUIs of all views
    stringList = {MLR.groups(:).name};
    for v = 1:length(MLR.views)
      if isview(MLR.views{v})
        mlrGuiSet(MLR.views{v},'groupPopup',stringList);
      end
    end
    % Save mrSession
    saveSession;

  case {'deletegroup'}
    % view = viewSet(view,'deleteGroup',groupNum);
    groupnum = val;
    groupName = viewGet(view,'groupName',groupnum);
    if strcmp(groupName,'Raw')
      mrWarnDlg('Cannot delete Raw group');
      return
    end
    % Remove it
    numgroups = viewGet(view,'numberofgroups');
    MLR.groups = MLR.groups(groupnum ~= [1:numgroups]);
    % Update the view
    groupNames = viewGet(view,'groupNames');
    if isempty(groupNames)
      groupNames = {'none'};
    end
    curgroup = viewGet(view,'currentGroup');
    if (curgroup > groupnum)
      view = viewSet(view,'currentGroup',curgroup-1);
    elseif (groupnum == curgroup)
      view = viewSet(view,'currentGroup',1);
    end
    mlrGuiSet(view,'groupPopup',groupNames);
    % Update all the other views
    for v = find(view.viewNum ~= [1:length(MLR.views)])
      if isview(MLR.views{v})
        curgroup = viewGet(v,'currentGroup');
        if (curgroup > groupnum)
          viewSet(v,'currentGroup',curgroup-1);
        elseif (groupnum == curgroup)
          viewSet(v,'currentGroup',1);
        end
        mlrGuiSet(MLR.views{v},'groupPopup',groupNames);
      end
    end
    % Save mrSession
    saveSession;

    % -------------------------------------------
    % Scan

  case {'stimfilename'}
    % view = viewSet(view,'stimFilename',stimFilename,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isempty(val)
        MLR.groups(g).auxParams(s).stimFileName = [];
      else
        MLR.groups(g).auxParams(s).stimFileName = val;
      end
    end

  case {'eyepos'}
    % view = viewSet(view,'eyepos',eyeposFilename,eyeposNum,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin(2:end),varargin{1});
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isempty(val)
        MLR.groups(g).auxParams(s).eyeposFileName = [];
        MLR.groups(g).auxParams(s).eyeposNum = [];
      else
        MLR.groups(g).auxParams(s).eyeposFileName = val;
        MLR.groups(g).auxParams(s).eyeposNum = varargin{1};
      end
    end

  case {'spikeinfo'}
    % view = viewSet(view,'spikeinfo',spikeinfo,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).auxParams(s).spikeinfo = val;
    end
  case {'niftihdr'}
    % view = viewSet(view,'niftiHdr',hdr,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr = val;
    end

  case {'scanparams'}
    % view = viewSet(view,'scanParams',scanParams,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    if isscan(val)
      nscans = viewGet(view,'nscans',g);
      if (nscans >= s) & (s > 0)
        MLR.groups(g).scanParams(s) = val;
      end
    else
      disp(sprintf('(viewSet) Invalid scanParams'));
    end

  case {'scanxform','sform'}
    % view = viewSet(view,'sform',sform,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr.sform44 = val;
    end
  case {'sformcode'}
    % view = viewSet(view,'sformcode',sformcode,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr.sform_code = val;
    end
  case {'newscan','updatescan'}
    % view = viewSet(view,'newScan',scanParams);
    % val must be a scanParams structure with scanParams.fileName set
    % to a nifti file in the tseries directory. Only the filename,
    % description, junkframes, nframes, originalFileName and
    % originalGroupName fields are used; the other fields in scanParams
    % are reset to insure consistency with the nifti file.
    %
    % view = viewSet(view,'updateScan',scanParams,scanNum);
    % Check for tseries file and (re-)build scanParams to insure
    % consistency with the nifti file.
    scanParams = orderfields(val);
    if ~isfield(scanParams,'fileName')
      mrErrorDlg(['scanParams.fileName must be specified']);
    end
    tseriesdir = viewGet(view,'tseriesdir');
    path = fullfile(tseriesdir,scanParams.fileName);
    if ~exist(path,'file')
      mrErrorDlg(['Tseries file not found: ', path]);
    end
    hdr = cbiReadNiftiHeader(path);
    if ~isfield(scanParams,'description')
      scanParams.description = '';
    end
    scanParams.fileType = 'Nifti';
    scanParams.niftiHdr = hdr;
    scanParams.voxelSize = hdr.pixdim([2,3,4])';
    scanParams.totalFrames = hdr.dim(5);
    if ~isfield(scanParams,'junkFrames')
      scanParams.junkFrames = 0;
    end
    if ~isfield(scanParams,'nFrames')
      scanParams.nFrames = scanParams.totalFrames;
    end
    scanParams.dataSize = hdr.dim([2,3,4])';
    scanParams.framePeriod = hdr.pixdim(5)/1000;
    % check for field that tells what the original
    % filename is. if it is not passed in
    % then just set it to the current name and group
    if ~isfield(scanParams,'originalFileName')
      scanParams.originalFileName{1} = scanParams.fileName;
      scanParams.originalGroupName{1} = viewGet(view,'groupName');
    elseif isstr(scanParams.originalFileName)
      originalFileName{1} = scanParams.originalFileName;
      scanParams.originalFileName = scanParams.originalFileName;
    end
    % Check for required fields and sort them
    [check scanParams] = isscan(scanParams);
    if ~check
      mrErrorDlg(['Invalid scanParams.']);
    end
    % Add scanParams to group
    curgroup = viewGet(view,'currentGroup');
    curGroupName = viewGet(view,'groupName',curgroup);
    nscans = viewGet(view,'nscans');
    if strcmp(lower(param),'newscan')
      if (nscans < 1)
        MLR.groups(curgroup).scanParams = scanParams;
        % create and auxparams with a default field, using
        % struct([]) to create an empty structure do not
        % behave well.
        MLR.groups(curgroup).auxParams = struct('auxParams',1);
      else
        MLR.groups(curgroup).scanParams(nscans+1) = scanParams;
        % get an empty auxParams
        if isfield(MLR.groups(curgroup),'auxParams') && (length(MLR.groups(curgroup).auxParams) >= nscans)
          auxParams = MLR.groups(curgroup).auxParams(nscans);
          auxParamFields = fieldnames(auxParams);
          for aNum = 1:length(auxParamFields)
            auxParams.(auxParamFields{aNum}) = [];
          end
          % and add it to the list
          MLR.groups(curgroup).auxParams(nscans+1) = auxParams;
        end
      end
      % Reconcile analysis params and overlay data/params with tseries
      % files, adding empty data and default params to new scan.
      view = reconcileAnalyses(view);
      % Update GUI
      mlrGuiSet(view,'nScans',nscans+1);
      % Do the same for all the other views
      for v = find(view.viewNum ~= [1:length(MLR.views)])
        if isview(MLR.views{v})
          group = viewGet(v,'currentGroup');
          if (group == curgroup)
            mlrGuiSet(v,'nScans',nscans+1);
          end
        end
      end
    elseif strcmp(lower(param),'updatescan')
      % get the scan number we are overwriting
      if length(varargin) ~= 1
        error('viewGet: updatescan, requires scan number.');
      else
        scanNum = varargin{1};
      end
      % check to make sure it exists
      if isempty(viewGet(view,'scanParams',scanNum))
        error('viewGet: Scan #%i does not exist',scanNum);
      end
      % change it
      MLR.groups(curgroup).scanParams(scanNum) = scanParams;
    end

    % Save mrSession
    saveSession;

  case {'deletescan'}
    % view = viewSet(view,'deleteScan',scanNum);
    scannum = val;
    numscans = viewGet(view,'nscans');
    curgroup = viewGet(view,'currentGroup');
    curGroupName = viewGet(view,'groupName',curgroup);
    curscan = viewGet(view,'currentScan');
    % Remove it and reset curscan
    MLR.groups(curgroup).scanParams = MLR.groups(curgroup).scanParams(scannum ~= [1:numscans]);
    % Reconcile analysis params and overlay data/params with tseries
    % files, removing overlay data corresponding to this scan.
    view = reconcileAnalyses(view);
    % Update GUI
    mlrGuiSet(view,'nScans',max(numscans-1,0));
    if (curscan >= scannum)
      mlrGuiSet(view,'scan',max(curscan-1,0));
    end
    % Do the same for all the other views
    for v = find(view.viewNum ~= [1:length(MLR.views)])
      if isview(MLR.views{v})
        group = viewGet(v,'currentGroup');
        if (group == curgroup)
          mlrGuiSet(v,'nScans',max(numscans-1,0));
        end
      end
    end
    % Save mrSession
    saveSession;

    % -------------------------------------------
    % Base anatomy

  case{'newbase'}
    % view = viewSet(view,'newbase',baseStructure);
    %
    % val must be a structure with fields
    % - name: string
    % - data: [x y z] array
    % - range: min and max values
    % - clip: min and max to be displayed
    % - hdr: nifti header
    % - permutation matrix: keeps track of slice orientation (see
    %   loadAnat)

    % Check that val has the required fields
    [check baseAnatomy] = isbase(val);
    if ~check
      mrErrorDlg('Invalid base anatomy');
    end
    % If baseVolume.name already exists then replace the existing one
    % with this new one. Otherwise, add it to the end of the baseVolumes
    % list.
    newBaseName = baseAnatomy.name;
    newBaseNum = [];
    numBaseVolumes = viewGet(view,'numberofBaseVolumes');
    for num = 1:numBaseVolumes
      baseName = viewGet(view,'basename',num);
      if strcmp(newBaseName,baseName)
        newBaseNum = num;
      end
    end
    if isempty(newBaseNum)
      newBaseNum = numBaseVolumes + 1;
    end
    % Add it to the list of baseVolumes
    if (newBaseNum == 1)
      view.baseVolumes = baseAnatomy;
    else
      view.baseVolumes(newBaseNum) = baseAnatomy;
    end
    % Set it to be the current base Volume
    view = viewSet(view,'curBase',newBaseNum);
    % clear the caches of any reference to a base with
    % the same name (if it is reloaded, the base may
    % have changed its xform for instance and so the
    % caches need to be recomputed)
    view = viewSet(view,'baseCache','clear',baseAnatomy.name);
    view = viewSet(view,'overlayCache','clear',baseAnatomy.name);
    view = viewSet(view,'roiCache','clear',baseAnatomy.name);
    % Update the gui
    stringList = {view.baseVolumes(:).name};
    mlrGuiSet(view,'basePopup',stringList);

  case {'deletebase'}
    % view = viewSet(view,'deletebase',baseNum);
    basenum = val;
    numbases = viewGet(view,'numberofbasevolumes');
    curbase = viewGet(view,'currentBase');
    % Remove it and reset currentBase
    view.baseVolumes = view.baseVolumes(basenum ~= [1:numbases]);
    if (numbases > 1)
      if (curbase > basenum)
        view = viewSet(view,'currentBase',curbase-1);
      elseif (basenum == curbase)
        view = viewSet(view,'currentBase',1);
      end
    else
      view.curBase = [];
    end
    % Update the gui
    stringList = {view.baseVolumes(:).name};
    if isempty(stringList)
      stringList = {'none'};
    end
    mlrGuiSet(view,'basePopup',stringList);

  case{'currentbase','curbase','curanat'}
    % view = viewSet(view,'currentbase',baseNum);
    baseNum = val;
    numBases = viewGet(view,'numberofBaseVolumes');
    % first save rotation, curSlice and sliceIndex
    % for this base image
    curBase = viewGet(view,'curBase');
    rotate = viewGet(view,'rotate');
    curSlice = viewGet(view,'curSlice');
    sliceOrientation = viewGet(view,'sliceOrientation');
    % set the current state of the gui in the base
    if (curBase > 0) & (curBase <= numBases)
      view.baseVolumes(curBase).rotate = rotate;
      view.baseVolumes(curBase).curSlice = curSlice;
      view.baseVolumes(curBase).sliceOrientation = sliceOrientation;
    end
    % now switch to new base
    if (baseNum > 0) & (baseNum <= numBases)
      view.curBase = baseNum;
      % update popup menu and sliders
      mlrGuiSet(view,'baseVolume',baseNum);
      % set basedims
      baseDims = viewGet(view,'baseDims',baseNum);
      mlrGuiSet(view,'baseDims',baseDims);
      % set base min and max sliders
      baseClip = viewGet(view,'baseClip',baseNum);
      baseRange = viewGet(view,'baseRange',baseNum);
      mlrGuiSet(view,'baseMinRange',baseRange);
      mlrGuiSet(view,'baseMaxRange',baseRange);
      mlrGuiSet(view,'baseMin',baseClip(1));
      mlrGuiSet(view,'baseMax',baseClip(2));
      % get the last know settings for this base
      baseCurSlice = viewGet(view,'baseCurSlice',baseNum);
      baseSliceOrientation = viewGet(view,'baseSliceOrientation',baseNum);
      baseRotate = viewGet(view,'baseRotate',baseNum);
      % set the slice orientation if there is a valid one saved
      if ~isempty(baseSliceOrientation)
	viewSet(view,'sliceOrientation',baseSliceOrientation);
      end
      % update nSlices and reset slice to be within range
      sliceIndex = viewGet(view,'basesliceindex',baseNum);
      nSlices = baseDims(sliceIndex);
      curSlice = viewGet(view,'curSlice');
      % if the base has a current slice set, then use that
      if isempty(baseCurSlice) || (baseCurSlice > nSlices)
	mlrGuiSet(view,'slice',max(1,min(curSlice,nSlices)));
      else
	mlrGuiSet(view,'slice',baseCurSlice);
      end
      if ~isempty(baseRotate)
	mlrGuiSet(view,'rotate',baseRotate);
      end
      mlrGuiSet(view,'nSlices',nSlices);
    end

  case 'basemin'
    % view = viewSet(view,'basemin',number,[baseNum]);
    curBase = viewGet(view,'currentBase');
    if ~isempty(varargin)
      baseNum = varargin{1};
    else
      baseNum = curBase;
    end
    if ~isempty(baseNum)
      view.baseVolumes(baseNum).clip(1) = val;
    end
    if (baseNum == curBase)
      mlrGuiSet(view,'baseMin',val);
    end

  case 'basemax'
    % view = viewSet(view,'basemax',number,[baseNum]);
    curBase = viewGet(view,'currentBase');
    if ~isempty(varargin)
      baseNum = varargin{1};
    else
      baseNum = curBase;
    end
    if ~isempty(baseNum)
      view.baseVolumes(baseNum).clip(2) = val;
    end
    if (baseNum == curBase)
      mlrGuiSet(view,'baseMax',val);
    end

  case 'baserange'
    % view = viewSet(view,'baserange',[min max],[baseNum]);
    range = val;
    curBase = viewGet(view,'currentBase');
    if ieNotDefined('varargin')
      baseNum = curBase;
    else
      baseNum = varargin{1};
    end
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).range = range;
      if (baseNum == curBase)
        mlrGuiSet(view,'baseMinRange',range);
        mlrGuiSet(view,'baseMaxRange',range);
        clip = view.baseVolumes(baseNum).clip;
        mlrGuiSet(view,'baseMin',clip(1));
        mlrGuiSet(view,'baseMax',clip(2));
      end
    end

  case{'basename'}
    % view = viewSet(view,'basename',nameString,[baseNum]);
    if ieNotDefined('varargin')
      baseNum = viewGet(view,'curBase');
    else
      baseNum = varargin{1};
    end
    if ~isempty(baseNum)
      view.baseVolumes(baseNum).name = val;
      mlrGuiSet(view,'basePopup',viewGet(view,'baseNames'));
    end

    % -------------------------------------------
    % Analysis

  case{'newanalysis'}
    % view = viewSet(view,'newanalysis',analysisStructure);
    %
    % val must be a structure with fields
    % - name: string
    % - groupName: string group name
    % - function: string function name by which it was computed
    % - reconcileFunction: string function name that reconciles params
    %   and data with tseries file names.
    %      [newparams,newdata] = reconcileFunction(groupName,params,data)
    % - guiFunction: string function name that allows user to specify
    %   or change params.
    %      params = guiFunction('groupName',groupName,'params',params)
    % - params: structure specifying arguments to function
    %      To recompute: view = function(view,params)
    % - overlays: struct array of overlays
    % - curOverlay: integer corresponding to currently selected overlay
    % - date: specifies when it was computed

    % Check that is has the required fields
    [check analysis] = isanalysis(val);
    if ~check
      mrErrorDlg('Invalid analysis');
    end
    % Error if groupNames don't match
    curGroupNum = viewGet(view,'currentGroup');
    curGroupName = viewGet(view,'groupName',curGroupNum);
    if ~strcmp(analysis.groupName,curGroupName)
      mrErrorDlg(['analysis is incompatible with group: ',curGroupName]);
    end
    % Reconcile analysis params with group
    if isfield(analysis,'d')
      [analysis.params analysis.d] = feval(analysis.reconcileFunction,curGroupName,analysis.params,analysis.d);
    else
      analysis.params = feval(analysis.reconcileFunction,curGroupName,analysis.params);
    end
    % If analysis.name already exists then replace the existing one with
    % this new one. Otherwise, add it to the end of the analysis list.
    newAnalysisName = analysis.name;
    newAnalysisNum = [];
    nAnalyses = viewGet(view,'numberofAnalyses');
    for num = 1:nAnalyses
      analysisName = viewGet(view,'analysisName',num);
      if strcmp(newAnalysisName,analysisName)
        newAnalysisNum = num;
      end
    end
    if isempty(newAnalysisNum)
      newAnalysisNum = nAnalyses + 1;
    end
    % Grab overlays to reinstall later
    overlays = analysis.overlays;
    analysis.overlays = [];
    % Add analysis to the list of analyses
    view.analyses{newAnalysisNum} = analysis;
    % Update the gui
    mlrGuiSet(view,'analysisPopup',viewGet(view,'analysisNames'));
    % Set it to be the current analysis
    view = viewSet(view,'curanalysis',newAnalysisNum);
    % Reinstall overlays
    for n = 1:length(overlays)
      view = viewSet(view,'newOverlay',overlays(n));
    end
    % update the interrogator
    if isfield(MLR,'interrogator') && (view.viewNum <=length(MLR.interrogator))
      mrInterrogator('updateInterrogator',view.viewNum,viewGet(view,'interrogator'));
    end
    % completely clear the overlay cache
    view = viewSet(view,'overlayCache','init');
    % Set current overlay
    if isfield(val,'curOverlay')
      view = viewSet(view,'currentOverlay',val.curOverlay);
    end

  case {'interrogator'}
    % view = viewSet(view,'interrogator',interrogator);
    numAnalyses = viewGet(view,'numberofAnalyses');
    curAnalysis = viewGet(view,'currentAnalysis');
    if ~isempty(curAnalysis) & (curAnalysis >= 1) & (curAnalysis <= numAnalyses)
      curOverlay = viewGet(view,'currentOverlay');
      if ~isempty(curOverlay) & (curOverlay >= 1) & (curOverlay <= viewGet(view,'numOverlays'))
	view.analyses{curAnalysis}.overlays(curOverlay).interrogator = val;
      end
    end
  case {'deleteanalysis'}
    % view = viewSet(view,'deleteAnalysis',analysisNum);
    analysisNum = val;
    numAnalyses = viewGet(view,'numberofAnalyses');
    curAnalysis = viewGet(view,'currentAnalysis');
    % Remove it and reset currentAnalysis
    view.analyses = {view.analyses{analysisNum ~= [1:numAnalyses]}};
    if (numAnalyses > 1)
      if (curAnalysis > analysisNum)
        view = viewSet(view,'currentAnalysis',curAnalysis-1);
      elseif (analysisNum == curAnalysis)
        view = viewSet(view,'currentAnalysis',1);
      end
    else
      view = viewSet(view,'currentAnalysis',[]);
    end
    % Update the gui
    stringList = viewGet(view,'analysisNames');
    if isempty(stringList)
      stringList = {'none'};
    end
    mlrGuiSet(view,'analysisPopup',stringList);

  case {'currentanalysis','curanalysis'}
    % view = viewSet(view,'currentAnalysis',analysisNum);
    analysisNum = val;
    numAnalyses = viewGet(view,'numberofAnalyses');
    % Set curAnalysis field and update the gui
    if (analysisNum > 0) & (analysisNum <= numAnalyses)
      view.curAnalysis = analysisNum;
      mlrGuiSet(view,'analysis',analysisNum);
    else
      view.curAnalysis = [];
      mlrGuiSet(view,'analysis',1);
    end
    % Update overlay popup
    stringList = viewGet(view,'overlayNames',analysisNum);
    if isempty(stringList)
      stringList = {'none'};
    end
    mlrGuiSet(view,'overlayPopup',stringList);
    % Set current overlay
    curOverlay = viewGet(view,'currentOverlay',analysisNum);
    view = viewSet(view,'currentOverlay',curOverlay);
    % update the interrogator
    if isfield(MLR,'interrogator') && (view.viewNum <=length(MLR.interrogator))
      mrInterrogator('updateInterrogator',view.viewNum,viewGet(view,'interrogator'));

    end

  case{'analysisname'}
    % view = viewSet(view,'analysisname',nameString,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      view.analyses{analysisNum}.name = val;
      mlrGuiSet(view,'analysisPopup',viewGet(view,'analysisNames'));
    end
  case{'analysisguifunction'}
    % view = viewSet(view,'analysisGUIFunction',GUIFunctionString,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      view.analyses{analysisNum}.guiFunction = val;
    end
  case{'analysismergefunction'}
    % view = viewSet(view,'analysisMergeFunction',mergeFunctionString,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      view.analyses{analysisNum}.mergeFunction = val;
    end
  case{'analysisreconcilefunction'}
    % view = viewSet(view,'analysisReconcileFunction',reconcileFunctionString,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      view.analyses{analysisNum}.reconcileFunction = val;
    end

  case{'analysistype'}
    % view = viewSet(view,'analysisType',typeString,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      view.analyses{analysisNum}.type = val;
    end

    % -------------------------------------------
    % Overlay (parameter map)

  case{'newd'}
    % view = viewSet(view,'newd',dStructure,scanNum,[analysisNum]);
    if ieNotDefined('varargin')
      mrErrorDlg('(viewSet:newd) No passed in scanNum');
    end
    scanNum = varargin{1};
    if length(varargin) == 1
      anum = viewGet(view,'curAnalysis');
    else
      anum = varargin{2};
    end
    if (anum > 0) && (anum <= viewGet(view,'numAnalyses'))
      if (scanNum > 0) && (scanNum <= viewGet(view,'nScans'))
	view.analyses{anum}.d{scanNum} = val;
      else
	mrErrorDlg(sprintf('(viewSet:newd) Scan %i out of range',scanNum));
      end
    else
      mrErrorDlg(sprintf('(viewSet:newd) Analysis %i out of range',anum));
    end
  case{'newoverlay'}
    % view = viewSet(view,'newoverlay',overlayStructure,[analysisNum]);
    %
    % val must be a structure with fields
    % - name: string
    % - groupName: string group name
    % - function: string function name by which it was computed
    % - reconcileFunction: string function name that reconciles params
    %   and data with tseries file names.
    %      [newparams,newdata] = reconcileFunction(groupName,params,data)
    % - params: structure specifying arguments to function
    %      To recompute: view = function(view,params)
    % - interrogator: function that gets called when you click on the
    %   overlay (e.g., to produce a time series plot).
    % - data: cell array of [x y z] arrays
    % - date: specifies when it was computed, typically generated using
    %   datestr(now)
    % - range: [min max] values
    % - clip: [min max] to be displayed/thresholded
    % - colormap: 256x3 array of RGB values
    % - alpha: transparency value for alphaSlider

    % Check that it has the required fields
    [check overlay] = isoverlay(val);
    if ~check
      mrErrorDlg('Invalid overlay');
    end
    % groupName, groupNum, nScans
    groupName = overlay.groupName;
    groupNum = viewGet(view,'groupnum',groupName);
    nScans = viewGet(view,'nScans',groupNum);
    % analysisNum and analysisName
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    analysisName = viewGet(view,'analysisName',analysisNum);
    analysis = viewGet(view,'analysis',analysisNum);
    % Error if groupNames don't match
    if ~strcmp(overlay.groupName,analysis.groupName)
      mrErrorDlg(['(viewSet) overlay is incompatible with group: ',groupName]);
    end
    % Reconcile overlay data and params with tseries files, reordering
    % them as necessary.
    [params,data] = feval(overlay.reconcileFunction,groupName,overlay.params,overlay.data);
    overlay.params = params;
    overlay.data = data;
    % Check that it has the correct number of scans
    if (length(overlay.data) ~= nScans)
      %mrErrorDlg('Invalid overlay, incorrect number of scans.');
      if isempty(data)
	disp('(viewSet) No data in overlay');
      else
	disp('(viewSet) Invalid overlay, incorrect number of scans.');
      end
      return;
    end
    % Check that the data arrays for each scan are the correct size
    for s = 1:nScans
      overlayDims = viewGet(view,'dims',s,groupNum);
      overlaySize = size(overlay.data{s});
      % if the overlay only has one slice, then we need
      % to add the 1 to the last dimension
      if length(overlaySize) == 2
	overlaySize(3) = 1;
      end
      if (~isempty(overlay.data{s})) & ~isequal(overlaySize,overlayDims)
        mrErrorDlg('(viewSet) Invalid overlay, incorrect data array size.');
      end
    end
    % If overlay.name already exists then replace the existing one with
    % this new one. Otherwise, add it to the end of the overlays list.
    newOverlayName = overlay.name;
    newOverlayNum = [];
    nOverlays = viewGet(view,'numberofOverlays',analysisNum);
    for num = 1:nOverlays
      overlayName = viewGet(view,'overlayName',num,analysisNum);
      if strcmp(newOverlayName,overlayName)
        newOverlayNum = num;
      end
    end
    if isempty(newOverlayNum)
      newOverlayNum = nOverlays + 1;
    end
    % Add it to the list of overlays
    if (newOverlayNum == 1) & isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays = overlay;
    else
      view.analyses{analysisNum}.overlays(newOverlayNum) = overlay;
    end
    % Update the gui
    mlrGuiSet(view,'overlayPopup',viewGet(view,'overlayNames',analysisNum));
    % Set it to be the current overlay
    view = viewSet(view,'curOverlay',newOverlayNum);

  case {'deleteoverlay'}
    % view = viewSet(view,'deleteoverlay',overlayNum,[analysisNum]);
    overlayNum = val;
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    numoverlays = viewGet(view,'numberofoverlays',analysisNum);
    curoverlay = viewGet(view,'currentOverlay',analysisNum);
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      % Remove it and reset currentOverlay
      view.analyses{analysisNum}.overlays = ...
        view.analyses{analysisNum}.overlays(overlayNum ~= [1:numoverlays]);
      if (numoverlays > 1)
        if (curoverlay > overlayNum)
          view = viewSet(view,'currentOverlay',curoverlay-1);
        elseif (overlayNum == curoverlay)
          view = viewSet(view,'currentOverlay',1);
        end
      else
        view = viewSet(view,'currentOverlay',[]);
      end
      % Update the gui
      stringList = viewGet(view,'overlayNames',analysisNum);
      if isempty(stringList)
        stringList = {'none'};
      end
      mlrGuiSet(view,'overlayPopup',stringList);
    end

  case{'currentoverlay','curoverlay'}
    % view = viewSet(view,'currentOverlay',overlayNum);
    overlayNum = val;
    analysisNum = viewGet(view,'currentAnalysis');
    numAnalyses = viewGet(view,'numberofAnalyses');
    numOverlays = viewGet(view,'numberofOverlays',analysisNum);
    if (analysisNum > 0) & (analysisNum <= numAnalyses)
      if (overlayNum > 0) & (overlayNum <= numOverlays)
        % Set the curOverlay field
        view.analyses{analysisNum}.curOverlay = overlayNum;
        % Update the gui
        mlrGuiSet(view,'overlay',overlayNum);
        % set overlay min and max sliders
        overlayClip = viewGet(view,'overlayClip',overlayNum,analysisNum);
        overlayRange = viewGet(view,'overlayRange',overlayNum,analysisNum);
        if isempty(overlayClip)
          overlayClip = [0 1];
        end
        if isempty(overlayRange)
          overlayRange = [0 1];
        end
        mlrGuiSet(view,'overlayMinRange',overlayRange);
        mlrGuiSet(view,'overlayMaxRange',overlayRange);
        mlrGuiSet(view,'overlayMin',overlayClip(1));
        mlrGuiSet(view,'overlayMax',overlayClip(2));
        % set overlay alpha slider
        alpha = viewGet(view,'overlayAlpha',overlayNum,analysisNum);
        if isempty(alpha)
          alpha = 1;
        end
        mlrGuiSet(view,'overlayAlpha',alpha);
      else
        view.analyses{analysisNum}.curOverlay = [];
        mlrGuiSet(view,'overlay',1);
      end
    else
      mlrGuiSet(view,'overlay',1);
    end
    % update the interrogator
    if isfield(MLR,'interrogator') && (view.viewNum <=length(MLR.interrogator))
      mrInterrogator('updateInterrogator',view.viewNum,viewGet(view,'interrogator'));

    end

  case 'overlayname'
    % view = viewSet(view,'overlayname',nameString,[overlayNum]);
    if ieNotDefined('varargin')
      overlayNum = viewGet(view,'currentOverlay');
    else
      overlayNum = varargin{1};
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays(overlayNum).name = val;
      mlrGuiSet(view,'overlaypopup',viewGet(view,'overlayNames',analysisNum));
    end

  case 'overlaycmap'
    % view = viewSet(view,'overlaycmap',cmapName,[overlayNum]);
    if ieNotDefined('varargin')
      overlayNum = viewGet(view,'currentOverlay');
    else
      overlayNum = varargin{1};
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      evalstr = [val,'(256)'];
      view.analyses{analysisNum}.overlays(overlayNum).colormap = eval(evalstr);
    end

  case 'overlaymin'
    % view = viewSet(view,'overlaymin',number,[overlayNum]);
    curOverlay = viewGet(view,'currentOverlay');
    if ieNotDefined('varargin')
      overlayNum = curOverlay;
    else
      overlayNum = varargin{1};
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays(overlayNum).clip(1) = val;
      if (overlayNum == curOverlay)
        mlrGuiSet(view,'overlayMin',val);
      end
    end

  case 'overlaymax'
    % view = viewSet(view,'overlaymax',number,[overlayNum]);
    curOverlay = viewGet(view,'currentOverlay');
    if ieNotDefined('varargin')
      overlayNum = curOverlay;
    else
      overlayNum = varargin{1};
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays(overlayNum).clip(2) = val;
      if (overlayNum == curOverlay)
        mlrGuiSet(view,'overlayMax',val);
      end
    end

  case 'overlayrange'
    % view = viewSet(view,'overlayrange',[min max],[overlayNum]);
    range = val;
    curOverlay = viewGet(view,'currentOverlay');
    if ieNotDefined('varargin')
      overlayNum = curOverlay;
    else
      overlayNum = varargin{1};
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays(overlayNum).range = range;
      if (overlayNum == curOverlay)
        mlrGuiSet(view,'overlayMinRange',range);
        mlrGuiSet(view,'overlayMaxRange',range);
        clip = view.analyses{analysisNum}.overlays(overlayNum).clip;
        mlrGuiSet(view,'overlayMin',clip(1));
        mlrGuiSet(view,'overlayMax',clip(2));
      end
    end

  case 'alpha'
    % view = viewSet(view,'alpha',number,[overlayNum]);
    curOverlay = viewGet(view,'currentOverlay');
    if ~isempty(varargin)
      overlayNum = varargin{1};
    else
      overlayNum = curOverlay;
    end
    analysisNum = viewGet(view,'currentAnalysis');
    if ~isempty(analysisNum) & ~isempty(overlayNum) & ...
        ~isempty(view.analyses{analysisNum}.overlays)
      view.analyses{analysisNum}.overlays(overlayNum).alpha = val;
      if (overlayNum == curOverlay)
        mlrGuiSet(view,'alpha',val);
      end
    end



    % -------------------------------------------
    % ROI

  case {'showrois'}
    % view = viewSet(view,'showrois',string);
    switch val
      case{'all','selected','all perimeter','selected perimeter','hide'}
        view.showROIs = val;
    end

  case {'newroi'}
    % view = viewSet(view,'newROI',roiStructure);
    % val must be a structure with the following required fields
    % - name: string
    % - viewType:
    % - color: color string
    % - coords: 4xN matrix (homogeneous coordinates)
    % - xform: 4x4 matrix
    % - voxelSize: 3 vector
    % - date: specifies when it was created or last modified

    % Check that is has the required fields
    [check ROI] = isroi(val);
    if ~check
      mrErrorDlg('Invalid ROI');
    end
    % Add it to view.ROIs
    pos = length(view.ROIs)+1;
    if (pos == 1)
      % First ROI
      view.ROIs = ROI;
    else
      view.ROIs(pos) = ROI;
    end
    % clear the cache of any old reference to this roi
    view = viewSet(view,'ROICache','clear',ROI.name);
    % Update the gui
    stringList = {view.ROIs(:).name};
    mlrGuiSet(view,'roiPopup',stringList);

  case {'deleteroi'}
    % view = viewSet(view,'deleteroi',roiNum);
    roinum = val;
    numrois = viewGet(view,'numberofrois');
    curroi = viewGet(view,'currentROI');
    % Remove it and reset currentROI
    view.ROIs = view.ROIs(roinum ~= [1:numrois]);
    if (numrois > 1)
      if (curroi > roinum)
        view = viewSet(view,'currentROI',curroi-1);
      elseif (roinum == curroi)
        view = viewSet(view,'currentROI',1);
      end
    else
      view.curROI = [];
    end
    % Update the gui
    stringList = {view.ROIs(:).name};
    if isempty(stringList)
      stringList = {'none'};
    end
    mlrGuiSet(view,'roiPopup',stringList);

  case {'currentroi','curroi','selectedroi'}
    % view = viewSet(view,'currentROI',roiNum);
    roiNum = val;
    numROIs = viewGet(view,'numberofROIs');
    if (roiNum > 0) & (roiNum <= numROIs)
      view.curROI = roiNum;
      % update popup menu
      mlrGuiSet(view,'roi',roiNum);
    end

  case {'roicoords'}
    % view = viewSet(view,'roiCoords',array,[roiNum]);
    if ~isempty(varargin)
      ROInum = varargin{1};
    else
      ROInum = viewGet(view,'currentROI');
    end
    if ~isempty(ROInum)
      view.ROIs(ROInum).coords = val;
      view.ROIs(ROInum).date = datestr(now);
    end

  case {'prevroicoords'}
    % view = viewSet(view,'prevROICoords',array);
    view.prevROIcoords = val;

  case 'roicolor'
    % view = viewSet(view,'roiColor',color,[roiNum]);
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).color = val;
    end

  case 'roiname'
    % view = viewSet(view,'roiName',nameString,[roiNum]);
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).name = val;
    end
    mlrGuiSet(view,'roipopup',{view.ROIs(:).name});

    % -------------------------------------------
    % Figure and GUI
  case 'roicache'
    % view = viewSet(view,'ROICache',roidata,[roiNum]);
    % view = viewSet(view,'ROICache','clear',roiName);
    if isstr(val)
      if strcmp(val,'clear')
        % clear the cache
        MLR.caches{view.viewNum}.roiCache = ...
          mrCache('clear',MLR.caches{view.viewNum}.roiCache,varargin{1});
      else
        % init the cache
        MLR.caches{view.viewNum}.roiCache = mrCache('init',mrGetPref('roiCacheSize'));
      end
      % add to the cache
    else
      if ~isempty(varargin)
	% with another argument, then we set the roi cache for
	% a specific ROI
	roiID = viewGet(view,'ROICacheID',varargin{1});
      else
	roiID = viewGet(view,'ROICacheID');
      end
      MLR.caches{view.viewNum}.roiCache = ...
        mrCache('add',MLR.caches{view.viewNum}.roiCache,roiID,val);
    end
  case 'basecache'
    % view = viewSet(view,'baseCache',basedata);
    % view = viewSet(view,'baseCache','clear',baseName);
    if isstr(val)
      if strcmp(val,'clear')
        % clear the cache
        MLR.caches{view.viewNum}.baseCache = ...
          mrCache('clear',MLR.caches{view.viewNum}.baseCache,varargin{1});
      else
        % init the cache
        MLR.caches{view.viewNum}.baseCache = mrCache('init',mrGetPref('baseCacheSize'));
      end
      % add to the cache
    else
      baseID = viewGet(view,'baseCacheID');
      MLR.caches{view.viewNum}.baseCache = ...
        mrCache('add',MLR.caches{view.viewNum}.baseCache,baseID,val);
    end
  case 'overlaycache'
    % view = viewSet(view,'overlayCache',overlaydata);
    % view = viewSet(view,'overlayCache','clear',overlayName);
    if isstr(val)
      if strcmp(val,'clear')
        % clear the cache
        MLR.caches{view.viewNum}.overlayCache = ...
          mrCache('clear',MLR.caches{view.viewNum}.overlayCache,varargin{1});
      else
        % init the cache
        MLR.caches{view.viewNum}.overlayCache = mrCache('init',mrGetPref('overlayCacheSize'));
      end
      % add to the cache
    else
      overlayID = viewGet(view,'overlayCacheID');
      MLR.caches{view.viewNum}.overlayCache = ...
        mrCache('add',MLR.caches{view.viewNum}.overlayCache,overlayID,val);
    end
  case {'figure'}
    % view = viewSet(view,'figure',handle);
    view.figure = val;

  case {'curslicebasecoords'}
    % view = viewSet(view,'curslicebasecoords',array);
    view.curslice.baseCoords = val;

  case {'cursliceoverlaycoords'}
    % view = viewSet(view,'cursliceoverlaycoords',array);
    view.curslice.overlayCoords = val;

  case {'sliceorientation','baseSliceIndex'}
    % view = viewSet(view,'sliceOrientation',n);
    if ~isscalar(val)
      switch val
        case {'sagittal'}
          sliceOrientation = 3;
        case {'coronal'}
          sliceOrientation = 2;
        case {'axial'}
          sliceOrientation = 1;
      end
    else
      sliceOrientation = val;
    end
    mlrGuiSet(view,'sliceOrientation',sliceOrientation);
    % Update slice and nSlices
    baseNum = viewGet(view,'currentBase');
    if ~isempty(baseNum)
      baseDims = viewGet(view,'basedims',baseNum);
      sliceIndex = viewGet(viewNum,'baseSliceIndex',baseNum);
      nSlices = baseDims(sliceIndex);
      coords = viewGet(view,'curCoords');
      slice = coords(sliceOrientation);
      mlrGuiSet(view,'nSlices',nSlices);
      mlrGuiSet(view,'slice',max(1,min(slice,nSlices)));
    end

  otherwise
    % viewType dependent
    switch(view.viewType)
      case 'Volume'
        view = volumeSet(view,param,val,varargin{:});
      case 'Surface'
        view = surfaceSet(view,param,val,varargin{:});
      case 'Flat'
        view = flatSet(view,param,val,varargin{:});
    end
end

% Side effect the global variable
MLR.views{viewNum} = view;

return;


%------------------------------
function view = volumeSet(view,param,val,varargin)

switch lower(param)
  case {'foo'}
    view.foo = val;
  otherwise
    error(['Invalid parameter for inplane view: ',param]);
end
return;

%------------------------------
function view = surfaceSet(view,param,val,varargin)

switch lower(param)
  case {'foo'}
    view.foo = val;
  otherwise
    error(['Invalid parameter for gray view: ',param]);
end
return;

%------------------------------
function view = flatSet(view,param,val,varargin)

switch lower(param)
  case {'foo'}
    view.foo = val;
  otherwise
    error(['Invalid parameter for flat view: ',param]);
end
return;

%----------------------------------------------------------

function view = reconcileAnalyses(view)
% This function is called by viewSet newscan and by viewSet deletescan to
% update the analyses and overlays, reflecting the change in the number of
% scans.

mrGlobals;

curgroup = viewGet(view,'currentGroup');
curGroupName = viewGet(view,'groupName',curgroup);

% Reconcile analysis params and overlay data/params with tseries
% files, adding empty data and default params to new scan.
for a = 1:viewGet(view,'numberofAnalyses')
  analysisReconcileFunction = viewGet(view,'reconcilefunction',a);
  analysisParams = viewGet(view,'analysisParams',a);
  view.analyses{a}.params = feval(analysisReconcileFunction,curGroupName,analysisParams);
  for ov = 1:viewGet(view,'numberofOverlays',a);
    overlay = viewGet(view,'overlay',ov,a);
    if ~isempty(overlay)
      overlayReconcileFunction = viewGet(view,'overlayReconcileFunction',ov,a);
      overlayParams = viewGet(view,'overlayParams',ov,a);
      overlayData = viewGet(view,'overlaydata',[],ov,a);
      [params,data] = feval(overlayReconcileFunction,curGroupName,overlayParams,overlayData);
      view.analyses{a}.overlays(ov).params = params;
      view.analyses{a}.overlays(ov).data = data;
    end
  end
end

% Do the same for all the other views
for v = find(view.viewNum ~= [1:length(MLR.views)])
  if isview(MLR.views{v})
    group = viewGet(v,'currentGroup');
    if (group == curgroup)
      for a = 1:viewGet(v,'numberofAnalyses')
	% jg-changed the viewGet(view... to viewGet(v... in
	% the lines below since that seems like what is wanted?
        analysisReconcileFunction = viewGet(v,'reconcilefunction',a);
        analysisParams = viewGet(v,'analysisParams',a);
        MLR.views{v}.analyses{a}.params = feval(analysisReconcileFunction,...
          curGroupName,analysisParams);
	% view->v
        for ov = 1:viewGet(v,'numberofOverlays',a);
          overlay = viewGet(v,'overlay',ov,a);
          if ~isempty(overlay)
	    % view->v
            overlayReconcileFunction = viewGet(v,'overlayReconcileFunction',ov,a);
            overlayParams = viewGet(v,'overlayParams',ov,a);
            overlayData = viewGet(v,'overlaydata',[],ov,a);
            [params,data] = feval(overlayReconcileFunction,...
              curGroupName,overlayParams,overlayData);
            MLR.views{v}.analyses{a}.overlays(ov).params = params;
            MLR.views{v}.analyses{a}.overlays(ov).data = data;
          end
        end
      end
    end
  end
end


%------------------------------
function dispViewSetHelp()

% open this file and scan the text
fid = fopen(which('viewSet'));
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
