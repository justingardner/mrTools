function [view tf] = viewSet(view,param,val,varargin)
%
%   view = viewSet(view,param,val,varargin)
%
% view can be either a view structure or a viewNum which is interpreted as
% MLR.views{viewNum}. Modifies the global variable MLR.views{viewNum} as
% well as the local and returned copy.
%
% ====================================================================
% Type viewSet w/out any arguments for the full set of parameters and for
% optional arguments. If you want help on a specific parameter, you 
% can also do:
% viewSet([],'parameterName','help')
% ====================================================================
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
%
% 6/2004 djh

mrGlobals

if ieNotDefined('param')
  dispViewSetHelp;
  return
elseif ~ieNotDefined('val') && ieNotDefined('varargin') && isstr(val) && (strcmp(lower(val),'help') || strcmp(lower(val),'?'))
  dispViewSetHelp(param);
  return
end

% the return value tf signifies whether the viewSet was successful
tf = 1;

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

  case{'currentscan','curscan'}
    % view = viewSet(view,'currentScan',n);
    nScans = viewGet(view,'nScans');
    if ((val > 0) && (val <= nScans))
      if viewGet(view,'curScan') ~= val
	view.curScan = val;
	mlrGuiSet(view,'scan',val);
      end
    else
      disp(sprintf('(viewSet) Scan %i out of range [%i:%i]',val,1,nScans));
    end
  case{'currentgroup','curgroup'}
    % view = viewSet(view,'currentGroup',n);
    % view = viewSet(view,'currentGroup','groupName')
    if isstr(val)
      % convert to a groupNum if we got passed in a groupname
      groupNum = viewGet(view,'groupNum',val);
      if isempty(groupNum)
	mrWarnDlg(sprintf('(viewSet:curGroup) No group %s found',val));
      end
      val = groupNum;
    end
    if isempty(val) || (val < 0) || (val > viewGet(view,'nGroups'))
      mrErrorDlg(sprintf('(viewSet) groupNum %i out of range: [1 %i]',val,viewGet(view,'nGroups')));
      return
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
      scanNum = max(min(viewGet(view,'groupScanNum',view.curGroup),nScans),1);
      mlrGuiSet(view,'scan',scanNum);
      view.curScan = scanNum;
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
 case{'renamegroup'}
    % rename a group. This will change the current group name (and change the directory name)
    % view = viewSet(view,'renameGroup',string);
    % Also, you can specify a groupnumber
    % view = viewSet(view,'renameGroup',string, groupNum);
    if length(varargin) == 0
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    % make sure groupNum is a number not a name
    groupNum = viewGet(view,'groupNum',groupNum);
    if isempty(groupNum)
      mrWarnDlg(sprintf('(viewSet:renameGroup) Could not find group %s',varargin{1}));
      return
    end
    % get old name
    oldGroupName = viewGet(view,'groupName',groupNum);
    oldGroupDirname = viewGet(view,'datadir',groupNum);
    % change to new name
    newGroupDirname = fullfile(viewGet(view,'homeDir'),val);
    % and move directory
    success = movefile(oldGroupDirname,newGroupDirname);
    % only change the group name if successful
    if success
      MLR.groups(groupNum).name = val;
      % Update the GUIs of all views
      stringList = {MLR.groups(:).name};
      for v = 1:length(MLR.views)
	if isview(MLR.views{v})
	  mlrGuiSet(MLR.views{v},'groupPopup',stringList);
	end
      end
    end
  case{'groupname'}
    % view = viewSet(view,'groupName',string);
    n = viewGet(view,'groupnum',val);
    view = viewSet(view,'currentGroup',n);

  case {'newgroup'}
    % view = viewSet(view,'newGroup',groupName);
    newgroup.name = val;
    % check to see if we already have the group
    if any(strcmp(newgroup.name,viewGet(view,'groupNames')))
      mrWarnDlg(sprintf('(viewSet:newGroup) Group %s already exists',newgroup.name));
      return
    end
    newgroup.scanParams = [];
    newgroup.auxParams = [];
    [tf newgroup] = isgroup(newgroup);
    % make sure all MLR.groups are valid groups
    if ~isempty(MLR.groups)
      for i = 1:length(MLR.groups)
	[tf oldgroups(i)] = isgroup(MLR.groups(i));
      end
      MLR.groups = oldgroups;
    end
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
    % make directories
    if ~isdir(fullfile(MLR.homeDir,val))
      mkdir(fullfile(MLR.homeDir,val));
    end
    if ~isdir(fullfile(MLR.homeDir,val,'TSeries'))
      mkdir(fullfile(MLR.homeDir,val,'TSeries'));
    end
    % Save mrSession
    saveSession;

  case {'deletegroup'}
    % view = viewSet(view,'deleteGroup',groupNum);
    groupnum = val;
    groupName = viewGet(view,'groupName',groupnum);
    nScans = viewGet(view,'nScans');
    % confirm with user
    if nScans > 0
      queststr = sprintf('There are %i scans in group %s. Are you sure you want to delete?',nScans,groupName);
    else
      queststr = sprintf('Are you sure you want to delete empty group: %s?',groupName);
    end
    if ~strcmp(questdlg(queststr,'Delete group'),'Yes')
      return
    end
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
    % Save mrSession
    saveSession;

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
  case {'scanvol2tal'}
    % view = viewSet(view,'scanvol2tal',vol2tal,[scanNum],[groupNum]);
    % set the vol2tal transform for the specified scan. 
    % This specifies the transformation from volume coordinates to talairach
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    % check that the passed in xform is 4x4
    if ~isequal(size(val),[4 4])
      mrErrorDlg('(viewSet) scanVol2tal transform should be a 4x4 matrix');
    end
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).vol2tal = val;
    end
  case {'scanvol2mag'}
    % view = viewSet(view,'scanvol2mag',vol2mag,[scanNum],[groupNum]);
    % set the vol2mag transform for the specified scan. 
    % This specifies the transformation from volume coordinates to magnet
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    % check that the passed in xform is 4x4
    if ~isequal(size(val),[4 4])
      mrErrorDlg('(viewSet) scanVol2mag transform should be a 4x4 matrix');
    end
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).vol2mag = val;
    end
  case {'scanxform'}
    % view = viewSet(view,'scanXform',sform,scanNum,groupNum);
    % This will set the scan sform to the passed in sform
    % and set the sform_code to 1
    % viewSet of scanXform is redundant with scanSform
    % and the latter is the preferred method because
    % it does not assume that the sform should be
    % in magnet coordinates.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr.sform44 = val;
      MLR.groups(g).scanParams(s).niftiHdr.sform_code = 1;
    end
  case {'base2scan'}
    % view = viewSet(view,'base2scan',base2scan)
    % This will set the base sform such that the base2scan
    % will come out as the input matrix. This is a very
    % specialized setting that should only really be used
    % as a last resort--usually mrAlign will give the right
    % transforms. But this allows for a quick way to get things
    % right. Will set the sform_code to 1.
    scan2mag = viewGet(view,'scan2mag');
    view = viewSet(view,'baseSform',scan2mag*val*inv(shiftOriginXform));
    view = viewSet(view,'baseSformCode',1);
  case {'scan2base'}
    % view = viewSet(view,'scan2base',scan2base)
    % This will set the scan sform such that the scan2base
    % will come out as the input matrix. This is a very
    % specialized setting that should only really be used
    % as a last resort--usually mrAlign will give the right
    % transforms. But this allows for a quick way to get things
    % right. Will set the sform_code to 1.
    base2mag = viewGet(view,'base2mag');
    view = viewSet(view,'scanSform',base2mag*val*inv(shiftOriginXform));
    view = viewSet(view,'scanSformCode',1);
  case {'scansform','sform'}
    % view = viewSet(view,'scanSform',sform,scanNum,groupNum);
    % The scanSform is the sform set in the nifti header.
    % If the sform_code is set to 1, then this field has
    % been set by mrAlign to be the tansformation from this
    % scan to the volume in magnet coordinates, 
    % If the sform_code is set to 3, then this field has
    % been set by mrAlign to be the tansformation from this
    % scan to the volume in talairach coordinates, 
    % You will need to viewSet the sform_code for the scan
    % accordingly if you change the scanSform.
    % viewSet of scanXform assumes that the sform is
    % in magnet coordinates.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr.sform44 = val;
      % now get the pix dimensions
      [q r] = qr(val);
      pixdim = diag(r(1:3,1:3));
      disp(sprintf('(viewSet:scanSform) The implied pixel dims are being reset to [%0.2f %0.2f %0.2f]',pixdim(1),pixdim(2),pixdim(3)));
      MLR.groups(g).scanParams(s).niftiHdr.pixdim(2:4) = pixdim;
      MLR.groups(g).scanParams(s).voxelSize = pixdim';
    end
  case {'sformcode','scansformcode'}
    % view = viewSet(view,'sformcode',sformcode,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).niftiHdr.sform_code = val;
    end
  case {'scanvoxelsize'}
    % view = viewSet(view,'scanVoxelSize',voxelSize,scanNum,groupNum);
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      MLR.groups(g).scanParams(s).voxelSize = val;
      MLR.groups(g).scanParams(s).niftiHdr.pixdim(2:4) = val;
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
    % scanParams.framePeriod = hdr.pixdim(5)/1000;
    % see definitions for NIFTI_UNITS in nifti1.h
    % space units are multiples of 1, time units multiples of 8 up to 64
    niftiSpaceUnit = rem(hdr.xyzt_units, 8); 
    niftiTimeUnit = rem(hdr.xyzt_units-niftiSpaceUnit, 64);
    if niftiTimeUnit == 8 % seconds
      scanParams.framePeriod = hdr.pixdim(5)./1;
    elseif niftiTimeUnit == 16 % milliseconds
      scanParams.framePeriod = hdr.pixdim(5)./1000;
    elseif niftiTimeUnit == 32 % microseconds
      scanParams.framePeriod = hdr.pixdim(5)./10e6;
    else % for old Analyze files
      scanParams.framePeriod = hdr.pixdim(5)./1000;
    end
    disp(sprintf('(viewSet) framePeriod set to: %f',scanParams.framePeriod))
    if strcmp(lower(mrGetPref('verbose')),'yes')
      % 8 -> 10^0, 16 -> 10^3, 32-> 10^6
      disp(sprintf('(viewSet) Timing. Pixdim(5) units: %d. Scaling by 10e%d',niftiTimeUnit, 3*(log2(niftiTimeUnit)-3)));
    end
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
    numBaseVolumes = viewGet(view,'numberofBaseVolumes'); %get the number of loaded bases
    if numBaseVolumes %if there are bases loaded, check for base with identical name
       % If baseVolume.name already exists then replace the existing one
       % with this new one. Otherwise, add it to the end of the baseVolumes
       % list.
       newBaseName = baseAnatomy.name;
       newBaseNum = [];
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
         view.baseVolumes(newBaseNum) = baseAnatomy;
    else    %if there is no base loaded, add the base in first position
       newBaseNum = 1;
       view.baseVolumes = baseAnatomy;
    end
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
    % Set it to be the current base Volume
    view = viewSet(view,'curBase',newBaseNum);

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

  case{'rotate'}
    % view = viewSet(view,'rotate',rotation);
    curBase = viewGet(view,'curBase');
    numBases = viewGet(view,'numberofBaseVolumes');
    baseType = viewGet(view,'baseType');
    if (curBase > 0) & (curBase <= numBases)
      % surfaces are rotated differently, the
      % rotate field causes surfaceRotate to
      % change rather than rotate. So we
      % need to set the appropriate field
      if baseType <= 1
	view.baseVolumes(curBase).rotate = val;
      else
	view.baseVolumes(curBase).coordMap.rotate = val;
      end
    end
      
  case{'tilt'}
    % view = viewSet(view,'tilt',tilt);
    curBase = viewGet(view,'curBase');
    numBases = viewGet(view,'numberofBaseVolumes');
    baseType = viewGet(view,'baseType');
    if (curBase > 0) & (curBase <= numBases)
      % surfaces are the ones that can be tilted.
      if baseType == 2
	view.baseVolumes(curBase).tilt = val;
	mlrGuiSet(view,'baseTilt',val);
      end
    end
      
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
      % set gamma sliders
      baseClip = viewGet(view,'baseClip',baseNum);
      baseRange = viewGet(view,'baseRange',baseNum);
      baseType = viewGet(view,'baseType',baseNum);
      baseGamma = viewGet(view,'baseGamma',baseNum);
      mlrGuiSet(view,'baseType',baseType);
      mlrGuiSet(view,'baseGamma',baseGamma);
      % get the last know settings for this base
      baseCurSlice = viewGet(view,'baseCurSlice',baseNum);
      baseSliceOrientation = viewGet(view,'baseSliceOrientation',baseNum);
      baseRotate = viewGet(view,'baseRotate',baseNum);
      % if this is a flat then we need to set the sliceIndex to 3
      % since we shouldn't be able to view in any other
      % configuration
      if baseType
	baseSliceOrientation = 1;
      end
      % set the slice orientation if there is a valid one saved
      if ~isempty(baseSliceOrientation)
	view = viewSet(view,'sliceOrientation',baseSliceOrientation);
      end
      % update nSlices and reset slice to be within range
      sliceIndex = viewGet(view,'basesliceindex',baseNum);
      nSlices = baseDims(sliceIndex);
      curSlice = viewGet(view,'curSlice');
      % if the base has a current slice set, then use that
      if isempty(baseCurSlice) || (baseCurSlice > nSlices)
	view = viewSet(view,'curSlice',max(1,min(curSlice,nSlices)));
      else
	view = viewSet(view,'curSlice',min(baseCurSlice,nSlices));
      end
      if ~isempty(baseRotate)
	mlrGuiSet(view,'rotate',baseRotate);
      end
      mlrGuiSet(view,'nSlices',nSlices);
      baseTilt = viewGet(view,'baseTilt',baseNum);
      if baseType == 2
	mlrGuiSet(view,'baseTilt',baseTilt);
	if ~mrInterrogator('isactive',viewGet(view,'viewNum'));
	  % turn on free rotation
	  mlrSetRotate3d(view,1);
	else
	  mlrSetRotate3d(view,0);
	end
      else
	% otherwise turn off free rotation
	mlrSetRotate3d(view,0);
      end
    end

  case{'basecoordmappath'}
    % view = viewSet(view,'basecoordmapdir',baseCoordMapPath,[baseNum]);
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum)
      baseType = viewGet(view,'baseType',baseNum);
      if isfield(view.baseVolumes(baseNum),'coordMap')
	view.baseVolumes(baseNum).coordMap.path = val;
      end
    end
  case{'basecoordmap'}
    % view = viewSet(view,'basecoordmap',baseCoordMap,[baseNum]);
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum)
      if isfield(view.baseVolumes(baseNum),'coordMap')
	view.baseVolumes(baseNum).coordMap = val;
      end
    end

  case{'basemin'}
    % view = viewSet(view,'basemin',number,[baseNum]);
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum)
      view.baseVolumes(baseNum).clip(1) = val;
    end

  case{'basemax'}
    % view = viewSet(view,'basemax',number,[baseNum]);
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum)
      view.baseVolumes(baseNum).clip(2) = val;
    end

  case{'baserange'}
    % view = viewSet(view,'baserange',[min max],[baseNum]);
    range = val;
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).range = range;
    end

  case{'basegamma'}
    % view = viewSet(view,'baseGamma',gamma,[baseNum]);
    gamma = val;
    curBase = viewGet(view, 'curBase');
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).gamma = gamma;
    end
    if (baseNum == curBase)
      mlrGuiSet(view,'baseGamma',gamma);
    end

  case {'basevol2mag'}
    % view = viewSet(view,'baseVol2mag',xform,[baseNum]);
    % This specifies the transformation from volume coordinates to magnet
    % coordinates of the base volume that this base was aligned to.
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).vol2mag = val;
    end
  case {'basevol2tal'}
    % view = viewSet(view,'baseVol2tal',xform,[baseNum]);
    % This specifies the transformation from volume coordinates to talairach
    % coordinates of the base volume that this base was aligned to.
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).vol2tal = val;
    end
  case {'basetalinfo'}
    % view = viewSet(view,'basetalinfo',talinfo,[baseNum]);
    % This specifies the information necessary to define the talairach transformation
    % in terms of set points (AC, PC, AAC, etc)
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      view.baseVolumes(baseNum).talInfo = val;
    end
  case {'basexform'}
    % view = viewSet(view,'baseXform',sform,[baseNum]);
    % set the base sform to the passed in value and
    % sets the base sform_code to 1
    % viewSet of baseXform is redundant with baseSform
    % and the latter is the preferred method because
    % it does not assume that the sform should be
    % in magnet coordinates.
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      % get base header
      hdr = view.baseVolumes(baseNum).hdr;
      % set sform
      hdr = cbiSetNiftiSform(hdr,val);
      hdr.sform_code = 1;
      % reset in structure
      view.baseVolumes(baseNum).hdr = hdr;
    end
  case {'basesform'}
    % view = viewSet(view,'basesform',sform,[baseNum]);
    % This returns the sform of the base, which should
    % be set by mrAlign.
    % This is the transformation from the base to
    % volume in magnet coordinates when the sform_code
    % of the base is set to 1.
    % This is the transformation from the base to
    % volume in talairach coordinates when the sform_code
    % of the base is set to 3.
    % You will need to viewSet the sform_code for the base
    % accordingly if you change the baseSform.
    % viewSet of baseXform assumes that the sform is
    % in magnet coordinates.
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes)
      % get base header
      hdr = view.baseVolumes(baseNum).hdr;
      % set sform, but don't necessarily set sform_code
      sform_code = hdr.sform_code;
      hdr = cbiSetNiftiSform(hdr,val);
      hdr.sform_code = sform_code;
      % reset in structure
      view.baseVolumes(baseNum).hdr = hdr;
      % now get the pix dimensions
      [q r] = qr(val);
      pixdim = diag(r(1:3,1:3));
      disp(sprintf('(viewSet:baseSform) The implied pixel dims are being reset to [%0.2f %0.2f %0.2f]',pixdim(1),pixdim(2),pixdim(3)));
      view.baseVolumes(baseNum).hdr.pixdim(2:4) = pixdim;
    end

  case {'basesformcode'}
    % view = viewSet(view,'baseSformCode',sformCode,[baseNum]);
    baseNum = getBaseNum(view,varargin);
    if ~isempty(baseNum) & ~isempty(view.baseVolumes) 
      view.baseVolumes(baseNum).hdr.sform_code = val;
    end
  case{'basename'}
    % view = viewSet(view,'basename',nameString,[baseNum]);
    baseNum = getBaseNum(view,varargin);
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
      mrWarnDlg('(viewSet:newAnalysis) Invalid analysis');
    end
    % Error if groupNames don't match
    curGroupNum = viewGet(view,'currentGroup');
    curGroupName = viewGet(view,'groupName',curGroupNum);
    if ~strcmp(analysis.groupName,curGroupName)
      mrWarnDlg(['(viewSet:newAnalysis) Analysis is incompatible with group: ',curGroupName]);
    end
    % Reconcile analysis params with group
    if isfield(analysis,'d') && ~isempty(analysis.d)
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
    if isempty(curAnalysis),return,end
    % Remove it and reset currentAnalysis
   view.analyses = {view.analyses{setdiff(1:numAnalyses,analysisNum)}}; %JB: allows deleting several analyses at once
   numAnalyses = viewGet(view,'numberofAnalyses');                      
   if (curAnalysis > numAnalyses)                                       % JB: does not change current analysis, unless not enough analyses         
     view = viewSet(view,'currentAnalysis',numAnalyses);                % in which case set the last analysis
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

  case{'analysisinterrogator'}
    % view = viewSet(view,'analysisInterrogator',interrogatorName,[analysisNum]);
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if ~isempty(analysisNum)
      for i = 1:viewGet(view,'nOverlays')
	view.analyses{analysisNum}.overlays(i).interrogator = val;
      end
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
      mrErrorDlg(['(viewSet:newOverlay) overlay is incompatible with group: ',groupName]);
    end
    % Reconcile overlay data and params with tseries files, reordering
    % them as necessary.
    [params,data] = feval(overlay.reconcileFunction,groupName,overlay.params,overlay.data);
    overlay.params = params;
    overlay.data = [];
    for scanNum = 1:length(data)
      % to save memory, we keep overlays as single precision
      % (commented out now)
      %overlay.data{scanNum} = single(data{scanNum});
      overlay.data{scanNum} = data{scanNum};
    end
    % Check that it has the correct number of scans
    if (length(overlay.data) ~= nScans)
      %mrErrorDlg('Invalid overlay, incorrect number of scans.');
      if isempty(data)
	disp('(viewSet:newOverlay) No data in overlay');
      else
	disp('(viewSet:newOverlay) Invalid overlay, incorrect number of scans.');
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
        mrErrorDlg('(viewSet:newOverlay) Invalid overlay, incorrect data array size.');
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
    % clear overlay cache
    view = viewSet(view,'overlayCache','clear',sprintf('_%i_%i_',analysisNum,newOverlayNum));
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
         view.analyses{analysisNum}.overlays(setdiff(1:numoverlays,overlayNum) ); %JB: allows deleting several overlays at once
      numoverlays = viewGet(view,'numberofoverlays',analysisNum);                %
      if (curoverlay > numoverlays)                                              % JB: does not change current overlay, unless not enough overlays         
         view = viewSet(view,'currentOverlay',numoverlays);                      % in which case set the last overlay
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
    overlayNum = viewGet(view,'overlayNum',val);
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

  case {'overlayname'}
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

  case {'overlaycmap'}
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

  case {'overlaymin'}
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

  case {'overlaymax'}
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

  case {'overlayrange'}
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

  case {'alpha'}
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
      case{'all','selected','all perimeter','selected perimeter','group','group perimeter','hide'}
        view.showROIs = val;
    end
    mlrGuiSet(view,'showROIs',val);
  case {'roigroup'}
    % view = viewSet(view,'showrois',string);
    % set the cell array of roi names that will be displayed
    % when ROI view is set to group or group perimeter
    view.roiGroup = val;
  case {'labelrois'}
    % view = viewSet(view,'labelrois',val);
    view.labelROIs = val;
    mlrGuiSet(view,'labelROIs',val);
  case {'newroi'}
    % view = viewSet(view,'newROI',roiStructure,<replaceDuplicate>);
    % val must be a structure with the following required fields
    % - name: string
    % - viewType:
    % - color: color string
    % - coords: 4xN matrix (homogeneous coordinates)
    % - xform: 4x4 matrix
    % - voxelSize: 3 vector
    % - date: specifies when it was created or last modified
    % if replaceDuplicate=1 then if the roi already exists
    % than it will replace the new one.

    % Check that is has the required fields
    [check ROI] = isroi(val);
    if ~check
      mrErrorDlg('Invalid ROI');
    end
    if length(varargin) > 0
      replaceDuplicates = varargin{1};
    else
      replaceDuplicates = 0;
    end
    % check to make sure the name is unique
    roiNames = viewGet(view,'roiNames');
    nameMatch = find(strcmp(ROI.name,roiNames));
    while ~isempty(nameMatch)
      if ~replaceDuplicates
	paramsInfo{1} = {'roiName',ROI.name,'Change the name to a unique ROI name'};
	paramsInfo{2} = {'replace',0,'type=checkbox','Check to replace the loaded ROI with the same name'};
	params = mrParamsDialog(paramsInfo,'Non unique ROI name, please change');
	if isempty(params),tf=0;return,end
      else
	params.replace = replaceDuplicates;
      end
      if params.replace
	% if replace, then delete the existing one
	view = viewSet(view,'deleteROI',viewGet(view,'roiNum',ROI.name));
	roiNames = viewGet(view,'roiNames');
      else
	% change the name
	ROI.name = params.roiName;
      end
      nameMatch = find(strcmp(ROI.name,roiNames));
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
    % add the roi to the current roigroup
    roiGroupNames = viewGet(view,'roiGroupNames');
    roiGroupNames{end+1} = ROI.name;
    view = viewSet(view,'roiGroup',roiGroupNames);
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

  case {'roicolor'}
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

  case {'roinotes'}
    % v = viewSet(v,'roiNotes',notes,[roiNum]);
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).notes = val;
    end

  case {'roivol2mag'}
    % v = viewSet(v,'roiVol2mag',xform,[roiNum]);
    % sets the xform of the volume coordinates to
    % the magnet coordinates for the canonical volume
    % for the called for roi
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).vol2mag = val;
    end

  case {'roivol2tal'}
    % v = viewSet(v,'roiVol2mag',xform,[roiNum]);
    % sets the xform of the volume coordinates to
    % the talairach coordinates for the canonical volume
    % for the called for roi
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).vol2tal = val;
    end
 
  case {'roisformcode'}
    % v = viewSet(v,'roisformcode',sformcode,[roiNum]);
    % sets the sformcode for the roiXform
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).sformCode = val;
    end

  case {'roixform'}
    % v = viewSet(v,'roixform',xform,[roiNum]);
    % sets the xform for the roi
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).xform = val;
    end
  case {'roivoxelsize'}
    % v = viewSet(v,'roiVoxelSize',voxelSize,[roiNum]);
    % sets the voxel size for the roi
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).voxelSize = val;
    end

  case {'roiname'}
    % view = viewSet(view,'roiName',nameString,[roiNum]);
    curRoi = viewGet(view,'currentRoi');
    if ~isempty(varargin)
      roiNum = varargin{1};
    else
      roiNum = curRoi;
    end
    % check to make sure the name is unique
    roiNames = viewGet(view,'roiNames');
    nameMatch = find(strcmp(val,roiNames));
    while ~isempty(nameMatch) && (nameMatch~=roiNum)
      paramsInfo{1} = {'roiName',val,'Change the name to a unique ROI name'};
      params = mrParamsDialog(paramsInfo,'Non unique ROI name, please change');
      if isempty(params),tf=0;return,end
      val = params.roiName;
      nameMatch = find(strcmp(val,roiNames));
    end
    if ~isempty(roiNum)
      view.ROIs(roiNum).name = val;
    end
    mlrGuiSet(view,'roipopup',{view.ROIs(:).name});

    % -------------------------------------------
    % Figure and GUI
  case {'roicache'}
    % view = viewSet(view,'ROICache',roidata,[roiNum]);
    % view = viewSet(view,'ROICache','clear',roiName);
    if isstr(val)
      if strcmp(val,'clear')
	if (length(MLR.caches) >= view.viewNum) && isfield(MLR.caches{view.viewNum},'roiCache')
	  % clear the cache
	  MLR.caches{view.viewNum}.roiCache = ...
	      mrCache('clear',MLR.caches{view.viewNum}.roiCache, varargin{1});
	end
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
  case {'basecache'}
    % view = viewSet(view,'baseCache',basedata);
    % view = viewSet(view,'baseCache','clear',baseName);
    if isstr(val)
      if strcmp(val,'clear')
        % clear the cache
	if (length(MLR.caches) >= view.viewNum)
	  MLR.caches{view.viewNum}.baseCache = ...
	      mrCache('clear',MLR.caches{view.viewNum}.baseCache,varargin{1});
	end
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
  case {'overlaycache'}
    % view = viewSet(view,'overlayCache',overlaydata);
    % view = viewSet(view,'overlayCache','clear',overlayName);
    if isstr(val)
      if strcmp(val,'clear')
        % clear the cache
	if (length(MLR.caches) >= view.viewNum)
	  MLR.caches{view.viewNum}.overlayCache = ...
	      mrCache('clear',MLR.caches{view.viewNum}.overlayCache,varargin{1});
	end
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

  case {'curslice'}
    % view = viewSet(view,'curSlice',sliceNum);
    baseDims = viewGet(view,'baseDims');
    sliceIndex = viewGet(view,'basesliceindex');
    if ~isempty(baseDims)
      nSlices = baseDims(sliceIndex);
    else
      nSlices = 0;
    end
    if (val > 0) && (val <= nSlices)
      if viewGet(view,'curSlice') ~= val
	view.curslice.sliceNum = val;
	mlrGuiSet(view,'slice',val);
      end
    else
      disp(sprintf('(viewSet) Slice %i out of range: [1 %i]',val,nSlices));
    end
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
        case 'sagittal'
          sliceOrientation = 3;
        case 'coronal'
          sliceOrientation = 2;
        case 'axial'
          sliceOrientation = 1;
      end
    else
      sliceOrientation = val;
    end
    if ((sliceOrientation > 0) && (sliceOrientation <= 3))
      % set slice orientation in view
      view.sliceOrientation = sliceOrientation;
      % Update slice and nSlices
      baseNum = viewGet(view,'currentBase');
      if ~isempty(baseNum) && ~isempty(viewGet(view,'fignum'))
	mlrGuiSet(view,'sliceOrientation',sliceOrientation);
	baseDims = viewGet(view,'basedims',baseNum);
	sliceIndex = viewGet(view,'baseSliceIndex',baseNum);
	nSlices = baseDims(sliceIndex);
	coords = viewGet(view,'curCoords');
	slice = coords(sliceOrientation);
	view = viewSet(view,'curSlice',min(slice,nSlices));
%	view = viewSet(view,'curSlice',slice);
	mlrGuiSet(view,'nSlices',nSlices);
	mlrGuiSet(view,'slice',max(1,min(slice,nSlices)));
      end
    end

 case {'defaultinterrogators'}
    % view = viewSet(view,'defaultInterrogators',defaultInterrogators,<replaceCurrentDefaultInterrogators>)
    % set a cell array of function names that are
    % defaultInterrogators that will show up in the interrogators list
    % this is similar to using the mrSetPref for
    % defaultInterrogators but is not presistent (i.e. will only
    % show up for the current session - it is used by mlrAdjustGUI
    % by default will add the interrogators to the current list,
    % set replaceCurrentDefaultInterrogators to false if you want
    % to replace the current list
    val = cellArray(val);
    % check that everything is an m-file
    defaultInterrogators = {};
    for i = 1:length(val)
      if exist(val{i}) ~= 2
	disp(sprintf('(viewSet:defaultInterrogators) %s is not an m-file on the path. Not adding to default interrogators list.',val{i}));
      else
	defaultInterrogators{end+1} = val{i};
      end
    end
    % if replaceCurrentDefaultInterrogators is set then
    % just save in MLR variables
    if ((length(varargin) > 0) && varargin{1}) || ~isfield(MLR,'defaultInterrogators')
      MLR.defaultInterrogators = defaultInterrogators;
    else
      %otherwise append to list
      MLR.defaultInterrogators = {MLR.defaultInterrogators{:} defaultInterrogators{:}};
    end
 case {'colormaps'}
    % view = viewSet(view,'colormaps',colormapList,<replaceCurrentColormaps>)
    % set a cell array of function names that are
    % colormaps that will show up in places like /edit/overlays
    % This is not presistent (i.e. will only
    % show up for the current session - it is used by mlrAdjustGUI)
    % By default will add the colormaps to the current list,
    % set replaceCurrentColormaps to false if you want
    % to replace the current list
    val = cellArray(val);
    % check that everything is an m-file
    colormaps = {};
    for i = 1:length(val)
      if exist(val{i}) ~= 2
	disp(sprintf('(viewSet:colormaps) %s is not an m-file on the path. Not adding to colormaps.',val{i}));
      else
	colormaps{end+1} = val{i};
      end
    end
    % if replaceCurrentColormaps is set then
    % just save in MLR variables
    if ((length(varargin) > 0) && varargin{1}) || ~isfield(MLR,'colormaps')
      MLR.colormaps = colormaps;
    else
      %otherwise append to list
      MLR.colormaps = {MLR.colormaps{:} colormaps{:}};
    end
 otherwise
   mrWarnDlg(sprintf('(viewSet) Unknown parameter %s',param));
end

% Side effect the global variable
MLR.views{viewNum} = view;

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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispViewSetHelp   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispViewSetHelp(param)

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

% if the user wants help on a particular command
if ~ieNotDefined('param')
  % find the command
  commandNum = [];
  % find the command
  for i = 1:length(commands)
    if ~isempty(findstr(lower(param),commands{i}))
      commandNum = i;
    end
  end
  % if found, print out the comments
  if ~isempty(commandNum) && (commandNum > 0) && (commandNum <= length(commandComments)) && ~isempty(commandComments{commandNum})
    disp(commandComments{commandNum});
  else
    disp(sprintf('(viewSet) No viewSet for %s',param));
  end
  return
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

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getScanAndGroup   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
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
% if group is a string, then convert it to a number
if isstr(g)
  g = viewGet(view,'groupNum',g);
end

%%%%%%%%%%%%%%%%%%%%
%%   getBaseNum   %%
%%%%%%%%%%%%%%%%%%%%
function baseNum = getBaseNum(view,varg)

curBase = viewGet(view,'currentBase');
if ieNotDefined('varg')
  baseNum = curBase;
else
  baseNum = varg{1};
end
