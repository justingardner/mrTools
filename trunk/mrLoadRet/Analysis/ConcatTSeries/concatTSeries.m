% concatTSeries.m
%
%        $Id$
%      usage: concatTSeries(view, params)
%         by: justin gardner
%       date: 10/12/06
%    purpose: concatenate time series together.
%
%             to just get a default parameter structure:
% 
%             v = newView;
%             [v params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1','scanList=[1 2]');
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = concatTSeries(v,[],'justGetParams=1');
%
function [view params] = concatTSeries(view,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help concatTSeries
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

interpTypes = {'nearest','linear','spline','cubic'};
if find(strcmp(mrGetPref('interpMethod'),interpTypes))
  interpTypes = putOnTopOfList(mrGetPref('interpMethod'),interpTypes);
end

% description of paramaters (used by mrParamsDialog functions)
paramsInfo = {...
    {'groupName',putOnTopOfList(viewGet(view,'groupName'),viewGet(view,'groupNames')),'Name of group from which to make concatenation'},...
    {'newGroupName','Concatenation','Name of group to which concatenation will be saved. If group does not already exist, it will be created.'},...
    {'description','Concatenation of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'filterType',1,'minmax=[0 1]','incdec=[-1 1]','Which filter to use, for now you can only turn off highpass filtering'},...
    {'filterCutoff',0.01,'minmax=[0 inf]','Highpass filter cutoff in Hz'},...
    {'percentSignal',1,'type=checkbox','Convert to percent signal change'},...
    {'warp',1,'type=checkbox','Warp images based on alignment. This can be used to concatenate together scans taken on different days. If the scans have the same slice prescription this will not do anything.'},...
    {'warpInterpMethod',interpTypes,'Interpolation method for warp','contingent=warp'}
	     };
% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo);
  end
  % no params means user hit cancel
  if isempty(params),return,end
  % select scans
  view = viewSet(view, 'groupName', params.groupName);
  if ~ieNotDefined('scanList')
    params.scanList = scanList;
  elseif defaultParams
    params.scanList = 1:viewGet(view,'nScans');
  else
    params.scanList = selectScans(view);
  end
  if isempty(params.scanList),return,end
  % if warp is set, then ask which scan to use as base scan for warp
  if params.warp
    % create a list of scans
    for i = 1:viewGet(view,'nScans')
      scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(view,'description',i),viewGet(view,'tSeriesFile',i));
    end
    warpParams = mrParamsDialog({{'warpBaseScan',scanNames,'The scan that will be used as the base scan to warp all the other scans to'}});
    if isempty(warpParams),return,end
    params.warpBaseScan = find(strcmp(warpParams.warpBaseScan,scanNames));
  end
  % check the parameters
  params = mrParamsReconcile(params.groupName,params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;
  params = mrParamsReconcile(params.groupName,params);
end
drawnow;

% Abort if params empty
if ieNotDefined('params'),return,end

% if just getting params then return
if justGetParams,return,end
  
% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',params.groupName);

if (groupNum == 0)
  mrErrorDlg('concatTSeries: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the average group name. Create the
% group if necessary.
viewConcat = newView(viewGet(view,'viewType'));
concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
if isempty(concatGroupNum)
  view = viewSet(view,'newgroup',params.newGroupName);
  concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
end
viewConcat = viewSet(viewConcat,'currentGroup',concatGroupNum);

% Check that all scans in scanList have the same 
% scanxform, scanvoxelsize, scandims
d.tr = viewGet(viewBase,'framePeriod',params.scanList(1));
d.voxelSize = viewGet(viewBase,'scanvoxelsize',params.scanList(1));
d.dim = viewGet(viewBase,'scandims',params.scanList(1));
d.nFrames = viewGet(viewBase,'nFrames',params.scanList(1));
d.dim(4) = d.nFrames;
d.sform = viewGet(viewBase,'scanSform',params.scanList(1));
if ~params.warp
  d.vol2mag = viewGet(viewBase,'scanVol2mag',params.scanList(1));
  d.vol2tal = viewGet(viewBase,'scanVol2tal',params.scanList(1));
else
  d.vol2mag = viewGet(viewBase,'scanVol2mag',params.warpBaseScan);
  d.vol2tal = viewGet(viewBase,'scanVol2tal',params.warpBaseScan);
end
for iscan = 1:length(params.scanList)
  % check frame period for mismatch
  if (viewGet(viewBase,'framePeriod',params.scanList(iscan)) ~= d.tr)
    mrWarnDlg(sprintf('concatTSeries: These scans have different TR. (%0.4f vs %0.4f)',viewGet(viewBase,'framePeriod',params.scanList(iscan)),d.tr));
  end
  % make sure that the voxel sizes are the same to within roundoff error
  baseVoxelSize = viewGet(viewBase,'scanvoxelsize',params.scanList(iscan));
  roundoff = 100000;
  if ~isequal(round(baseVoxelSize*roundoff)/roundoff,round(d.voxelSize*roundoff)/roundoff)
    disp(sprintf('(concatTSeries) Scans have different voxel sizes %i:[%s]~=[%s]',params.scanList(iscan),num2str(baseVoxelSize),num2str(d.voxelSize)));
  end
  % check the scan dims
  if ~isequal(viewGet(viewBase,'scandims',params.scanList(iscan)),d.dim(1:3))
    disp('(concatTSeries) Scans have different dimensions.');
  end
  % check the scan sforms
  if ~isequal(viewGet(viewBase,'scanSform',params.scanList(iscan)),d.sform)
    % this is only an issue if warp is  not set
    if ~params.warp
      mrWarnDlg(sprintf('(concatTSeries) Sform for scan %s:%i does not match %s:%i. This means that they have different slice prescriptions. Usually you should select warp in this case so that the different scans are all warped together to look like the base scan. You have not selected warp.',params.groupName,params.scanList(iscan),params.groupName,params.scanList(1)));
    end
  end
  % check the vol2mag and vol2tal
  if (~isequal(viewGet(viewBase,'scanVol2mag',params.scanList(iscan)),d.vol2mag) || ...
      ~isequal(viewGet(viewBase,'scanVol2tal',params.scanList(iscan)),d.vol2tal))
    % this is only an issue if warp is  not set
    if ~params.warp
      mrWarnDlg(sprintf('(concatTSeries) The scanVol2mag/scanVol2tal for scan %s:%i does not match %s:%i. This means that they have been aligned to different volume anatomies. Usually you should select warp in this case so that the different scans are all warped together to look like the base scan. You have not selected warp.',params.groupName,params.scanList(iscan),params.groupName,params.scanList(1)));
    end
  end
end
disp(sprintf('(concatTSeries) FramePeriod for scan is: %0.2f',d.tr));

% initialize some things
concatInfo.n = 0;
concatInfo.whichScan = [];
concatInfo.whichVolume = [];

set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
tic
% Compute output volume
waitHandle = mrWaitBar(0,'Concatenating tSeries.  Please wait...');
for iscan = 1:length(params.scanList)
  scanNum = params.scanList(iscan);
  
  % Load it
  mrDisp(sprintf('\nLoading scan %i from %s\n',scanNum,viewGet(viewBase,'groupName')));
  d.data = loadTSeries(viewBase,scanNum,'all');
	
  % Dump junk frames
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  d.dim(4) = nFrames;
  d.data = d.data(:,:,:,junkFrames+1:junkFrames+nFrames);

  % get the total junked frames. This is the number of frames
  % we have junked here, plus whatever has been junked in previous ones
  totalJunkedFrames(iscan) = junkFrames+viewGet(viewBase,'totalJunkedFrames',scanNum);
  
  % Compute transform
  if params.warp
    % get the scan2scan xform
    scan2scan = viewGet(view,'scan2scan',params.warpBaseScan,groupNum,scanNum,groupNum);
    
    if ~isequal(scan2scan,eye(4))
      % swapXY seems to be needed here, presumably becuase of the way that 
      % warpAffine3 works.
      swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];

      % compute transformation matrix
      M = swapXY * scan2scan * swapXY;

      % display transformation
      disp(sprintf('Transforming %s:%i to match %s:%i with transformation: ',params.groupName,scanNum,params.groupName,params.warpBaseScan));
      for rownum = 1:4
	disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
      end
      % Warp the frames
      for frame = 1:d.nFrames
	d.data(:,:,:,frame) = warpAffine3(d.data(:,:,:,frame),M,NaN,0,params.warpInterpMethod);
      end   
    end
  end
  
  % do other  processing here
  if params.filterType == 1
    d = eventRelatedHighpass(d,params.filterCutoff);
  elseif params.filterType == 2
    d = detrendTSeries(d)
  end
    
  % convert to percent signal change
  warning off                           % to avoid divide by zero warnings...
  if params.percentSignal
    d.mean = mean(d.data,4);
    % for means that are zero, divide by nan
    d.mean(d.mean==0) = nan;

    disppercent(-inf, 'Converting to percent signal change');
    for i = 1:d.dim(4)
      d.data(:,:,:,i) = (d.data(:,:,:,i)./d.mean);
      if params.percentSignal == 2           % scale it to mean of 1,000
          params.scaleFactor = 10000;
          d.data(:,:,:,i) = d.data(:,:,:,i) * params.scaleFactor;
      end
      disppercent(i/d.dim(4));
    end
    disppercent(inf);
  end
  warning on
  d.data(isnan(d.data))=0;              % b/c nan's can be annoying 
  
  % get the path and filename
  [path,filename,ext,versn] = fileparts(viewGet(viewBase,'tseriesPath',scanNum));
  baseGroupName = viewGet(viewBase,'groupName');

  % Save TSeries (using header of 1st scan on scanList as the template for
  % the nifti header), and add it as a new scan.
  if iscan == 1
    scanParams.fileName = [];
    scanParams.junkFrames = 0;
    scanParams.nFrames = d.nFrames;
    scanParams.description = params.description;
    scanParams.originalFileName{1} = filename;
    scanParams.originalGroupName{1} = baseGroupName;
    scanParams.totalJunkedFrames = totalJunkedFrames;
    scanParams.vol2mag = d.vol2mag;
    scanParams.vol2tal = d.vol2tal;

    hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',params.scanList(1)));
    % data *MUST* be written out as float32 b/c of the small values-epm
    hdr.datatype = 16;
    % if we are warping, then we need to change the sform to the
    % sform of the scan we warped to
    if params.warp
      hdr.sform44 = viewGet(view,'scanXform',params.warpBaseScan,groupNum);
      hdr.sform_code = viewGet(view,'scanSformCode',params.warpBaseScan,groupNum);
     end
    [viewConcat,tseriesFileName] = saveNewTSeries(viewConcat,d.data,scanParams,hdr);
    % get new scan number
    saveScanNum = viewGet(viewConcat,'nScans');
    
    % now load up channels
    stimfile = viewGet(viewBase,'stimFile',scanNum);
%    keyboard
    % This code is no longer necessary--should be removed once tested
    %if (length(stimfile) == 1) && isfield(stimfile{1},'myscreen') && isfield(stimfile{1}.myscreen,'traces')
    %  concatInfo.traces = stimfile{1}.myscreen.traces;
    %else
    %  disp(sprintf('Missing stimfile for scan %i',scanNum));
    %end
  % for subsequent scans, we are going to append
  else
    oldScanParams = viewGet(viewConcat,'scanParams',saveScanNum);
    oldScanParams.originalFileName{end+1} = filename;
    oldScanParams.originalGroupName{end+1} = baseGroupName;
    oldScanParams.totalJunkedFrames = totalJunkedFrames;
    viewConcat = saveTSeries(viewConcat,d.data,saveScanNum,oldScanParams,[],1);

    % This code is no longer necessary--should be removed once tested
    % now append traces (but only if we have a traces field
    % this way we only create traces if all of the stimfiles exist
    %if isfield(concatInfo,'traces')
    %  stimfile = viewGet(viewBase,'stimFile',scanNum);
    %if (length(stimfile) == 1) && isfield(stimfile{1},'myscreen') && isfield(stimfile{1}.myscreen,'traces')
    %	concatInfo = concatTraces(concatInfo,stimfile{1}.myscreen.traces);
    %  else
    %	concatInfo = rmfield(concatInfo,'traces')
    %	disp(sprintf('Missing stimfile for scan %i',scanNum));
    %  end
    %end
  end   
  % remember some info about where data comes from
  concatInfo.n = concatInfo.n+1;
  concatInfo.nifti(concatInfo.n) = hdr;

  % save the original path and filename
  concatInfo.filename{concatInfo.n} = filename;
  concatInfo.path{concatInfo.n} = path;

  % save which scan and volume every frame is from
  concatInfo.whichScan = [concatInfo.whichScan iscan*ones(1,d.dim(4))];
  concatInfo.whichVolume = [concatInfo.whichVolume 1:d.dim(4)];

  % save the junk frames
  concatInfo.junkFrames(concatInfo.n) = junkFrames;
  
  % keep the highpass filter used
  if isfield(d,'hipassfilter')
    concatInfo.hipassfilter{concatInfo.n} = d.hipassfilter;
  end
  
  % update the runtransition field, which keeps the beginning
  % and end volume of each run
  totalVolumes = length(concatInfo.whichVolume);
  concatInfo.runTransition(iscan,1) = totalVolumes-d.dim(4)+1;
  concatInfo.runTransition(iscan,2) = totalVolumes;

  % Update waitbar
  mrWaitBar(iscan/length(params.scanList),waitHandle);
end
mrCloseDlg(waitHandle);
toc;

% Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = concatTSeries(view,params);'];
tseriesdir = viewGet(viewConcat,'tseriesdir');
[pathstr,filename,ext,versn] = fileparts(fullfile(tseriesdir,tseriesFileName));
save(fullfile(pathstr,filename),'evalstr','params','concatInfo');

% Delete temporary viewBase and viewAverage
deleteView(viewBase);
deleteView(viewConcat);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat channels together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function concatInfo = concatTraces(concatInfo, newTraces)

% find last acquistition pulse in data set
lastacq = last(find(concatInfo.traces(1,:)));
% find spacing of acq pulses
acqspace = median(diff(find(concatInfo.traces(1,:))));
% find out where the new acq pulses start
firstacq = first(find(newTraces(1,:)));
% find loc of new acq pulses
newacqpos = firstacq:length(newTraces(1,:));
% find length of new acq pulses
newacqlen = length(newacqpos);
% figure out where to put new acq pulses
putnewacqpos = (lastacq+acqspace):(lastacq+acqspace+newacqlen-1);
% put the stim traces there too
concatInfo.traces(:,putnewacqpos) = newTraces(:,newacqpos);


function[d] = detrendTseries(d)
fprint('not functional yet')
 
