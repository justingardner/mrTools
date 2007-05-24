% concatTSeries.m
%
%        $Id$
%      usage: concatTSeries(view, params)
%         by: justin gardner
%       date: 10/12/06
%    purpose: 
%
function view = concatTSeries(view,params)

% check arguments
if ~any(nargin == [1 2])
  help concatTSeries
  return
end

% description of paramaters (used by mrParamsDialog functions)
paramsInfo = {...
    {'groupName',viewGet(view,'groupNames'),'Name of group from which to make concatenation'},...
    {'newGroupName','Concatenation','Name of group that will be created'},...
    {'description','Concatenation of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'filterType',1,'minmax=[0 1]','incdec=[-1 1]','Which filter to use, for now you can only turn off highpass filtering'},...
    {'filterCutoff',0.01,'minmax=[0 inf]','Highpass filter cutoff in Hz'},...
    {'percentSignal',1,'type=checkbox','Convert to percent signal change'},...
    {'warp',0,'type=checkbox','Warp images based on alignment (not implemented yet)'},...
    {'warpInterpMethod',{'nearest','bilinear'},'Interpolation method for warp (not implemented yet)','contingent=warp'}
	     };
% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = mrParamsDialog(paramsInfo);
  % no params means user hit cancel
  if isempty(params),return,end
  % select scans
  view = viewSet(view, 'groupName', params.groupName);
  params.scanList = selectScans(view);
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
for iscan = 1:length(params.scanList)
  if (viewGet(viewBase,'framePeriod',params.scanList(iscan)) ~= d.tr)
    mrWarnDlg(sprintf('concatTSeries: These scans have different TR. (%0.4f vs %0.4f)',viewGet(viewBase,'framePeriod',params.scanList(iscan)),d.tr));
  end
  if ~isequal(viewGet(viewBase,'scanvoxelsize',params.scanList(iscan)),d.voxelSize)
    disp('(concatTSeries) Scans have different voxel sizes.');
  end
  if ~isequal(viewGet(viewBase,'scandims',params.scanList(iscan)),d.dim(1:3))
    disp('(concatTSeries) Scans have different dimensions.');
  end
end

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
	
  % Compute transform
  if params.warp
    % get base and scan xforms
    baseXform = viewGet(view,'scanXform',params.warpBaseScan,groupNum);
    scanXform = viewGet(view,'scanXform',scanNum,groupNum);
    % Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the
    % origin.
    shiftXform = shiftOriginXform;
    swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
    M = swapXY * inv(shiftXform) * inv(scanXform) * baseXform * shiftXform * swapXY;

    % display transformation
    disp(sprintf('Transforming time series with transformation: '));
    for rownum = 1:4
      disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
    end
    % Warp the frames
    for frame = 1:d.nFrames
      d.data(:,:,:,frame) = warpAffine3(d.data(:,:,:,frame),M,NaN,0,params.warpInterpMethod);
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
  baseGroupName = viewGet(view,'groupName');

  % Save TSeries (using header of 1st scan on scanList as the template for
  % the nifti header), and add it as a new scan.
  if iscan == 1
    scanParams.fileName = [];
    scanParams.junkFrames = 0;
    scanParams.nFrames = d.nFrames;
    scanParams.description = params.description;
    scanParams.originalFileName{1} = filename;
    scanParams.originalGroupName{1} = baseGroupName;
    hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',params.scanList(1)));
    % data *MUST* be written out as float32 b/c of the small values-epm
    hdr.datatype = 16;
    % if we are warping, then we need to change the sform to the
    % sform of the scan we warped to
    if params.warp
      hdr.sform44 = viewGet(view,'scanXform',params.warpBaseScan,groupNum);
    end
    [viewConcat,tseriesFileName] = saveNewTSeries(viewConcat,d.data,scanParams,hdr);
    % get new scan number
    saveScanNum = viewGet(viewConcat,'nScans');
    
    % now load up channels
    stimfile = viewGet(viewBase,'stimFile',scanNum);
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
 
