% concatTSeries.m
%
%      usage: concatTSeries()
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

% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
%  params = concatTSeriesGUI('groupName',viewGet(view,'groupName'));
  params = concatTSeriesReconcileParams(viewGet(view,'groupName'));
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params = concatTSeriesReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('concatTSeries cancelled');
  return
end

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
d.tr = viewGet(viewBase,'framePeriod',params.baseScan);
d.voxelSize = viewGet(viewBase,'scanvoxelsize',params.baseScan);
d.dim = viewGet(viewBase,'scandims',params.baseScan);
d.nFrames = viewGet(viewBase,'nFrames',params.baseScan);
d.dim(4) = d.nFrames;
for iscan = 1:length(params.scanList)
  if (viewGet(viewBase,'framePeriod',params.scanList(iscan)) ~= d.tr)
    mrWarnDlg('concatTSeries: These scans have different TR.');
  end
  if ~isequal(viewGet(viewBase,'scanvoxelsize',params.scanList(iscan)),d.voxelSize)
    mrErrorDlg('concatTSeries: Scans have different voxel sizes.');
  end
  if ~isequal(viewGet(viewBase,'scandims',params.scanList(iscan)),d.dim(1:3))
    mrErrorDlg('concatTSeries: Scans have different dimensions.');
  end
end

% initialize some things
concatInfo.n = 0;
concatInfo.whichScan = [];
concatInfo.whichVolume = [];

% Compute output volume
waitHandle = mrWaitBar(0,'Concatenating tSeries.  Please wait...');
for iscan = 1:length(params.scanList)
  scanNum = params.scanList(iscan);
  
  % Load it
  mydisp(sprintf('\nLoading scan %i from %s\n',scanNum,viewGet(viewBase,'groupName')));
  d.data = loadTSeries(viewBase,scanNum,'all');
	
  % Dump junk frames
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);
	
  % Compute transform
  if params.warp
    % *** Not fully tested yet ***
    baseXform = viewGet(view,'scanXform',baseScan,groupNum);
    scanXform = viewGet(view,'scanXform',scanNum,groupNum);
    % Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the
    % origin.
    shiftXform = eye(4);
    shiftXform(1:3,4) = - 1;
    M = inv(shiftXform) * inv(scanXform) * baseXform * shiftXform;
	
	
    % Warp the frames
    for frame = 1:d.nFrames
      d.data(:,:,:,frame) = warpAffine3(d.data(:,:,:,frame),M,NaN,0,params.interpMethod);
    end   
  end
  
  % do other  processing here
  if params.filterType == 1
    d = myhighpass(d,params.filterCutoff);
  end

  % convert to percent signal change
  if params.percentSignal
    d.mean = mean(d.data,4);
    disppercent(-inf, 'Converting to percent signal change');
    for i = 1:d.dim(4)
      d.data(:,:,:,i) = d.data(:,:,:,i)./d.mean;
      disppercent(i/d.dim(4));
    end
    disppercent(inf);
  end
  
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
    hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',params.baseScan));
    [viewConcat,tseriesFileName] = saveNewTSeries(viewConcat,d.data,scanParams,hdr);
    % get new scan number
    saveScanNum = viewGet(viewConcat,'nScans');
    
    % now load up channels
    stimfile = viewGet(viewBase,'stimFile',scanNum);
    if length(stimfile) == 1
      concatInfo.traces = stimfile{1}.traces;
    else
      disp(sprintf('Missing stimfile for scan %i',scanNum));
    end
  % for subsequent scans, we are going to append
  else
    oldScanParams = viewGet(viewConcat,'scanParams',saveScanNum);
    oldScanParams.originalFileName{end+1} = filename;
    oldScanParams.originalGroupName{end+1} = baseGroupName;
    viewConcat = saveTSeries(viewConcat,d.data,saveScanNum,oldScanParams,[],1);

    % now append traces (but only if we have a traces field
    % this way we only create traces if all of the stimfiles exist
    if isfield(concatInfo,'traces')
      stimfile = viewGet(viewBase,'stimFile',scanNum);
      if length(stimfile) == 1
	concatInfo = concatTraces(concatInfo,stimfile{1}.traces);
      else
	concatInfo = rmfield(concatInfo,'traces')
	disp(sprintf('Missing stimfile for scan %i',scanNum));
      end
    end
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


% Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = concatTSeries(view,params);'];
tseriesdir = viewGet(viewConcat,'tseriesdir');
[pathstr,filename,ext,versn] = fileparts(fullfile(tseriesdir,tseriesFileName));
save(fullfile(pathstr,filename),'evalstr','params','concatInfo');

% Delete temporary viewBase and viewAverage
deleteView(viewBase);
deleteView(viewConcat);

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
