% eventRelated.m
%
%      usage: view = eventRelated(view,params)
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%
function view = eventRelated(view,params)

% check arguments
if ~any(nargin == [1 2])
  help eventRelated
  return
end

paramsInfo = {...
    {'groupName',viewGet(view,'groupNames'),'Name of group from which to do eventRelated analysis'},...
    {'description','Event related analysis of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'hdrlen',25,'Length of response in seconds to calculate'}...
    {'preprocess','','String of extra commands for preprocessing (see wiki for details)'}...
};

% First get parameters
if ieNotDefined('params')
  % Get parameter values
  params = mrDefaultParamsGUI(paramsInfo);
  % if empty user hit cancel
  if isempty(params),return,end
  % get scans
  view = viewSet(view,'groupName',params.groupName);
  params.scanNum = selectScans(view);
  if isempty(params.scanNum),return,end
  % check parameters
  params = mrDefaultParamsReconcile(viewGet(view,'groupName'),params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;  
  params = mrDefaultParamsReconcile(viewGet(view,'groupName'),params);
end

% Abort if params empty
if ieNotDefined('params'),return,end

% make sure the group is set properly
view = viewSet(view,'groupName',params.groupName);

% make sure we are running on a set with a stimfile
stimfile = viewGet(view,'stimfile',params.scanNum);

if isempty(stimfile)
  mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum,params.groupName));
  return
end

% decide how many slices to do at a time, this is done
% simply to save memory -- currently our system is limited
% to 2G of memory and for large concatenations, you need
% to break up the analysis into smaller portions of the data
numSlices = viewGet(view,'nSlices',params.scanNum);
numVolumes = viewGet(view,'nFrames',params.scanNum);
dims = viewGet(view,'dims',params.scanNum);
% choose how many slices based on trying to keep a certain
% amount of data in the memory
numSlicesAtATime = floor(250000000/(8*numVolumes*prod(dims(1:2))));
currentSlice = 1;
ehdr = [];r2 = [];

for i = 1:ceil(numSlices/numSlicesAtATime)
  % load the scan
  d = loadScan(view,params.scanNum,[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);;
  % keyboard
  % do any called for preprocessing
  % d =       (d,params.preprocess);
  % get the stim volumes
  d = getStimvol(d);
  % make a stimulation convolution matrix
  d = makescm(d,ceil(params.hdrlen/d.tr));
  % compute the hemodynamic responses
  d = getr2(d);
  % update the curernt slice we are working on
  currentSlice = currentSlice+numSlicesAtATime;
  % cat with what has already been computed for other slices
  ehdr = cat(3,ehdr,d.ehdr);
  r2 = cat(3,r2,d.r2);
end

% now put all the data from all the slices into the structure
d.ehdr = ehdr;
d.r2 = r2;
d.dim(3) = size(d.r2,3);
clear r2;

% create the r2 overlay
dateString = datestr(now);
r2.name = 'co';
r2.function = 'corAnal';
r2.groupName = params.groupName;
r2.reconcileFunction = 'mrDefaultParamsReconcile';
r2.data = cell(1,viewGet(view,'nScans'));
r2.data{params.scanNum} = d.r2;
r2.date = dateString;
r2.params = params;
r2.range = [0 1];
r2.clip = [0 1];
r2.colormap = hot(256);
r2.alpha = 1;
r2.interrogator = 'eventRelatedPlot';

% install analysis
erAnal.name = 'erAnal';  % This can be reset by editAnalysisGUI
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'mrDefaultParamsReconcile';
erAnal.guiFunction = 'eventRelatePlot';
erAnal.params = params;
erAnal.overlays = r2;
erAnal.curOverlay = 1;
erAnal.date = dateString;
erAnal.ehdr = d.ehdr;
erAnal.nhdr = d.nhdr
d = rmfield(d,'data');
erAnal.d = d;
view = viewSet(view,'newAnalysis',erAnal);

% Save it
saveAnalysis(view,erAnal.name);

%%%%%%%%%%%%%%%%%%%
% loadScan
%%%%%%%%%%%%%%%%%%%
function d = loadScan(view,scanNum,sliceNum)

% default to loading all slices
if ~exist('sliceNum','var'),sliceNum = [];end
  
% load parameters
d.ver = 4.5;
d.tr = viewGet(view,'framePeriod',scanNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum);
d.nFrames = viewGet(view,'nFrames',scanNum);
d.dim(4) = d.nFrames;
d.filename = viewGet(view,'tseriesfile',scanNum);
d.filepath = viewGet(view,'tseriespathstr',scanNum);
[d.fullpath d.expname] = fileparts(pwd);

% Load data
mrDisp(sprintf('Loading scan %i from %s slices=[%s]\n',scanNum,viewGet(view,'groupName'),num2str(sliceNum)));
d.data = loadTSeries(view,scanNum,sliceNum);
	
% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum);
d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);

% load dicom header
d.dicom = viewGet(view,'dicom',scanNum);

% load stimfile and set traces
d.stimfile = viewGet(view,'stimfile',scanNum);

% get any concat info
d.concatInfo = viewGet(view,'concatInfo',scanNum);

% get dimensions
d.dim = size(d.data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getstimvol.m
%
%      usage: d = getStimvol(d)
%         by: justin gardner
%       date: 12/21/05
%    purpose: gets the stimulus vols from the traces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getStimvol(d)

% check arguments
if ~any(nargin == [1])
  help getstimvol
  return
end

% check to make sure we have a stimfile
if ~isfield(d,'stimfile')
  disp(sprintf('(getStimvol) stimFile is not loaded'));
  return
end

% depending on what kind of stimfile we have, we get the
% stimvolumes differently
for i = 1:length(d.stimfile)
  % get the stimvol for this particular stimfile
  switch d.stimfile{i}.filetype,
   case 'mgltraces',
     stimvol = getStimvolFromTraces(d.stimfile{i});
   case 'eventtimes',
     stimvol = getStimvolFromEventTimes(d.stimfile{i},d.tr);
  end
  % if we have more than one stimfile, than we have to concatenate
  % together the stimvol. For this we are going to need to have
  % a concatInfo field
  if (i > 1)
    % check for valid concatInfo
    if ~isfield(d,'concatInfo')
      disp(sprintf('(getStimvol) No concatInfo found for multiple stimfiles'));
      return
    end
    % now add the number of volumes encountered from the previous runs
    % to the stimvol and concatenate on to general stimvol
    for nhdr = 1:length(stimvol)
      d.stimvol{nhdr} = [d.stimvol{nhdr} (stimvol{nhdr}+d.concatInfo.runTransition(i,1)-1)];
    end
  % on first file, we just set stimvol in the d field
  else
    d.stimvol = stimvol;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old way of getting stim vols
% from traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvol = getStimvolFromTraces(stimfile)

% get acquisition pulses
acq = [0 2*(diff(stimfile.myscreen.traces(1,:))==1)];

% get the stim times
stimraw = stimfile.myscreen.traces(stimfile.myscreen.stimtrace,:);
stimraw(stimraw < 0) = 0;
stimtimes = find([0 (diff(stimraw)~=0)]);

% get the image number
acqnum = cumsum(acq>1);

% set the beginning acqnum to 1, so that
% any event that happens before the first
% acquistion pulse is assumed to happen
% during the first acquisition pulse.
acqnum(1:first(find(acqnum == 1))) = 1;

% sort into stimuli
nhdr = max(stimraw);
for i = 1:nhdr
  thisStimtimes = stimtimes(stimraw(stimtimes) == i);
  stimvol{i} = acqnum(thisStimtimes);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimvol from Farshad's file type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvol = getStimvolFromEventTimes(stimfile,tr)

% sort into stimuli
nhdr = length(stimfile.stimtimes_s);
for i = 1:nhdr
  stimtimes{i} = stimfile.stimtimes_s{i};
  stimvol{i} = round(stimtimes_s{i} / tr);
end
