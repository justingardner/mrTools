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

% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
%  params = eventRelatedGUI('groupName',viewGet(view,'groupName'));
  params.scanNum = selectScans(view);
  params = eventRelatedReconcileParams(viewGet(view,'groupName'),params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params = eventRelatedReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('eventRelated cancelled');
  return
end

% make sure we are running on a set with a stimfile
stimfile = viewGet(view,'stimfile',params.scanNum,params.groupNum);

if isempty(stimfile)
  mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum,params.groupName));
  return
end

% load the scan
d = loadScan(view,params.scanNum);
% keyboard
% do any called for preprocessing
% d =       (d,params.preprocess);
% get the stim volumes
d = getStimVol(d);
% make a stimulation convolution matrix
d = makescm(d,ceil(params.hdrlen/d.tr));
% compute the hemodynamic responses
d = getr2(d);

% create the r2 overlay
dateString = datestr(now);
r2.name = 'co';
r2.function = 'corAnal';
r2.groupName = params.groupName;
r2.reconcileFunction = 'eventRelatedReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
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
erAnal.reconcileFunction = 'eventRelatedReconcileParams';
erAnal.guiFunction = 'eventRelatedGUI';
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
function d = loadScan(view,scanNum)

% load parameters
d.tr = viewGet(view,'framePeriod',scanNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum);
d.dim = viewGet(view,'scandims',scanNum);
d.nFrames = viewGet(view,'nFrames',scanNum);
d.dim(4) = d.nFrames;

% Load data
mrDisp(sprintf('Loading scan %i from %s\n',scanNum,viewGet(view,'groupName')));
d.data = loadTSeries(view,scanNum,'all');
	
% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum);
d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);

% load dicom header
d.dicom = viewGet(view,'dicom',scanNum);

% load stimfile and set traces
% keyboard
stimfile = viewGet(view,'stimfile',scanNum);
if length(stimfile) == 1
  d.stimfile = stimfile{1};
  if ~isfield(d.stimfile, 'filetype'),
      d.stimfile.filetype = 'traces';
  end;
  switch d.stimfile.filetype,
      case 'traces',
        d.traces = stimfile{1}.traces;
        d.stimtrace = stimfile{1}.stimtrace;
        d.stimfile = rmfield(d.stimfile,'traces');
        % get acquisition times. 2 means volume acq (we don't know slice acq)
        d.acq = [0 2*(diff(d.traces(1,:))==1)];
        d.filetype = 'traces';
      case 'eventtimes',
        % stimtimes_s is the stimulus onset times in seconds wrt starting
        % the scan
        d.stimtimes_s = stimfile{1}.stimtimes_s;
        d.filetype = 'eventtimes';
      otherwise,
        mrMsgBox(sprintf('Unknown of invalid stimfile type %s',stimfile.filetype));
  end;
end
d.baseline = ones(d.dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getstimvol.m
%
%      usage: d = getStimVol(d)
%         by: justin gardner
%       date: 12/21/05
%    purpose: gets the stimulus vols from the traces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getStimVol(d,keepstim)

% check arguments
if ~any(nargin == [1 2])
  help getstimvol
  return
end
if exist('keepstim')~=1,keepstim = [];,end;
switch d.filetype,
    case 'traces',
        d = getStimVolFromTraces(d);
    case 'eventtimes',
        d = getStimVolFromEventTimes(d);
end;
% get the first volume
d = getfirstvol(d);

% check if we need to only keep some vols
if ~isempty(keepstim)
  if (min(keepstim) < 1) | (max(keepstim) > length(d.stimvol))
    mrDisp(sprintf('UHOH: Keepstim out of range. Ignoring'));
  else
    % go through and only keep the stim values asked for
    stimvol = [];
    for i = 1:length(keepstim)
      stimvol{i} = d.stimvol{keepstim(i)};
    end
    d.stimvol = stimvol;
  end
end


function d = getStimVolFromTraces(d)
% get the stim times
stimraw = d.traces(d.stimtrace,:);
stimraw(stimraw < 0) = 0;
stimtimes = find([0 (diff(stimraw)~=0)]);

% get the image number
acqnum = cumsum(d.acq>1);

% set the beginning acqnum to 1, so that
% any event that happens before the first
% acquistion pulse is assumed to happen
% during the first acquisition pulse.
acqnum(1:first(find(acqnum == 1))) = 1;

% sort into stimuli
nhdr = max(stimraw);
for i = 1:nhdr
  d.stimtimes{i} = stimtimes(stimraw(stimtimes) == i);
  d.pulselens(i) = i;
  d.stimvol{i} = acqnum(d.stimtimes{i});
end


function d = getStimVolFromEventTimes(d)
% sort into stimuli
nhdr = length(d.stimtimes_s);
for i = 1:nhdr
  d.stimtimes{i} = d.stimtimes_s{i};
  d.pulselens(i) = i;
  d.stimvol{i} = round(d.stimtimes_s{i} / d.tr);
end
