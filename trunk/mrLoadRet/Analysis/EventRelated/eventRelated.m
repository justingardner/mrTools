% eventRelated.m
%
%      usage: view = eventRelated(view,params)
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%
function [view d] = eventRelated(view,params)

d = [];

% check arguments
if ~any(nargin == [1 2])
  help eventRelated
  return
end

mrGlobals;

% First get parameters
if ieNotDefined('params')
  % put up the gui
  params = eventRelatedGUI;
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params = mrParamsReconcile([],params);
end

% Abort if params empty
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

% create the parameters for the overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'eventRelated';
r2.reconcileFunction = 'mrParamsReconcile';
r2.data = cell(1,viewGet(view,'nScans'));
r2.date = dateString;
r2.params = params;
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';

tic
for scanNum = params.scanNum
  % decide how many slices to do at a time, this is done
  % simply to save memory -- currently our system is limited
  % to 2G of memory and for large concatenations, you need
  % to break up the analysis into smaller portions of the data
  numSlices = viewGet(view,'nSlices',scanNum);
  numVolumes = viewGet(view,'nFrames',scanNum);
  dims = viewGet(view,'dims',scanNum);
  % choose how many slices based on trying to keep a certain
  % amount of data in the memory
  if isfield(MLR.prefs,'maxBlocksize')
    maxBlocksize = MLR.prefs.maxBlocksize;
  else
    maxBlocksize = 250000000;
  end
  numSlicesAtATime = floor(maxBlocksize/(8*numVolumes*prod(dims(1:2))));
  currentSlice = 1;
  ehdr = [];ehdrste = [];thisr2 = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % load the scan
    d = loadScan(view,scanNum,[],[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);;
    % get the stim volumes, if empty then abort
    d = getStimvol(d,params.scanParams{scanNum}.varname);
    if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    % do any called for preprocessing
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
    % make a stimulation convolution matrix
    d = makescm(d,ceil(params.scanParams{scanNum}.hdrlen/d.tr));
    % compute the estimated hemodynamic responses
    d = getr2(d);
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    ehdrste = cat(3,ehdrste,d.ehdrste);
    thisr2 = cat(3,thisr2,d.r2);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = ehdr;
  d.ehdrste = ehdrste;
  d.r2 = thisr2;
  d.dim(3) = size(d.r2,3);

  % save the r2 overlay
  r2.data{scanNum} = d.r2;

  % save other eventRelated parameters
  erAnal.d{scanNum}.ver = d.ver;
  erAnal.d{scanNum}.filename = d.filename;
  erAnal.d{scanNum}.filepath = d.filepath;
  erAnal.d{scanNum}.dim = d.dim;
  erAnal.d{scanNum}.ehdr = d.ehdr;
  erAnal.d{scanNum}.ehdrste = d.ehdrste;
  erAnal.d{scanNum}.nhdr = d.nhdr;
  erAnal.d{scanNum}.hdrlen = d.hdrlen;
  erAnal.d{scanNum}.tr = d.tr;
  erAnal.d{scanNum}.stimvol = d.stimvol;
  erAnal.d{scanNum}.stimNames = d.stimNames;
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;
end
toc

% install analysis
erAnal.name = params.saveName;
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'mrParamsReconcile';
erAnal.guiFunction = 'eventRelateGUI';
erAnal.params = params;
erAnal.overlays = r2;
erAnal.curOverlay = 1;
erAnal.date = dateString;
view = viewSet(view,'newAnalysis',erAnal);

% Save it
saveAnalysis(view,erAnal.name);

% for output
if nargout > 1
  for i = 1:length(d)
    erAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(erAnal.d) == 1
    d = erAnal.d{1}
  else
    d = erAnal.d;
  end
end

