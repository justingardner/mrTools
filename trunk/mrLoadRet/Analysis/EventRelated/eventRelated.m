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

paramsInfo = {...
    {'groupName',viewGet(view,'groupNames'),'Name of group from which to do eventRelated analysis'},...
    {'description','Event related analysis of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'hdrlen',25,'Length of response in seconds to calculate'}...
    {'preprocess','','String of extra commands for preprocessing (see wiki for details)'}...
};

mrGlobals;

% First get parameters
if ieNotDefined('params')
  % Get parameter values
  params = mrParamsDialog(paramsInfo);
  % if empty user hit cancel
  if isempty(params),return,end
  % get scans
  view = viewSet(view,'groupName',params.groupName);
  params.scanNum = selectScans(view);
  if isempty(params.scanNum),return,end
  % check parameters
  params = mrParamsReconcile(viewGet(view,'groupName'),params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;  
  params = mrParamsReconcile(viewGet(view,'groupName'),params);
end

% Abort if params empty
if ieNotDefined('params'),return,end

% make sure the group is set properly
view = viewSet(view,'groupName',params.groupName);

%g et the name of the variable to do the evented related analysis on
params = getEventRelatedVarname(view,params);

% create the parameters for the overlay
dateString = datestr(now);
r2.name = 'r2';
r2.function = 'eventRelated';
r2.groupName = params.groupName;
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
    d = getStimvol(d,params.eventRelatedVarname{scanNum});
    if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    % do any called for preprocessing
    d = eventRelatedPreProcess(d,params.preprocess);
    % make a stimulation convolution matrix
    d = makescm(d,ceil(params.hdrlen/d.tr));
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
erAnal.name = 'erAnal';  % This can be reset by editAnalysisGUI
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'mrParamsReconcile';
erAnal.guiFunction = 'eventRelatePlot';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the variable name that the user wants
% to do the event related analysis on, puts up a gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = getEventRelatedVarname(view,params);

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(params.scanNum)
  % make sure we are running on a set with a stimfile
  stimfile = viewGet(view,'stimfile',params.scanNum(scanNum));
  
  if isempty(stimfile)
    mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum(scanNum),params.groupName));
    return
  end

  % see if we have a stimfile from mgl, in which case we should
  % ask the user what the variable name is that they want ot use for the analysis
  if strfind(stimfile{1}.filetype,'mgl')

    % check to see what style this is, if the task variable does
    % not have a segmentTrace then it mus be an old style, in which
    % we used channels
    if isfield(stimfile{1}.myscreen,'traces') && ...
	  ((iscell(stimfile{1}.task{1}) && ~isfield(stimfile{1}.task{1}{1},'segmentTrace')) ||...
	   (~iscell(stimfile{1}.task{1}) && ~isfield(stimfile{1}.task{1},'segmentTrace')))
      % this is the old style, get the stimtrace number
      taskVarParams{1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
    else
      % this is the new tyle, ask for a variable name
      [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
      taskVarParams = {};
      % if there is more than one task, then ask the user for that
      if length(stimfile{1}.task)>1
	taskVarParams{end+1} = {'taskNum',num2cell(1:length(stimfile{1}.task)),'The task you want to use'};
      end
      % if there are multiple phases, then ask for that
      maxPhaseNum = 0;
      for tnum = 1:length(stimfile{1}.task)
	phaseNum{tnum} = num2cell(1:length(stimfile{1}.task{tnum}));
	maxPhaseNum = max(maxPhaseNum,length(stimfile{1}.task{tnum}));
      end
      if maxPhaseNum > 1
	if length(stimfile{1}.task) == 1
	  taskVarParams{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
	else
	  taskVarParams{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
	end
      end
      % set up to get the variable name from the user
      taskVarParams{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)};
    end
    
    % give the option to use the same variable for all
    if (scanNum == 1) && (length(params.scanNum)>1)
      taskVarParams{end+1} = {'sameForAll',1,'type=checkbox','Use the same variable name for all analyses'};
    end
    % either ask the user for the single variable name, or if
    taskVarParams = mrParamsDialog(taskVarParams);
    % user hit cancel
    if isempty(taskVarParams),return, end
    % check if we were passed a cell array, if we were
    % then convert it
    if (isfield(taskVarParams,'varname') && isstr(taskVarParams.varname) && ...
	(length(taskVarParams.varname) > 1) && (taskVarParams.varname(1) == '{'))
      taskVarParams.varname = eval(taskVarParams.varname);
    end
    % check the variables--deal with sameForAll
    taskVarParamsFieldnames = fieldnames(taskVarParams);
    for tnum = 1:length(taskVarParamsFieldnames)
      if ~strcmp(taskVarParamsFieldnames{tnum},'sameForAll') & ~strcmp(taskVarParamsFieldnames{tnum},'paramInfo')
	% set it in this params
	params.eventRelatedVarname{params.scanNum(scanNum)} = taskVarParams.(taskVarParamsFieldnames{tnum});
	% and set it for all scans if called for
	if isfield(taskVarParams,'sameForAll') && taskVarParams.sameForAll
	  for i = 1:length(params.scanNum)
	    params.eventRelatedVarname{i} = taskVarParams.(taskVarParamsFieldnames{tnum});
	  end
	end
      end
    end
  end
  % break out of for look if we have sameForAll set
  if isfield(taskVarParams,'sameForAll') && taskVarParams.sameForAll
    break
  end
  taskVarParams = {};
end
