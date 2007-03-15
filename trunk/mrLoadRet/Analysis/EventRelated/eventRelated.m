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

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(params.scanNum)
  % make sure we are running on a set with a stimfile
  stimfile = viewGet(view,'stimfile',params.scanNum(scanNum));
  
  if isempty(stimfile)
    mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum(scanNum),params.groupName));
    return
  end

  % default is to not have an mgl variable name
  params.varname{params.scanNum(scanNum)} = '';

  % see if we have a stimfile from mgl, in which case we should
  % ask the user what the variable name is that they want ot use for the analysis
  if strfind(stimfile{1}.filetype,'mgl')
    [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
    % either ask the user for the single variable name, or if
    if (scanNum == 1) && (length(params.scanNum)>1)
      vals = mrDefaultParamsGUI({{'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)},...
			 {'sameForAll',1,'type=checkbox','Use the same variable name for all analyses'}});
      % user hit cancel
      if isempty(vals)
	return
      end
      % otherwise set the variable name
      params.varname{params.scanNum(scanNum)} = vals.varname;
      % and set it for all scans if called for
      if vals.sameForAll
	for i = 1:length(params.scanNum)
	  params.varname{params.scanNum(i)} = vals.varname;
	  scanNum = length(params.scanNum);
	end
      end
    else
      vals = mrDefaultParamsGUI({{'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)}});
      % user hit cancel
      if isempty(vals)
	return
      end
      % otherwise set the variable name
      params.varname{params.scanNum(scanNum)} = vals.varname;
      % if we have multiple scans, give the user the option of selecting
      % the same varaible for all scans
    end
  end
end
drawnow;

% create the parameters for the overlay
dateString = datestr(now);
r2.name = 'r2';
r2.function = 'eventRelated';
r2.groupName = params.groupName;
r2.reconcileFunction = 'mrDefaultParamsReconcile';
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
  numSlicesAtATime = floor(250000000/(8*numVolumes*prod(dims(1:2))));
  currentSlice = 1;
  ehdr = [];thisr2 = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % load the scan
    d = loadScan(view,scanNum,[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);;
    % keyboard
    % do any called for preprocessing
    % d =       (d,params.preprocess);
    % get the stim volumes, if there is a variable name used (for mgl)
    % pass that along as well
    d = getStimvol(d,params.varname{scanNum});
    % make a stimulation convolution matrix
    d = makescm(d,ceil(params.hdrlen/d.tr));
    % compute the estimated hemodynamic responses
    d = getr2(d);
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    thisr2 = cat(3,thisr2,d.r2);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = ehdr;
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
  erAnal.d{scanNum}.nhdr = d.nhdr;
  erAnal.d{scanNum}.hdrlen = d.hdrlen;
  erAnal.d{scanNum}.tr = d.tr;
  erAnal.d{scanNum}.stimvol = d.stimvol;
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;
end

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
function d = getStimvol(d,varname)

% check arguments
if ~any(nargin == [1 2])
  help getstimvol
  return
end

if ~exist('varname','var'),varname = '';,end

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
   case 'mgl',
    % if we have a variable name, then get the stimvols from that
    if ~isempty(varname)
      stimvol = getStimvolFromVarname(varname,d.stimfile{i}.myscreen,d.stimfile{i}.task);
    % otherwise get it the old style from the traces
    else
      stimvol = getStimvolFromTraces(d.stimfile{i});
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getStimvolFromVarname uses mgl/task encoded var names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvolOut = getStimvolFromVarname(varnameIn,myscreen,task)

% check arguments
if ~any(nargin == [3])
  help getStimvolFromVarname
  return
end

% first make varname into a cell array of cell arrays.
if isstr(varnameIn)
  varname{1}{1} = varnameIn;
else
  % see if this is a cell array of strings
  isArrayOfStr = 1;
  for i = 1:length(varnameIn)
    if ~isstr(varnameIn{i})
      isArrayOfStr = 0;
    end
  end
  %  now if it is an array of strings, then
  % we only have one condition
  if isArrayOfStr
    varname{1} = varnameIn;
  else
    for i = 1:length(varnameIn)
      if isstr(varnameIn{i})
	varname{i}{1} = varnameIn{i};
      else
	varname{i} = varnameIn{i};
      end
    end
  end
end
% get the task parameters
e = getTaskParameters(myscreen,task);
% now cycle over all array elements
for i = 1:length(varname)
  stimvol{i} = {};
  stimnames{i} = {};
  % cycle over all strings in varname
  for j = 1:length(varname{i})
    % get the value of the variable in question
    % on each trial
    [thisvarname varval enum] = getVarFromParameters(varname{i}{j},e);
    % see if it is a strict variable name
    if isempty(strfind(varname{i}{j},'='))
      % if it is then for each particular setting
      % of the variable, we make a stim type
      vartypes = unique(sort(varval));
      for k = 1:length(vartypes)
	stimvol{i}{end+1} = varval==vartypes(k);
	stimnames{i}{end+1} = sprintf('%s=%s',varname{i}{j},num2str(vartypes(k)));
      end
    end
  end
  % now we cycle through again, looking for any modifiers
  % that is, if we have something like varname=[1 2 3] then
  % it means that we only add the variable in, if it has
  % that var set appropriately. 
  for j = 1:length(varname{i})
    % get the value of the variable in question
    % on each trial
    [thisvarname varval enum] = getVarFromParameters(varname{i}{j},e);
    % see if it is a conditional variable, that is,
    % one that is like var=[1 2 3].
    if ~isempty(strfind(varname{i}{j},'='))
      [t,r] = strtok(varname{i}{j},'=');
      varcond = strtok(r,'=');
      % if it is then get the conditions
      varval = eval(sprintf('ismember(varval,%s)',varcond));
      % if we don't have any applied conditions applied, then
      % this is the condition
      if isempty(stimvol{i})
	stimvol{i}{1} = varval;
      else
	for k = 1:length(stimvol{i})
	  stimvol{i}{k} = stimvol{i}{k} & varval;
	  stimnames{i}{k} = sprintf('%s %s=%s',stimnames{i}{k},thisvarname,varcond);
	end
      end
    end
  end
end
% find stimvols and put into structure
stimvolOut = {};
for i = 1:length(stimvol)
  for j = 1:length(stimvol{i})
    % if we have a cell array of e that means there were multiple
    % tasks running, so we have to select the correct task
    if iscell(e)
      stimvolOut{end+1} = e{enum}.trialVolume(stimvol{i}{j});
    else
      stimvolOut{end+1} = e{enum}.trialVolume(stimvol{i}{j});
    end
  end
end

% check for non-unique conditions
if length(cell2mat(stimvolOut)) ~= length(unique(cell2mat(stimvolOut)))
  disp(sprintf('(getstimvol) Same trial in multiple conditions.'));
end

%%%%%%%%%%%%%%%%%%%%%%%%
% get a variable name from task parameters, could be
% in parameter or randVars
%%%%%%%%%%%%%%%%%%%%%%%%
function [varname varval enum] = getVarFromParameters(varname,e)

if ~iscell(e)
  olde = e;
  clear e;
  e{1} = olde;
end

varval = '';
enum = [];

% remove any equal sign
varname = strtok(varname,'=');
for i = 1:length(e)
  % see if it is a parameter
  if isfield(e{i}.parameter,varname)
    varval = e{i}.parameter.(varname);
    enum = i;
    % or a rand var
  elseif isfield(e{i},'randVars') && isfield(e{i}.randVars,varname)
    varval = e{i}.randVars.(varname);
    enum = i;
  end
end

if isempty(varval)
  disp(sprintf('(getstimvol) Could not find variable name %s',varname));
end 