%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getstimvol.m
%
%      usage: d = getStimvol(d,varname)
%         by: justin gardner
%       date: 12/21/05
%    purpose: gets the stimulus vols. varname
%             can be a variable name or can be 
%             a structure with fields: 
%             taskNum,phaseNum,segmentNum, and varname
%       e.g.: 
%             d = getStimvol(d,'orientation');
%             d = getStimvol(d,'orientation','taskNum=2');
%             d = getStimvol(d,'orintation','taskNum=2','phaseNum=2','segmentNum=1');
% 
%             d must have the field stimfile, tr, dim
%             optionally concatInfo, junkFrames, impulse, supersampling, 
% 
%             and optional arguments (impulse and supersampling for using eventtimes);
%
%             can also be used to get stimvol and stimNames directly w/out the
%             use of a d structure. Make sure to set the views curGroup and curScan
%             to the groupNum/scanNum that you want to get stimvols for.
%
%             v = newView;
%             v = viewSet(v,'curGroup',3);
%             v = viewSet(v,'curScan',1);
%             [stimvol stimNames var] = getStimvol(v,'orientation');
%            
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d stimNames var] = getStimvol(d,stimVariable,varargin)

% check arguments
if (nargin < 1)
  help getStimvol
  return
end

% evaluate other arguments
taskNum=[];phaseNum=[];segmentNum=[];stimfile=[];tr=[];nFrames=[];concatInfo=[];
junkFrames=[];impulse=[];supersampling=[];
getArgs(varargin,{'taskNum=[]','phaseNum=[]','segmentNum=[]','stimfile=[]','tr=[]','nFrames',[],'concatInfo',[],'junkFrames',[],'impulse',[],'supersampling',[],'verbose',true});

if ~exist('stimVariable','var'),stimVariable = ''; end

% if we are passed in a view instead
if isview(d)
  % change the names of the variable appropriately
  v = d;d = [];
  % set up d structure
  d.dim = viewGet(v,'scanDims');
  d.dim(end+1) = viewGet(v,'nFrames');
  d.stimfile = viewGet(v,'stimFile');
  d.tr = viewGet(v,'framePeriod');
  d.concatInfo = viewGet(v,'concatInfo');
  d.junkFrames = viewGet(v,'totalJunkedFrames');
  % and optional arguments (only necessary for eventtimes)
  if ~isempty(supersampling),d.supersampling = supersampling;end
  if ~isempty(impulse),d.impulse=impulse;end
  returnOnlyStimvol = 1;
else
  returnOnlyStimvol = 0;
end

% check to make sure we have a stimfile
if ~isfield(d,'stimfile')
  disp(sprintf('(getStimvol) stimFile is not loaded'));
  return
end

% check to make sure stimfile is not empty
if isempty(d.stimfile)
  disp(sprintf('(%s) No stimfiles found',mfilename));
  return
end
% convert simple string stimVariable into a structure
if ~isstruct(stimVariable)
  var.varname = stimVariable;
else
  var = stimVariable;
end

% add on any of the other parameters that were passed in
if ~ieNotDefined('taskNum'),var.taskNum = taskNum;end
if ~ieNotDefined('phaseNum'),var.phaseNum = phaseNum;end
if ~ieNotDefined('segmentNum'),var.segmentNum = segmentNum;end
if ~isfield(var,'verbose') var.verbose = verbose;end

% keep the varname that this was called with
d.varname = var;

% TR supersampling
if isfield(d,'supersampling')
  samplingf = d.supersampling;
  if ( samplingf~=floor(samplingf) ) || ( samplingf<1 )
    disp(sprintf('invalid TR supersampling factor (%g). using default (1.0)', samplingf));
    samplingf = 1;
  end
else
  samplingf = 1;
end

if isfield(d, 'impulse')
  impulse = d.impulse;
else
  impulse = 1;
end

%check if supersampling is supported for all runs
for i = 1:length(d.stimfile)
  if ~strcmp(d.stimfile{i}.filetype,'eventtimes') && ( samplingf~=1 )
    disp(sprintf('TR supersampling not supported for mgl or afni files. using default (1.0)'));
    samplingf = 1;
  end
end

% depending on what kind of stimfile we have, we get the
% stimvolumes differently
for i = 1:length(d.stimfile)
  % get the stimvol for this particular stimfile
  switch d.stimfile{i}.filetype,
   case 'mgl',
    % if we have a stimtrace then get the variables from that
    if isfield(var,'stimtrace') || isnumeric(var)
      stimvol = getStimvolFromTraces(d.stimfile{i},var);
    else
      % otherwise get it using the name of the variable
      if exist('getStimvolFromVarname')~=2
	mrErrorDlg('(getStimvol) The function getStimvol is missing from your path. Make sure that mgl is in your path');
      end
      [stimvol d.stimNames] = getStimvolFromVarname(var,d.stimfile{i}.myscreen,d.stimfile{i}.task);
    end
   case 'eventtimes',
    stimvol = getStimvolFromEventTimes(d.stimfile{i}.mylog, d.tr/samplingf, impulse);
    if isfield(d.stimfile{i}, 'stimNames')
      d.stimNames = d.stimfile{i}.stimNames;
    end
   case 'afni',
    stimvol = getStimvolFromStimTS(d.stimfile{i});
    if isfield(d.stimfile{i}, 'stimNames')
      d.stimNames = d.stimfile{i}.stimNames;
    end
   case 'stimvol',
    stimvol = d.stimfile{i}.stimvol;
    if isfield(d.stimfile{i}, 'stimNames')
      d.stimNames = d.stimfile{i}.stimNames;
    end
   otherwise
    disp(sprintf('(getStimvol) Unknown stimfile type: %s',d.stimfile{i}.filetype));
    d.stimvol = {};
    return
  end
  % set the stimnames if we don't have them already
  if ~isfield(d,'stimNames')
    d.stimNames = {};
  end
  % if we don't have a name for the stimulus, then give it a number
  for stimNameNum = 1:length(stimvol)
    if length(d.stimNames) < stimNameNum
      d.stimNames{stimNameNum} = num2str(stimNameNum);
    end
  end

  % get how many junk frames we have
  if isfield(d,'junkFrames'),junkFrames = d.junkFrames; else junkFrames = zeros(length(d.stimfile),1);end
  % then shift all the stimvols by that many junkvols
  for nhdr = 1:length(stimvol)
    % subtract junkFrames
    stimvol{nhdr} = stimvol{nhdr}-junkFrames(i)*samplingf;
    % and get rid of anything less than 0
    stimvol{nhdr} = stimvol{nhdr}(stimvol{nhdr}>0);
    % check for stimvol overrun
    if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
        runlen = d.dim(end)*samplingf;
    else
        runlen = diff(d.concatInfo.runTransition(i,:))*samplingf+1;
    end
    if ~isempty(find(stimvol{nhdr}>runlen,1))
      if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
         disp(sprintf('(getStimvol) Removing %i event(s) from scan since they happen after the last volume of the scan ',length(find(stimvol{nhdr}>runlen))));
      else
         disp(sprintf('(getStimvol) Removing %i event(s) from concatenated scan %i:%s since they happen after the last volume (%i) of the scan ',length(find(stimvol{nhdr}>runlen)),i,d.concatInfo.filename{i},runlen));
      end
      stimvol{nhdr} = stimvol{nhdr}(stimvol{nhdr}<=runlen);
    end
  end

  % if we have more than one stimfile, than we have to concatenate
  % together the stimvol. For this we are going to need to have
  % a concatInfo field
  if (i > 1)
    % check for valid concatInfo
    if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
      disp(sprintf('(getStimvol) No concatInfo found for multiple stimfiles'));
      return
    end
    % now add the number of volumes encountered from the previous runs
    % to the stimvol and concatenate on to general stimvol
    for nhdr = 1:length(stimvol)
      if length(d.stimvol) >= nhdr
	
	d.stimvol{nhdr} = [d.stimvol{nhdr} (stimvol{nhdr}+(d.concatInfo.runTransition(i,1)-1)*samplingf)];
      else
	d.stimvol{nhdr} = (stimvol{nhdr}+(d.concatInfo.runTransition(i,1)-1)*samplingf);
      end
    end
    % on first file, we just set stimvol in the d field
  else
    d.stimvol = stimvol;
  end
end

% return the actual TR supersampling factor
d.supersampling = samplingf;

% update the eventRelatedVarname
if isfield(d,'eventRelatedVarname')
  if isfield(d.varname,'varname')
    d.eventRelatedVarname = d.varname.varname;
  else
    d.eventRelatedVarname = d.varname;
  end
end

% change into return names
stimNames = d.stimNames;
if returnOnlyStimvol
  d = d.stimvol;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old way of getting stim vols
% from traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvol = getStimvolFromTraces(stimfile,stimtrace)

% if passed in a structure then get the stimtrace field
if isstruct(stimtrace)
  varname = stimtrace;
  stimtrace = varname.stimtrace;
  % if we are passed in how many response types there are, get it
  if isfield(varname,'nhdr')
    nhdr = varname.nhdr;
  end
end

if exist('stimtrace','var')
  stimfile.myscreen.stimtrace = stimtrace;
end

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
if ieNotDefined('nhdr'),nhdr = max(stimraw); end

stimvol = cell(1, nhdr);
for i = 1:nhdr
  thisStimtimes = stimtimes(stimraw(stimtimes) == i);
  stimvol{i} = acqnum(thisStimtimes);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimvol from Farshad's file type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvol = getStimvolFromEventTimes(stimfile,tr,impulse)

% sort into stimuli
nhdr = length(stimfile.stimtimes_s);
stimvol = cell(1, nhdr);

if ~impulse && isfield(stimfile, 'stimdurations_s')
  for i = 1:nhdr
    z = zeros(1, 1+ceil(max(stimfile.stimtimes_s{i})/tr));
    for event=1:length(stimfile.stimtimes_s{i})
      onset = max(stimfile.stimtimes_s{i}(event), 0)/tr;
      offset = max(onset, onset+stimfile.stimdurations_s{i}(event)/tr-0.49);
      z(1+[floor(onset):floor(offset)])=1;
    end
    stimvol{i} = unique(find(z));
  end
else
  for i = 1:nhdr
    stimvol{i} = 1+floor(stimfile.stimtimes_s{i}(:) / tr)';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimvol from AFNI file type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvol = getStimvolFromStimTS(stimfile)

% stimTimes = find(stimfile.stimts == 1);
[nvols, nhdr] = size(stimfile.stimts);

stimvol = cell(1, nhdr);

for i = 1:nhdr
  stimvol{i}=find(stimfile.stimts(:,i))';
end
