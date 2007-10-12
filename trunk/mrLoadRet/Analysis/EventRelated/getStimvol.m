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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getStimvol(d,stimVariable,varargin)

% check arguments
if (nargin < 1)
  help getstimvol
  return
end

% evaluate other arguments
eval(evalargs(varargin));

if ~exist('stimVariable','var'),stimVariable = ''; end

% check to make sure we have a stimfile
if ~isfield(d,'stimfile')
  disp(sprintf('(getStimvol) stimFile is not loaded'));
  return
end

% convert simple string stimVariable into a structure
if isstr(stimVariable)
  v.varname = stimVariable;
else
  v = stimVariable;
end

% add on any of the other parameters that were passed in
if ~ieNotDefined('taskNum'),v.taskNum = taskNum;end
if ~ieNotDefined('phaseNum'),v.phaseNum = phaseNum;end
if ~ieNotDefined('segmentNum'),v.segmentNum = segmentNum;end

% keep the varname that this was called with
d.varname = v;

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
    if isfield(v,'stimtrace') || isnumeric(v)
      stimvol = getStimvolFromTraces(d.stimfile{i},v);
    else
      % otherwise get it using the name of the variable
      if exist('getStimvolFromVarname')~=2
	mrErrorDlg('(getStimvol) The function getStimvol is missing from your path. Make sure that mgl is in your path');
      end
      [stimvol d.stimNames] = getStimvolFromVarname(v,d.stimfile{i}.myscreen,d.stimfile{i}.task);
    end
   case 'eventtimes',
    stimvol = getStimvolFromEventTimes(d.stimfile{i}.mylog, d.tr/samplingf, impulse);
    if isfield(d.stimfile{i}, 'stimNames')
      d.stimNames = d.stimfile{i}.stimNames;
    end
   case 'afni',
    stimvol = getStimvolFromStimTS(d.stimfile{i});
    d.stimNames = d.stimfile{i}.stimNames;
  end
  % set the stimnames if we don't have them already
  if ~isfield(d,'stimNames')
    for stimNameNum = 1:length(stimvol)
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
    runlen = diff(d.concatInfo.runTransition(i,:))+1;
    if ~isempty(find(stimvol{nhdr}>runlen))
      disp(sprintf('(getStimvol) Removing %i event(s) from concatenated scan %i:%s since they happen after the last volume (%i) of the scan ',length(find(stimvol{nhdr}>runlen)),i,d.concatInfo.filename{i},runlen));
      stimvol{nhdr} = stimvol{nhdr}(stimvol{nhdr}<=runlen);
    end
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
      d.stimvol{nhdr} = [d.stimvol{nhdr} (stimvol{nhdr}+(d.concatInfo.runTransition(i,1)-1)*samplingf)];
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

stimTimes = find(stimfile.stimts == 1);
[nvols, nhdr] = size(stimfile.stimts);

stimvol = cell(1, nhdr);

for i = 1:nhdr
  stimvol{i}=find(stimfile.stimts(:,i))';
end
