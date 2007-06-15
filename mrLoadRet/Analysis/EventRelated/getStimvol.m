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
%        
%             v.taskNum = 2;
%             v.phaseNum = 2;
%             v.segmentNum = 1;
%             v.varname = 'orientation';
%             d = getStimvol(d,'orientation');
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

% keep the varname that this was called with
d.varname = varname;

% depending on what kind of stimfile we have, we get the
% stimvolumes differently
for i = 1:length(d.stimfile)
  % get the stimvol for this particular stimfile
  switch d.stimfile{i}.filetype,
   case 'mgl',
    % if we have a stimtrace then get the variables from that
    if isfield(varname,'stimtrace') || isnumeric(varname)
      stimvol = getStimvolFromTraces(d.stimfile{i},varname);
    % otherwise get it using the varname
    else
      [stimvol d.stimNames] = getStimvolFromVarname(varname,d.stimfile{i}.myscreen,d.stimfile{i}.task);
    end
   case 'eventtimes',
      stimvol = getStimvolFromEventTimes(d.stimfile{i}.mylog, d.tr);
  end
  % set the stimnames if we don't have them already
  if ~isfield(d,'stimNames')
    for stimNameNum = 1:length(stimvol)
      d.stimNames{stimNameNum} = num2str(stimNameNum);
    end
  end
  % get how many junk frames we have
  if isfield(d,'junkFrames'),junkFrames = d.junkFrames;,else,junkFrames = zeros(length(d.stimfile),1);end
  % then shift all the stimvols by that many junkvols
  for nhdr = 1:length(stimvol)
    % subtract junkFrames
    stimvol{nhdr} = stimvol{nhdr}-junkFrames(i);
    % and get rid of anything less than 0
    stimvol{nhdr} = stimvol{nhdr}(stimvol{nhdr}>0);
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
if ieNotDefined('nhdr'),nhdr = max(stimraw);,end

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
stimvol = cell(1, nhdr);
for i = 1:nhdr
  stimvol{i} = round(stimfile.stimtimes_s{i}(:) / tr)';
end

