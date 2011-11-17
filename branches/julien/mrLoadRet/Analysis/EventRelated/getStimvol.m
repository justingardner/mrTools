%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getstimvol.m
%
%      usage: d = getStimvol(d or view,varname)
%         by: justin gardner, modified by julien besle (22/03/2010)
%       date: 12/21/05
%    purpose: gets the stimulus vols. varname
%             can be a variable name or can be 
%             a structure with fields: 
%             taskNum,phaseNum,segmentNum, and varname
%              $Id$
%       e.g.: 
%             d = getStimvol(d,'orientation');
%             d = getStimvol(d,'orientation','taskNum=2');
%             d = getStimvol(d,'orientation','taskNum=2','phaseNum=2','segmentNum=1');
% 
%             d must have the field stimfile, tr, dim
%             optionally concatInfo, junkFrames 
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
junkFrames=[];returnOnlyStimvol=[];
getArgs(varargin,{'taskNum=[]','phaseNum=[]','segmentNum=[]','stimfile=[]','tr=[]','nFrames',[],'concatInfo',[],'junkFrames',[]});

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
  if isempty(returnOnlyStimvol),returnOnlyStimvol = 1;end
else
  if isempty(returnOnlyStimvol),returnOnlyStimvol = 0;end
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

% keep the varname that this was called with
d.varname = var;

% Design super-sampling for cases where the duration/time of the stimulation is not a multiple of TR (only supported for eventtimes file types)
time_remainders = [];
if isfield(var,'forceStimOnSampleOnset') && ~var.forceStimOnSampleOnset
   %compute the remainders of all stim times divided by d.tr
   for i_file = 1:length(d.stimfile)
      if strfind(d.stimfile{i_file}.filetype,'eventtimes') 
         for i_type = 1:length(d.stimfile{i_file}.mylog.stimtimes_s)
            time_remainders = [time_remainders rem(d.stimfile{i_file}.mylog.stimtimes_s{i_type},d.tr)];
         end
      else
         disp(sprintf('(getStimvol) off TR stim onset is not supported for mgl or afni files (%s). Forcing stim on (sub)sample onset',d.stimfile{i_file}.filename));
         time_remainders = [];
         break;
      end
   end
end
duration_remainders = [];
if isfield(var,'stimDuration')
   if ischar(var.stimDuration) && strcmp(var.stimDuration,'fromFile')
      %compute the remainders of all stim durations divided by tr
      for i_file = 1:length(d.stimfile)
         if strfind(d.stimfile{i_file}.filetype,'eventtimes') 
            if isfield(d.stimfile{i_file}.mylog,'stimdurations_s') 
               for i_type = 1:length(d.stimfile{i_file}.mylog.stimtimes_s)
                  duration_remainders = [duration_remainders rem(d.stimfile{i_file}.mylog.stimdurations_s{i_type},d.tr)];
               end
            else
               disp(sprintf('(getStimvol) Field stimdurations_s missing in file %s. Using stimDuration = 1 (sub)frame',d.stimfile{i_file}.filename));
               duration_remainders = [];
               var.stimDuration = d.tr;
               break;
            end
         else
            disp(sprintf('(getStimvol) Stim duration from files not supported for mgl or afni files (%s). Using stimDuration = 1 (sub)frame',d.stimfile{i_file}.filename));
            duration_remainders = [];
            var.stimDuration = d.tr;
            break;
         end
      end
   else
      if ischar(var.stimDuration) 
         var.stimDuration = eval(var.stimDuration);
      end
      duration_remainders = [duration_remainders rem(var.stimDuration,d.tr)];
   end
else
   var.stimDuration = d.tr;
end

%find the least supersampling factor (= greatest common factor of all possible durations/times of the event)
if isfield(var,'estimationSupersampling') %if sub-sample estimation is required
  estimationSupersampling = var.estimationSupersampling;
else
  estimationSupersampling=1;
end

resolutionLimit = .05; %temporal resolution limit in s
remainders = [([time_remainders duration_remainders])*estimationSupersampling d.tr];
remainders = round(remainders/resolutionLimit); %temporal resolution limited to .01 s
remainders(remainders==0) = d.tr/resolutionLimit;    %replace zeros by TRs
remainders = unique(remainders);
greatest_common_factor = 0;
for i = 1:length(remainders)
   greatest_common_factor = gcd(greatest_common_factor,remainders(i));
   if greatest_common_factor ==1
      break
   end
end
designSupersampling = round(d.tr*estimationSupersampling/(greatest_common_factor*resolutionLimit));

%check if designSupersampling is supported for all runs
if designSupersampling~=1 || ischar(var.stimDuration) || var.stimDuration ~= d.tr
   for i = 1:length(d.stimfile)
      if ~strcmp(d.stimfile{i}.filetype,'eventtimes')
         disp(sprintf('(getStimvol) TR supersampling not supported for mgl or afni files. using default (1.0)'));
         designSupersampling = 1;
      end
   end
end


% depending on what kind of stimfile we have, we get the
% stimvolumes differently
for i = 1:length(d.stimfile)
  % get the stimvol for this particular stimfile
  stimNames = {};
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
      [stimvol stimNames d.trialNum{i}] = getStimvolFromVarname(var,d.stimfile{i}.myscreen,d.stimfile{i}.task);
    end
   case 'eventtimes',
    [stimvol,stimDurations] = getStimvolFromEventTimes(d.stimfile{i}.mylog, d.tr/designSupersampling, var.stimDuration);
    if isfield(d.stimfile{i}, 'stimNames')
      stimNames = d.stimfile{i}.stimNames;
    end
   case 'afni',
    stimvol = getStimvolFromStimTS(d.stimfile{i});
    if isfield(d.stimfile{i}, 'stimNames')
      stimNames = d.stimfile{i}.stimNames;
    end
   case 'stimvol',
    stimvol = d.stimfile{i}.stimvol;
    % validate the stimvols (i.e. they need to be 1xn)
    for j = 1:length(stimvol)
      arraysize = size(stimvol{j});
      if arraysize(1)>1
	if arraysize(2)==1
	  disp(sprintf('(getStimvol) Stimvol{%i} converted from %ix1 to 1x%i array',j,arraysize(1),arraysize(1)));
	  stimvol{j} = stimvol{j}';
	else
	  mrErrorDlg(sprintf('(getStimvol) Stimvol{%i} is %ix%i must be 1xn',j,arraysize(1),arraysize(2)));
	end
      end
    end
    % check for stimNames
    if isfield(d.stimfile{i}, 'stimNames')
      stimNames = d.stimfile{i}.stimNames;
    end
   otherwise
    disp(sprintf('(getStimvol) Unknown stimfile type: %s',d.stimfile{i}.filetype));
    d.stimvol = {};
    return
  end
  % if we don't have a name for the stimulus, then give it a number
  for stimNameNum = 1:length(stimvol)
    if length(stimNames) < stimNameNum
      stimNames{stimNameNum} = num2str(stimNameNum);
    end
  end
  
  %if stimDurations does not exist, set all durations to value user
  %entered in the GUI in TRs
  if ieNotDefined('stimDurations')
    for iStim = 1:length(stimvol)
      stimDurations{iStim} = round(stimVariable.stimDuration/d.tr)*ones(size(stimvol{iStim}));;
    end
  end

 
%   % compare to make sure we have the same stimulus names as last time
%   if ~isempty(lastStimNames)
%     if ~isequal(d.stimNames,lastStimNames)
%       % display what doesn't match
%       missingLastTime = setdiff(d.stimNames,lastStimNames);
%       for iMissing = 1:length(missingLastTime)
% 	disp(sprintf('(getStimvol) **** Scan %i is missing stimulus type: %s ****',i-1,missingLastTime{iMissing}));
%       end
%       missingThisTime = setdiff(lastStimNames,d.stimNames);
%       for iMissing = 1:length(missingThisTime)
% 	disp(sprintf('(getStimvol) **** Scan %i is missing stimulus type: %s ****',i,missingThisTime{iMissing}));
%       end
%       disp(sprintf('(getStimvol) This likely happened because you have not specified all the possible values a variable can take. For a randVar, you should add a field with the same name as the variable except with an _ after it (i.e. for randVar.calculated.myVar add randVar.calculated.myVar_ ) and set the variable to all the possible values that the variable can be set to. See the function addCalculatedVar for more info.'));
%       keyboard
%     end
%   end
%   lastStimNames = d.stimNames;

  % get how many junk frames we have
  if isfield(d,'junkFrames'),junkFrames = d.junkFrames; else junkFrames = zeros(length(d.stimfile),1);end
  % then shift all the stimvols by that many junkvols
  for iStim = 1:length(stimvol)
    if ~isempty(stimvol{iStim})
      % subtract junkFrames
      stimvol{iStim} = stimvol{iStim}-junkFrames(i)*designSupersampling;
      % and get rid of anything less than 0
      stimDurations{iStim} = stimDurations{iStim}(stimvol{iStim}>0);
      stimvol{iStim} = stimvol{iStim}(stimvol{iStim}>0);
      % check for stimvol overrun
      if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
          runlen = d.dim(end)*designSupersampling;
      else
          runlen = diff(d.concatInfo.runTransition(i,:))*designSupersampling+1;
      end
      if ~isempty(find(stimvol{iStim}>runlen,1))
        if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
           disp(sprintf('(getStimvol) Removing %i event(s) from scan since they happen after the last volume of the scan ',length(find(stimvol{iStim}>runlen))));
        else
           disp(sprintf('(getStimvol) Removing %i event(s) from concatenated scan %i:%s since they happen after the last volume (%i) of the scan ',length(find(stimvol{iStim}>runlen)),i,d.concatInfo.filename{i},runlen));
        end
        % update trialNum and stimvol to remove trials that occur after the end
        if isfield(d,'trialNum'),d.trialNum{i}{iStim} = d.trialNum{i}{iStim}(stimvol{iStim}<=runlen);end
        stimDurations{iStim} = stimDurations{iStim}(stimvol{iStim}<=runlen);
        stimvol{iStim} = stimvol{iStim}(stimvol{iStim}<=runlen);
      end
    end
  end

  if i==1
    d.stimNames = stimNames;
    d.stimvol = stimvol;
    d.stimDurations = stimDurations;
  else
    % if we have more than one stimfile, than we have to concatenate
    % together the stimvol. For this we are going to need to have
    % a concatInfo field
    % check for valid concatInfo
    if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
      disp(sprintf('(getStimvol) No concatInfo found for multiple stimfiles'));
      return
    end
    
    %I dont' think we should require all files to have exactly the same stim names
    %we'll just put stims with identical names together
    oldStimvol = d.stimvol;
    oldStimDurations = d.stimDurations;
    oldStimNames = d.stimNames;
    %put names together (this also sorts names in alphabetical order)
    d.stimNames = union(d.stimNames,stimNames);
    %create new cell arrays
    d.stimvol = cell(size(d.stimNames));
    d.stimDurations = cell(size(d.stimNames));
    %find the indices of the old names in the cell array of all names
    [dump,whichStims] = ismember(oldStimNames,d.stimNames);
    %and put the previous scan stim info in th right place
    d.stimDurations(whichStims)=oldStimDurations;
    d.stimvol(whichStims)=oldStimvol;
    %now find the indices of the new names in the cell array of all names
    [dump,whichStims] = ismember(stimNames,d.stimNames);
    for iStim = 1:length(whichStims)
      thisStimIndex = whichStims(iStim);
      d.stimvol{thisStimIndex}= [d.stimvol{thisStimIndex} ...
        stimvol{iStim}+(d.concatInfo.runTransition(i,1)-1)*designSupersampling];
      d.stimDurations{thisStimIndex}=[d.stimDurations{thisStimIndex} stimDurations{iStim}];
    end
  end

%   if (i > 1)
%     % now add the number of volumes encountered from the previous runs
%     % to the stimvol and concatenate on to general stimvol
%     for iStim = 1:length(stimvol)
%       if length(stimvol) >= iStim
% 	% if we already have some stimvols for this stimulus type then concatenate
% 	d.stimvol{iStim} = [d.stimvol{iStim} (stimvol{iStim}+(d.concatInfo.runTransition(i,1)-1)*designSupersampling)];
%       else
% 	% first time we have encoutered this stimvol just add it to d.stimvol
% 	d.stimvol{iStim} = (stimvol{iStim}+(d.concatInfo.runTransition(i,1)-1)*designSupersampling);
%       end
%   d.stimDurations{iStim} = [d.stimDurations{iStim} stimDurations{iStim}];
%     end
%     % on first file, we just set stimvol in the d field
%   else
%     d.stimvol = stimvol;
%     d.stimDurations = stimDurations;
%   end
end

% return the actual TR supersampling factor
d.designSupersampling = designSupersampling;

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
function [stimOnsets,stimDurations] = getStimvolFromEventTimes(stimfile,tr,duration)

nhdr = length(stimfile.stimtimes_s);
stimOnsets = cell(1, nhdr);
stimDurations = cell(1, nhdr);

%just round to time resolution
for i = 1:nhdr
  stimOnsets{i} = round(stimfile.stimtimes_s{i}/tr)+1;
  if ischar(duration) && strcmp(duration,'fromFile')
    stimDurations{i} = max(1,round(stimfile.stimdurations_s{i}/tr));
  else
    stimDurations{i} =  max(ones(size(stimOnsets{i})),round(duration/tr));
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
