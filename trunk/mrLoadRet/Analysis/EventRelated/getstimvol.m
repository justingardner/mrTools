% getstimvol.m
%
%      usage: d = getstimvol(d)
%         by: justin gardner
%       date: 12/21/05
%    purpose: gets the stimulus vols from the traces
%
function d = getstimvol(d,keepstim)

% check arguments
if ~any(nargin == [1 2])
  help getstimvol
  return
end
if exist('keepstim')~=1,keepstim = [];,end;

% get the stim times
stimraw = d.channels(d.stimchannel,:);
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

% get the first volume
d = getfirstvol(d);

% check if we need to only keep some vols
if ~isempty(keepstim)
  if (min(keepstim) < 1) | (max(keepstim) > length(d.stimvol))
    disp(sprintf('UHOH: Keepstim out of range. Ignoring'));
  else
    % go through and only keep the stim values asked for
    stimvol = [];
    for i = 1:length(keepstim)
      stimvol{i} = d.stimvol{keepstim(i)};
    end
    d.stimvol = stimvol;
  end
end

  
