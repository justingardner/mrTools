% makescm.m
%
%      usage: makescm(d,hdrlen)
%         by: justin gardner
%       date: 07/28/04
%       e.g.: makescm(d)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it. if d.hdrlen is set then you can
%             just pass in the data structure
%
function d = makescm(d,hdrlen)

if (nargin == 1)
  if (isfield(d,'hdrlen'))
    hdrlen = d.hdrlen;
  else
    hdrlen = 25;
  end
elseif (nargin ~= 2)
  help makescm;
  return
end

% if we have only a single run then we set
% the runTransitions for that single run
if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
  runTransition = [1 d.dim(4)];
else
  runTransition = d.concatInfo.runTransition;
end

% go through each run of the experiment
allscm = [];
for runnum = 1:size(runTransition,1)
  % default values
  scm = [];
  % make stimulus convolution matrix
  for stimnum = 1:length(d.stimvol)
    % make an array containing the stimulus times
    stimarray = zeros(1,runTransition(runnum,2)-runTransition(runnum,1)+1);
    % only use stimvols that are within this runs volume numbers
    stimarray(d.stimvol{stimnum}(find((d.stimvol{stimnum}>=runTransition(runnum,1)) & (d.stimvol{stimnum}<=runTransition(runnum,2))))-runTransition(runnum,1)+1) = 1;
    % stack stimcmatrices horizontally
    scm = [scm stimconv(stimarray,hdrlen)];
  end
  % stack this run's stimcmatrix on to the last one
  allscm = [allscm;scm];
end

% set values
d.nhdr = length(d.stimvol);
d.scm = allscm;
d.hdrlen = hdrlen;
d.volumes = 1:d.dim(4);

