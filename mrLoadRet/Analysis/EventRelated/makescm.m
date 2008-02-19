% makescm.m
%
%      usage: makescm(d,hdrlen)
%         by: justin gardner
%       date: 07/28/04
%       e.g.: makescm(d)
%    purpose: makes a stimulation convolution matrix
%             for data series. this correctly handles
%             run boundarys. it uses the stimulus volumes
%             found in d.stimvol. if hdrlen is not passed in
%             then it uses d.hdrlen (if it exists).
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

if isfield(d, 'hipassfilter')
    hipassfilter = d.hipassfilter;
else
    hipassfilter = [];
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
    m= stimconv(stimarray,hdrlen);
    % apply the same filter as original data
    if ~isempty(hipassfilter)
        m = real(ifft(fft(m) .* repmat(hipassfilter{runnum}', 1, size(m,2)) ));
        m = m-repmat(mean(m,1),size(m,1),1);
    end

    scm = [scm, m];
  end
  % stack this run's stimcmatrix on to the last one
  
  allscm = [allscm; scm];
end

% set values
d.nhdr = length(d.stimvol);
d.scm = allscm;
d.hdrlen = hdrlen;
d.volumes = 1:d.dim(4);

