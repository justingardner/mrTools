% makescm.m
%
%      usage: makeglm(d,hrf)
%         by: farshad moradi
%       date: 06/14/07
%       e.g.: makeglm(d, hrf)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it. if d.hrf is set then you can
%             just pass in the data structure
%
function d = makeglm(d,hrf)

if (nargin == 1)
  if (isfield(d,'hrf'))
      hrf = d.hrf;
  else
      help makeglm;
      return
  end
elseif (nargin ~= 2)
  help makeglm;
  return
end

% make sure hrf starts from zero
hrf = hrf-repmat(hrf(1,:), size(hrf,1), 1);
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
    stimarray = zeros(1,(runTransition(runnum,2)-runTransition(runnum,1))*d.supersampling+1);
    % only use stimvols that are within this runs volume numbers
    stimarray(d.stimvol{stimnum}(find((d.stimvol{stimnum}>=runTransition(runnum,1)*d.supersampling) & ...
        (d.stimvol{stimnum}<=runTransition(runnum,2)*d.supersampling)))-runTransition(runnum,1)*d.supersampling+1) = 1;
    m = convn(stimarray', hrf);
    m = m(1:length(stimarray),:);
    % remove mean 
    m = m-repmat(mean(m), size(m,1), 1);
    % downsample
    m = cumsum(m,1);
    mid_tr = floor(d.supersampling/2)+1;
    m = diff([zeros(1, size(m,2));m(mid_tr:d.supersampling:end, :)]);
    % stack stimcmatrices horizontally
    scm = [scm, m];
  end
  % stack this run's stimcmatrix on to the last one
  allscm = [allscm;scm];
end

% set values
d.nhdr = length(d.stimvol);
d.scm = allscm;
d.hdrlen = size(hrf,2);
d.volumes = 1:d.dim(4);

