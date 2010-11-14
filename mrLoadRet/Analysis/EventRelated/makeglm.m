% makeglm.m
%
%        $Id$
%      usage: makeglm(d,params,verbose)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 11/02/2010
%       e.g.: makeglm(d,params,verbose)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it, as well as a model hrf
%              optional parameters can be passed in the structure params (params
%
function d = makeglm(d,params,verbose)

if ~any(nargin == [1 2 3])
   help makeglm;
   return
end

if ieNotDefined('params')
  params=struct;
end
if ieNotDefined('verbose')
  verbose = 1;
end

% check if the hrf starts from zero (except if it is the identity matrix, the deconvolution case)
if verbose && d.hrf(1,1)>1e-6 && (size(d.hrf,1)~=size(d.hrf,2) || ~isempty(find(d.hrf^-1-d.hrf>1e-6, 1)))
   mrWarnDlg(['(makeglm) HRF does not start from zero (hrf(0) = ' num2str(d.hrf(1,1)) ')']);
end

if isfield(params,'nonLinearityCorrection') && params.nonLinearityCorrection && isfield(params.hrfParams,'maxModelHrf')
   saturationThreshold = params.saturationThreshold*params.hrfParams.maxModelHrf;
else
   saturationThreshold = repmat(Inf,1,size(d.hrf,2));
end

if isfield(params,'testParams') && isfield(params.testParams,'stimToEVmatrix') && ~isempty(params.testParams.stimToEVmatrix)
  stimToEVmatrix = params.testParams.stimToEVmatrix;
  if size(stimToEVmatrix,1)~=length(d.stimvol)
    mrErrorDlg('(makeglm) EV combination matrix is incompatible with number of event types');
  end
else
  stimToEVmatrix = eye(length(d.stimvol));
end

% if we have only a single run then we set
% the runTransitions for that single run
if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
  runTransition = [1 d.dim(4)];
else
  runTransition = d.concatInfo.runTransition;
end

runTransition(:,1) = ((runTransition(:,1)-1)*round(d.supersampling)+1);
runTransition(:,2) = runTransition(:,2)*round(d.supersampling);

sample_value = 1/d.supersampling;   %we need to correct the amplitude of the hrf if supersampling>1

% go through each run of the experiment
allscm = [];
for runnum = 1:size(runTransition,1)
  scm = [];
   % make an array containing the stimulus times
   stimArray = zeros(runTransition(runnum,2)-runTransition(runnum,1)+1,length(d.stimvol));
   for stimnum = 1:length(d.stimvol)
      % only use stimvols that are within this runs volume numbers
      stimvol = d.stimvol{stimnum};
      stimArray( stimvol(stimvol>=runTransition(runnum,1) & stimvol<=runTransition(runnum,2) )- runTransition(runnum,1)+1,stimnum ) = sample_value;
   end
   % apply EV combination matrix
   stimArray = stimArray*stimToEVmatrix;
   
   % make stimulus convolution matrix
   for iEV = 1:size(stimArray,2)
      m = convn(stimArray(:,iEV), d.hrf);
      m = m(1:size(stimArray,1),:);
      %apply saturation
      m = min(m,repmat(saturationThreshold,size(stimArray,1),1));
      % remove mean 
      m = m-repmat(mean(m), size(m,1), 1); %DOES IT CHANGE ANYTHING IF I REMOVE THIS ?
      % downsample
      m = downsample(m, d.supersampling);
      % apply the same filter as original data
      if isfield(d,'concatInfo') 
         % apply hipass filter
         if isfield(d.concatInfo,'hipassfilter') && ~isempty(d.concatInfo.hipassfilter{runnum})
           m = real(ifft(fft(m) .* repmat(d.concatInfo.hipassfilter{runnum}', 1, size(m,2)) ));
         end
         % project out the mean vector
         if isfield(d.concatInfo,'projection') && ~isempty(d.concatInfo.projection{runnum})
           projectionWeight = d.concatInfo.projection{runnum}.sourceMeanVector * m;
           m = m - d.concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
         end
      end
      % stack stimmatrices horizontally
      scm =  [scm m];
   end
   % stack this run's stimcmatrix on to the last one
   allscm = [allscm;scm];
end

% set values
d.nhdr = size(stimToEVmatrix,2);
d.scm = allscm;
d.hdrlen = size(d.hrf,2);

