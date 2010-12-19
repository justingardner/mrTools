% makeDesignMatrix.m
%
%        $Id: makeDesignMatrix.m 1928 2010-12-14 23:22:11Z julien $
%      usage: makeDesignMatrix(d,params,verbose)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 11/02/2010
%       e.g.: makeDesignMatrix(d,params,verbose)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it, as well as a model hrf
%              optional parameters can be passed in the params structure 
%
function d = makeDesignMatrix(d,params,verbose, scanNum)

if ~any(nargin == [1 2 3 4 5])
   help makeDesignMatrix;
   return
end

if ieNotDefined('params')
  params=struct;
end
if ieNotDefined('verbose')
  verbose = 1;
end

% check if the hrf starts from zero (except if it is the identity matrix, the deconvolution case)
if verbose && any(d.hrf(1,:)>1e-6) && (size(d.hrf,1)~=size(d.hrf,2) || ~isempty(find(d.hrf^-1-d.hrf>1e-6, 1)))
   mrWarnDlg(['(makeDesignMatrix) HRF does not start from zero (hrf(0) = ' mat2str(d.hrf(1,:)) ')']);
end

if isfield(params,'nonLinearityCorrection') && params.nonLinearityCorrection && isfield(params.hrfParams,'maxModelHrf')
   saturationThreshold = params.saturationThreshold*params.hrfParams.maxModelHrf;
else
   saturationThreshold = Inf(1,size(d.hrf,2));
end

if isfield(params.scanParams{scanNum},'estimationSupersampling') && ~isempty(params.scanParams{scanNum}.estimationSupersampling)
  estimationSupersampling = params.scanParams{scanNum}.estimationSupersampling;
else
  estimationSupersampling = 1;
end
if isfield(params.scanParams{scanNum},'acquisitionSubsample') && ~isempty(params.scanParams{scanNum}.acquisitionSubsample)
  acquisitionSubsample = params.scanParams{scanNum}.acquisitionSubsample;
else
  acquisitionSubsample = 1;
end
if isfield(params.scanParams{scanNum},'stimToEVmatrix') && ~isempty(params.scanParams{scanNum}.stimToEVmatrix)
  stimToEVmatrix = params.scanParams{scanNum}.stimToEVmatrix;
  if size(stimToEVmatrix,1)~=length(d.stimvol)
    mrWarnDlg('(makeDesignMatrix) EV combination matrix is incompatible with number of event types');
    d.scm = [];
    return;
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

runTransition(:,1) = ((runTransition(:,1)-1)*round(d.designSupersampling)+1);
runTransition(:,2) = runTransition(:,2)*round(d.designSupersampling);

%apply duration and convert to matrix form
stimMatrix = stimCell2Mat(d.stimvol,d.stimDurations,runTransition);
%if design sampling is larger than estimation sampling, we need to correct the amplitude of the hrf 
stimMatrix = stimMatrix*estimationSupersampling/d.designSupersampling; 

% apply EV combination matrix
d.EVmatrix = stimMatrix*stimToEVmatrix;
allscm = [];
for iRun = 1:size(runTransition,1)
  scm = [];
  thisEVmatrix = d.EVmatrix(runTransition(iRun,1):runTransition(iRun,2),:);
  % make stimulus convolution matrix
  for iEV = 1:size(thisEVmatrix,2)
      m = convn(thisEVmatrix(:,iEV), d.hrf);
      m = m(1:size(thisEVmatrix,1),:);
      %apply saturation
      m = min(m,repmat(saturationThreshold,size(thisEVmatrix,1),1));
      % remove mean 
      m = m-repmat(mean(m), size(m,1), 1); %DOES IT CHANGE ANYTHING IF I REMOVE THIS ?
      % downsample with constant integral to estimation sampling rate
      m = downsample(m, d.designSupersampling/estimationSupersampling);
      %only keep acquisition samples
      m = m(acquisitionSubsample:estimationSupersampling:end,:);
      % apply the same filter as original data
      if isfield(d,'concatInfo') 
         % apply hipass filter
         if isfield(d.concatInfo,'hipassfilter') && ~isempty(d.concatInfo.hipassfilter{iRun})
           m = real(ifft(fft(m) .* repmat(d.concatInfo.hipassfilter{iRun}', 1, size(m,2)) ));
         end
         % project out the mean vector
         if isfield(d.concatInfo,'projection') && ~isempty(d.concatInfo.projection{iRun})
           projectionWeight = d.concatInfo.projection{iRun}.sourceMeanVector * m;
           m = m - d.concatInfo.projection{iRun}.sourceMeanVector'*projectionWeight;
         end
      end
      % stack stimmatrices horizontally
      scm =  [scm m];
   end
   % stack this run's stimcmatrix on to the last one
   allscm = [allscm;scm];
%    d.EVmatrix = [d.EVmatrix;thisEVmatrix];
end

% set values
d.nhdr = size(stimToEVmatrix,2);
d.scm = allscm;
d.hdrlen = size(d.hrf,1);
d.nHrfComponents = size(d.hrf,2);
d.runTransitions = runTransition;
