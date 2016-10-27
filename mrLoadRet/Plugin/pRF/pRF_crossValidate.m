%%  pRF_crossValidate.m
%%
%%        $Id:$
%%      usage: pRF_crossValidate(view, roiName, analysis)
%%         by: akshay jagadeesh
%%       date: 10/04/16
%%    purpose: Given n MotionComp runs, this function computes a n-fold
%%             cross validation, running the pRF analysis n times.
%%             It generates and saves time series, model response, residuals,
%%             and covariance matrices for each fold.
%%
%%
%%     input:  roiName  - Name of the ROI we want to run this analysis on
%%             analysis - Name of analysis file to load
%%                      - This must be the Concat of Average of N scans

function fit = pRF_crossValidate(roiName, analysis)

%roiName = 'lV1';
%roiName = 'xVal';
%roiName = 'lineV1';
analysis = 'pRF_lV1_123456.mat';

% Set current group to Concat and load the Analysis file
v = newView;
v = viewSet(v, 'currentGroup', 'Concatenation');
v = loadAnalysis(v, ['pRFAnal/' analysis]);
analParams = v.analyses{1}.params;
analScanNum = analParams.scanNum;

ogn = viewGet(v, 'originalgroupname', analScanNum); % Get the average group scan
osn = viewGet(v, 'originalscannum', analScanNum);
ogn2 = viewGet(v, 'originalgroupname', osn, ogn{1}); % Get the motioncomp scans
osn2 = viewGet(v, 'originalscannum', osn, ogn{1});

% Get analysis params for pRFFit
concatInfo = viewGet(v, 'concatinfo', analScanNum);
d = viewGet(v, 'd', analScanNum);

% Get the N original scans & their tSeries
scans = loadROITSeries(v, roiName, osn2, ogn2{1}, 'straightXform=1');

for i=1:1
  
  % Initialize vars for use later
  avg = zeros(size(scans{1}.tSeries));
  leftOut = scans{i}.tSeries;
  numVoxels = size(leftOut, 1);
  modelResponse = zeros(numVoxels, size(leftOut, 2));
  residual = zeros(numVoxels, size(leftOut, 2));
  
  % compute average time series across N-1 scans
  for j = 1:length(scans)
    if i ~= j
      avg = avg + scans{j}.tSeries;
    end
  end
  avg = avg / (length(scans)-1);

  % Apply Concat Filtering to averages & left out
  for k = 1:numVoxels
    filteredAvg(k, :) = applyConcatFiltering(avg(k, :), concatInfo, 1);
    leftOut(k, :) = applyConcatFiltering(leftOut(k, :), concatInfo, 1);
  end

  disppercent(-inf, '(crossVal) Fitting pRF to voxels'); 
  % run pRFFit on the averaged tSeries
  coords = scans{i}.scanCoords.';
  for h = 1:numVoxels
    x = coords(h, 1); y = coords(h, 2); z = coords(h, 3);
    modelFit = pRFFit(v,[],x,y,z, 'tSeries', filteredAvg(h, :).', 'stim', d.stim, 'getModelResponse=1', 'params', d.params(:, h), 'concatInfo', d.concatInfo, 'fitTypeParams', analParams.pRFFit, 'paramsInfo', d.paramsInfo);
    modelResponse(h, :) = modelFit.modelResponse;
    residual(h, :) = modelFit.tSeries - modelFit.modelResponse;

    disppercent(h/numVoxels);
  end
  disppercent(inf);

  %Calculate covariance matrix of voxels on our calculated residual
  covMat = residual*residual';

  % Calculate covariance matrix of voxels using pRFNoise's method
  %[residual, covMat, tSeries, modelResponse] = pRFNoise(v, analysisScanNum, coordinates, analysis);
  %[resid, covMat2, tS, modResp] = pRFNoise([], analScanNum, coords, analysis);
  %noise.resid = resid;
  %noise.covMat = covMat2;
  %noise.tSeries = tS;
  %noise.modelResponse = modResp;

  % Compare model response to left-out timeseries
  probTable = testPRF(leftOut, modelResponse, diag(diag(covMat)));
end
fit.modelResponse = modelResponse;
fit.residual = residual;
fit.leftOut = leftOut;
fit.covMat = covMat;
fit.probTable = probTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applyConcatFiltering    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tSeries = applyConcatFiltering(tSeries,concatInfo,runnum)

% apply the same filter as original data
% check for what filtering was done
tSeries = tSeries(:);

% apply detrending (either if concatInfo does not say what it did or if
% the filterType field has detrend in it)
if ~isfield(concatInfo,'filterType') || ~isempty(findstr('detrend',lower(concatInfo.filterType)))
  tSeries = eventRelatedDetrend(tSeries);
end

% apply hipass filter
if isfield(concatInfo,'hipassfilter') && ~isempty(concatInfo.hipassfilter{runnum})
  % check for length match
  if ~isequal(length(tSeries),length(concatInfo.hipassfilter{runnum}))
    disp(sprintf('(pRFFit:applyConcatFiltering) Mismatch dimensions of tSeries (length: %i) and concat filter (length: %i)',length(tSeries),length(concatInfo.hipassfilter{runnum})));
  else
    tSeries = real(ifft(fft(tSeries) .* repmat(concatInfo.hipassfilter{runnum}', 1, size(tSeries,2)) ));
  end
end

% project out the mean vector
if isfield(concatInfo,'projection') && ~isempty(concatInfo.projection{runnum})
  projectionWeight = concatInfo.projection{runnum}.sourceMeanVector * tSeries;
  tSeries = tSeries - concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
end

% now remove mean
tSeries = tSeries-repmat(mean(tSeries,1),size(tSeries,1),1);

% make back into the right dimensions
tSeries = tSeries(:)';
