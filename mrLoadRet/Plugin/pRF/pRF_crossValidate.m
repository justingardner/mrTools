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
%%             newCoords - in the form [x1 y1 z1 1; x2 y2 z2 1].'

function fits = pRF_crossValidate(newCoords, roiName, analysis)

%% Plot Figures flag (set to [] if you don't want to plot)
plotFigs = 1;

roiName = 'V1';
%analysis = 'pRF_lV1_123456.mat';
analysis = 'pRF.mat';

% Set current group to Concat and load the Analysis file
v = newView;
v = viewSet(v, 'currentGroup', 'Concatenation');
%v = viewSet(v, 'currentGroup', 'Averages');
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
scanDims = viewGet(v, 'scanDims');

% Get the N original scans & their tSeries
%keyboard
if ~ieNotDefined('newCoords')
  if(strmatch(newCoords, 'best'))
    disp(sprintf('(crossValidate) Getting 90 best defined voxels in ROI %s', roiName));
    coords = getBestVoxels(analScanNum, 85, roiName, analysis);
    tempRoi = makeEmptyROI(v, sprintf('scanNum=%i',osn2(1)), 'groupNum', ogn2{1});
    tempRoi.coords = coords;
  else
    disp(sprintf('(crossValidate) Getting tSeries for provided coordinates'));
    tempRoi = makeEmptyROI(v, sprintf('scanNum=%i',osn2(1)), 'groupNum', ogn2{1});
    tempRoi.coords = newCoords;
  end
  scans = loadROITSeries(v, tempRoi, osn2, ogn2{1}, 'straightXform=1');
elseif ~ieNotDefined('roiName')
  disp(sprintf('(crossValidate) Coordinates not given, getting tSeries for ROI %s', roiName));
  scans = loadROITSeries(v, roiName, osn2, ogn2{1}, 'straightXform=1');
else
  disp(sprintf('(crossValidate) Neither ROI name nor coordinates provided. Exiting.'));
  return
end

numFolds = length(scans);
for i=1:numFolds;
  disp(sprintf('(crossValidate) Fold %d of %d', i, numFolds));
  % Initialize vars for use later
  fit = [];
  avg = zeros(size(scans{1}.tSeries));
  unfilteredLeftOut = scans{i}.tSeries;
  numVoxels = size(unfilteredLeftOut, 1);
  modelResponse = zeros(numVoxels, size(unfilteredLeftOut, 2));
  residual = zeros(numVoxels, size(unfilteredLeftOut, 2));
  
  % compute average time series across N-1 scans
  for j = 1:length(scans)
    if i ~= j
      avg = avg + scans{j}.tSeries;
    end
  end
  avg = avg / (length(scans)-1);

  % Apply Concat Filtering to averages & left out
  %keyboard
  for k = 1:numVoxels
    filteredAvg(k, :) = applyConcatFiltering(avg(k, :), concatInfo, 1);
    leftOut(k, :) = applyConcatFiltering(unfilteredLeftOut(k, :), concatInfo, 1);
  end

  disppercent(-inf, sprintf('\t(crossVal) Fitting pRF to %d voxels', numVoxels)); 
  % run pRFFit on the averaged tSeries
  coords = scans{i}.scanCoords.';
  for h = 1:numVoxels
    x = coords(h, 1); y = coords(h, 2); z = coords(h, 3);
    linCoords = sub2ind(scanDims, x, y, z);
    linVox = find(d.linearCoords == linCoords);
    if isempty(linVox)
      disp(sprintf('Encountered voxel, (%d, %d, %d), not included in d.linearCoords', x, y, z));
      continue 
    end
    modelFit = pRFFit(v,[],x,y,z, 'tSeries', filteredAvg(h, :).', 'stim', d.stim, 'getModelResponse=1', 'params', d.params(:, linVox), 'fitTypeParams', analParams.pRFFit, 'paramsInfo', d.paramsInfo, 'concatInfo', concatInfo);
    modelResponse(h, :) = modelFit.modelResponse;
    residual(h, :) = modelFit.tSeries - modelFit.modelResponse;
    fitParams(h, :) = [x, y, z, modelFit.p.x, modelFit.p.y, modelFit.p.std];
    disppercent(h/numVoxels);
  end
  leftOut(all(modelResponse==0, 2), :) = [];
  modelResponse(all(modelResponse==0, 2), :) = [];
  residual(all(residual==0,2), :) = [];
  numSuccess = size(modelResponse, 1);
  disp(sprintf('(crossVal) Successfully fit model to %d voxels', numSuccess)); 
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
  fit.modelResponse = modelResponse;
  fit.residual = residual;
  fit.leftOut = leftOut;
  fit.covMat = covMat;
  fit.probTable = probTable;
  fit.fitParams = fitParams;

  eval(sprintf('fits(%d)=fit;', i))
  timeLen = size(fit.modelResponse, 2);
  if ~ieNotDefined('plotFigs')
    figure; subplot(2, 1, 1); imagesc(log(fit.probTable)); axis ij; colorbar; title(sprintf('Fold %d of %d', i, numFolds));
    subplot(2, 1, 2); plot((1:timeLen), fit.leftOut, 'k');
    hold on; plot((1:timeLen), fit.modelResponse, 'r'); xlim([1 timeLen]); title('Left out time series + Model response');
  end

end

% Average together probability tables
pTable = zeros(timeLen, timeLen);
modResp = zeros(numSuccess, timeLen);
resid = zeros(numSuccess, timeLen);
leftOut = zeros(numSuccess, timeLen);
covMat = zeros(numSuccess, numSuccess);
fp = zeros(
for i = 1:length(fits)
  fit = fits(i);
  pTable = pTable + fit.probTable;
  modResp = modResp + fit.modelResponse;
  leftOut = leftOut + fit.leftOut;
  covMat = covMat + fit.covMat;
  resid = resid + fit.residual;
end
fit = [];
fit.modelResponse = modResp / length(fits);
fit.residual = resid / length(fits);
fit.leftOut = leftOut / length(fits);
fit.covMat = covMat / length(fits);
fit.probTable = pTable / length(fits);
fit.fitParams = [];
eval(sprintf('fits(%d)=fit;', numFolds+1));
%if ~ieNotDefined('plotFigs')
  %fit = fits(numFolds+1);
  %plotFigs(fit, d);
%end
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
  %disp(sprintf('(pRFFit:applyConcatFiltering) Apply detrending'));
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
