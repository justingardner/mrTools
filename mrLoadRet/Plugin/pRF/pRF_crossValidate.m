%%  pRF_crossValidate.m
%%
%%      usage: fits = pRF_crossValidate(newCoords, [], analysis) --> run on specified coordinates
%%             fits = pRF_crossValidate([], roiName, analysis) --> run on specified ROI
%%             fits = pRF_crossValidate('best', roiName, analysis) --> run on best N voxels in specified ROI
%%         by: akshay jagadeesh
%%       date: 10/04/16
%%    purpose: Given n MotionComp runs, this function computes a n-fold cross validation, running the pRF 
%%             analysis n times. It generates and saves time series, model response, residuals, and covariance 
%%             matrices for each fold.
%%             We then calculate the probability that each point in the observed time series is predicted
%%             by the pRF model response.
%%
%%
%%     input:  roiName  - Name of the ROI we want to run this analysis on
%%             analysis - Name of analysis file to load
%%                      - This must be the Concat of Average of N scans
%%             newCoords - array of scan coordinates in the form [x1 y1 z1 1; x2 y2 z2 1]'

function fits = pRF_crossValidate(newCoords, roiName, analysis)

%%%%%%% Default inputs %%%%%%%%%
newCoords = 'best';
roiName = 'bothV1';
analysis = 'pRF_v1.mat';
nBest = 5;
plotFigs = []; % set to [] to turn off plots, set to 1 to turn on plots


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
scanDims = viewGet(v, 'scanDims');

% Get the N original scans & their tSeries
if ~ieNotDefined('newCoords')
  tempRoi = makeEmptyROI(v, sprintf('scanNum=%i',osn2(1)), 'groupNum', ogn2{1});
  if(strmatch(newCoords, 'best'))
    disp(sprintf('(crossValidate) Getting %d best defined voxels in ROI %s', nBest, roiName));
    coords = getBestVoxels(v, analScanNum, nBest, roiName, analysis);
    tempRoi.coords = coords;
  else
    disp(sprintf('(crossValidate) Getting tSeries for provided coordinates'));
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
    modelFit = pRFFit(v,[],x,y,z, 'tSeries', filteredAvg(h, :).', 'stim', d.stim, 'getModelResponse=1',...
                      'params', d.params(:, linVox), 'fitTypeParams', analParams.pRFFit, 'paramsInfo', ...
                      d.paramsInfo, 'concatInfo', concatInfo);
    modelResponse(h, :) = modelFit.modelResponse;
    residual(h, :) = modelFit.tSeries - modelFit.modelResponse;
    fitParams(h, :) = [x, y, z, modelFit.p.x, modelFit.p.y, modelFit.p.std];
    disppercent(h/numVoxels);
  end
  leftOut(all(modelResponse==0, 2), :) = [];
  modelResponse(all(modelResponse==0, 2), :) = [];
  residual(all(residual==0,2), :) = [];
  numSuccess = size(modelResponse, 1);
  disppercent(inf);
  disp(sprintf('\t(crossVal) Successfully fit model to %d voxels', numSuccess));

  %Calculate covariance matrix of voxels on our calculated residual
  covMat = residual*residual';

  % Compare model response to left-out timeseries
  probTable = calcProb(leftOut, modelResponse, diag(diag(covMat)));
  fit.modelResponse = modelResponse;
  fit.residual = residual;
  fit.leftOut = leftOut;
  fit.covMat = covMat;
  fit.probTable = probTable;
  fit.fitParams = fitParams;

  eval(sprintf('fits(%d)=fit;', i))
  timeLen = size(fit.modelResponse, 2);
  if ~ieNotDefined('plotFigs')
    figure; 
    subplot(2, 1, 1); imagesc(log(fit.probTable)); axis ij; colorbar; title(sprintf('Fold %d of %d', i, numFolds));
    subplot(2, 1, 2); plot((1:timeLen), fit.leftOut, 'k'); hold on; 
    plot((1:timeLen), fit.modelResponse, 'r'); xlim([1 timeLen]); title('Left out time series + Model response');
  end

end

% Average across folds and store in fits struct array
fit = [];
fit.modelResponse = meanStruct(fits, 'modelResponse');
fit.residual = meanStruct(fits, 'residual');
fit.leftOut = meanStruct(fits, 'leftOut');
fit.covMat = meanStruct(fits, 'covMat');
fit.probTable = meanStruct(fits, 'probTable');
fit.fitParams = meanStruct(fits, 'fitParams');
fits(numFolds+1) = fit;


              %  -- end main program --  %


%%%%%%%%%%%%%%%%%%%%%% Helper Methods %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         meanStruct         %
%                            %
function avg = meanStruct(structArr, field, notI)

if(ieNotDefined('notI'))
  notI = -1;
else
  disp(sprintf('Averaging across all but the %d th fold', notI));
end
eval(sprintf('avg = zeros(size(structArr(1).%s));', field));

for i = 1:length(structArr)
  if i~=notI
    avg = avg + eval(sprintf('structArr(i).%s', field));
  end
end
avg = avg / length(structArr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applyConcatFiltering    %
%                            %
function tSeries = applyConcatFiltering(tSeries,concatInfo,runnum)

tSeries = tSeries(:);

% apply detrending 
if ~isfield(concatInfo,'filterType') || ~isempty(findstr('detrend',lower(concatInfo.filterType)))
  tSeries = eventRelatedDetrend(tSeries);
end

% apply hipass filter
if isfield(concatInfo,'hipassfilter') && ~isempty(concatInfo.hipassfilter{runnum})
  if ~isequal(length(tSeries),length(concatInfo.hipassfilter{runnum}))
    disp(sprintf('(applyConcatFiltering) Mismatch dimensions of tSeries (length: %i) and concat filter (length: %i)',length(tSeries),length(concatInfo.hipassfilter{runnum})));
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
tSeries = tSeries(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        getBestVoxels        %
%                             %
function cords = getBestVoxels(v, scanNum, bestN, roiName, analysis)

groupNum = viewGet(v, 'currentGroup');
roi = loadROITSeries(v, roiName, scanNum, groupNum, 'straightXform=1', 'loadType=none');

r2 = viewGet(v, 'overlaydata', scanNum, 1, 1);
roi = getSortIndex(v, roi, r2);

% Get sort index based on r2 and list of linear coords
sortIndex = roi{1}.sortindex;
scanLinCoords = roi{1}.scanLinearCoords;
scanDims = viewGet(v, 'scanDims');
% get the linear coords with highest sortIndex and convert to 3d coords
[i,j,k] = ind2sub(scanDims, scanLinCoords(sortIndex(1:bestN)));
cords = [i;j;k;ones(1, bestN)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         calcProb            %

function probTable = calcProb(testTSeries, modelResponse, covarMat)

numTimePoints = size(modelResponse, 2);
probTable = zeros(numTimePoints, numTimePoints);

disppercent(-inf, sprintf('\t(calcProb) Calculating prediction likelihoods'));
for i = 1:numTimePoints
  model_i = modelResponse(:, i);
  tSeries_i = testTSeries(:, i);

  for j = 1:numTimePoints
    tSeries_j = testTSeries(:, j);
    probTable(j, i) = mvnpdf(tSeries_j, model_i, covarMat);
    if probTable(j, i) == 0
      disp(sprintf('\n(calcProb) mvnpdf at index(%i, %i) returning 0; Exiting.', j, i));
      return;
    end
  end
  disppercent(i/numTimePoints);
end
disppercent(inf);
