% glmAnalysis.m
%
%      usage: thisView = glmAnalysis(thisView,params)
%         by: farshad moradi,  modified by julien besle 12/01/2010 to 25/10/2010 to perform an ANOVAs(statistic-tests) and statistic-tests
%       date: 06/14/07
%    purpose: GLM analysis using design matrix convolved with any HRF model (including deconvolution)
%              $Id$
%
% 
%     Fits a GLM to the data using ordinary or generalized Least Squares
%     GLM is composed of EVs which can be stimulus times or combinations thereof
%     outputs r2 overlay
%     computes any number of contrasts
%     statistic-tests are performed on any linear combinations of Explanatory Variables (EVs) against zero (H0 = contrast'*beta
%     statistic-tests tests any set of contrasts, not necessarily identical to above contrasts to subsets of  (H0 = no contrast differs from zero)

function [thisView,params] = glmAnalysis(thisView,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help glmAnalysis
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end
if ieNotDefined('params'),params = [];end

% First get parameters
if isempty(params) || justGetParams
  params = glmAnalysisGUI('thisView',thisView,'params',params,'defaultParams',defaultParams,'scanList',scanList);
end

% Abort if params empty
if ieNotDefined('params')
  disp('(glmAnalysis) GLM analysis cancelled');
  return
% just return parameters
elseif justGetParams
  return
end


% set the group
thisView = viewSet(thisView,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);  %this seems quite useless to me
%is it just here to link parameters to file names ?

%just to have shorter variable names
scanParams = params.scanParams;
numberFtests = length(params.restrictions);
numberContrasts = size(params.contrasts,1);
numberTests = numberFtests+numberContrasts;
computePermutations = numberTests && (params.permutationTests || (params.parametricTests && params.permutationFweAdjustment));

if params.covCorrection   %number of voxels to get around the ROI/subset box in case the covariance matrix is estimated
   voxelsMargin = floor(params.covEstimationAreaSize/2);
else
   voxelsMargin = 0;
end
if params.spatialSmoothing  %we'll also need a margin if we're spatially smoothing
  voxelsMargin = max(voxelsMargin,params.spatialSmoothing);
end

%--------------------------------------------------------- Main loop over scans ---------------------------------------------------
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
%initialize the data we're keeping for output overlays
precision = mrGetPref('defaultPrecision');
r2 = cell(1,params.scanNum(end));
if numberContrasts
  contrast = r2;
end
if numberTests
  if params.parametricTests
    parametricP = r2;  
    if params.outputStatistic
      statistic = r2;
    end
    fdrParametricP=r2;
    fweParametricP=r2;
    if params.bootstrapFweAdjustment
      bootstrapFweParametricP = r2;
    end
    if params.permutationFweAdjustment
      permuteFweParametricP = r2;
    end
    if params.TFCE 
      if params.bootstrapFweAdjustment
        bootstrapFweTfceP = r2;
      end
      if params.permutationFweAdjustment
        permuteFweTfceP = r2;
      end
    end
  end
  if params.permutationTests
    permuteP = r2;
    fdrPermuteP=r2;
    fwePermuteP=r2;
    if params.TFCE
      tfcePermuteP = r2;
      fdrTfcePermuteP=r2;
      fweTfceRandT=r2;
    end
  end
  if params.bootstrapTests
    bootstrapP = r2;
    fdrBootstrapP=r2;
    fweBootstrapP=r2;
    if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
      bootstrapFweBootstrapP = r2;
    end
    if params.TFCE  
      tfceBootstrapP = r2;
      fdrTfceBootstrapP=r2;
      fweTfceBootstrapP=r2;
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        bootstrapFweTfceBootstrapP = r2;
      end
    end
  end
end


for iScan = params.scanNum
  numVolumes = viewGet(thisView,'nFrames',iScan);
  scanDims = viewGet(thisView,'dims',iScan);

  %compute the dimensions of the subset of voxels to load
  switch(params.analysisVolume)
   case {'Whole volume','Subset box'}
       subsetBox = eval(scanParams{iScan}.subsetBox);
   case {'Loaded ROI(s)','Visible ROI(s)'}
      %get the smallest box containing all the voxels from all the (visible) ROIs 
      if strcmp(params.analysisVolume,'Visible ROI(s)')
        roiList = viewGet(thisView,'visibleRois');
      else
        roiList = 1:viewGet(thisView,'numberOfRois');
      end
      [subsetBox, whichRoi, marginVoxels] = getRoisBox(thisView,iScan,[voxelsMargin voxelsMargin 0],roiList);
      usedVoxelsInBox = marginVoxels | any(whichRoi,4);
      %clear('whichRoi','marginVoxels');
      if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None')
        [dump,brainMaskRoiNum] = ismember(params.covEstimationBrainMask,viewGet(thisView,'roiNames'));
        brainMaskScanCoords = getROICoordinates(thisView,brainMaskRoiNum,iScan);
        %keep only those voxels that are in the subset box
        brainMaskScanCoords(:,any(brainMaskScanCoords-repmat(subsetBox(:,2),1,size(brainMaskScanCoords,2))>0))=[];
        brainMaskScanCoords=brainMaskScanCoords-repmat(subsetBox(:,1),1,size(brainMaskScanCoords,2))+1;
        brainMaskScanCoords(:,any(brainMaskScanCoords<1))=[];
        if ~isempty(brainMaskScanCoords)
          %create a mask the same size as the subsetbox
          covEstimationBrainMask = false(size(usedVoxelsInBox));
          covEstimationBrainMask(sub2ind(size(usedVoxelsInBox),brainMaskScanCoords(1,:)',brainMaskScanCoords(2,:)',brainMaskScanCoords(3,:)'))=true;
        else
          covEstimationBrainMask = [];
        end
      end
  end
  subsetDims = diff(subsetBox,1,2)'+1;
  %Compute the number of slices to load at once in order to minimize memory usage
  if params.TFCE && computePermutations
   %in this particular case, we force loading the whole dataset at once
   slicesToLoad = subsetDims(3);
   loadCallsPerBatch = {1};
   rawNumSlices = subsetDims(3);
   %compare the size of the array to load to the memory preference and ask for confirmation if larger than maxBlockSize
   switch (precision)
     case {'double'}
      bytesPerNum = 8;
     case{'single'}
      bytesPerNum = 4;
   end
   sizeArray = bytesPerNum*numVolumes*prod(subsetDims);
   if sizeArray> mrGetPref('maxBlocksize')
      if ~askuser(sprintf('(glmAnalysis) This will load an array of %.2f Gb in memory. Are you sure you want to proceed ?', sizeArray/1024^3));
        return;
      end
   end
  else
    switch(params.analysisVolume)

      case {'Whole volume','Subset box'}
        [maxNSlices rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,subsetDims,precision);
        maxNSlices = min(subsetDims(3),maxNSlices);
        slicesToLoad = maxNSlices*ones(1,floor(subsetDims(3)/maxNSlices));
        if rem(subsetDims(3),maxNSlices)
          slicesToLoad(end+1) = rem(subsetDims(3),maxNSlices);
        end
        loadCallsPerBatch = num2cell(1:length(slicesToLoad));

      case {'Loaded ROI(s)','Visible ROI(s)'}
        [maxNSlices rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,subsetDims,precision);
        maxNSlices = min(subsetDims(3),maxNSlices);
        slicesToLoad = maxNSlices*ones(1,floor(subsetDims(3)/maxNSlices));
        if rem(subsetDims(3),maxNSlices)
          slicesToLoad(end+1) = rem(subsetDims(3),maxNSlices);
        end
        % since the box contains voxels we won't use, 
        % choose how many loads we can do while still computing all voxels at once 
        %(we assume that we can at least load one entire slice by calls to loadScan)
        loadCallsPerBatch = {[]};
        nVoxels = 0;
        currentFirstSlice = 0;
        nBatches = 1;
        for iLoad = 1:length(slicesToLoad)
          nVoxels =  nVoxels+nnz(usedVoxelsInBox(:,:,currentFirstSlice+(1:slicesToLoad(iLoad))) > 0);
          [roiSlicesAtATime rawRoiSlices] = getNumSlicesAtATime(numVolumes,[nVoxels 1 1]);
          if rawRoiSlices<1
            nBatches = nBatches+1;
            currentFirstSlice = sum(slicesToLoad(1:iLoad));
            loadCallsPerBatch{nBatches} = [];
          end
          loadCallsPerBatch{nBatches} = [loadCallsPerBatch{nBatches} iLoad];
        end
    end
  end

  if rawNumSlices<1
      mrWarnDlg(['Too many data points to perform the analysis on whole slices. Implement row analysis support, increase memory block size or reduce subset box/ROI(s) size (by a factor ' num2str(subsetDims(3)/rawNumSlices) ')']);
      return;
  end

  ehdr = []; s2 = []; autoCorrelationParameters = [];
  ehdrBootstrapCIs = [];  contrastBootstrapCIs = [];


  % calculate which slices we will be working on
  lastSlices = cumsum(slicesToLoad); 
  firstSlices = lastSlices - slicesToLoad + 1; 
  %loop
  for iBatch = 1:length(loadCallsPerBatch)
    switch(params.analysisVolume)
      case {'Whole volume','Subset box'}
        %Load data
        d = loadScan(thisView,iScan,[],subsetBox(3,1) + [firstSlices(iBatch) lastSlices(iBatch)] -1, precision,subsetBox(1,:),subsetBox(2,:));
      case {'Loaded ROI(s)','Visible ROI(s)'}
        for iLoad= loadCallsPerBatch{iBatch}
          dummy = loadScan(thisView,iScan,[],subsetBox(3,1) + [firstSlices(iLoad) lastSlices(iLoad)] -1, precision,subsetBox(1,:),subsetBox(2,:));
          dummy.data = reshape(dummy.data,[prod(dummy.dim(1:3)) dummy.dim(4)]);
          dummy.data = dummy.data(usedVoxelsInBox(:,:,firstSlices(iLoad):lastSlices(iLoad))>0,:);
          if iLoad==1
            d=dummy;
          else
            d.data = cat(1,d.data, dummy.data);
          end
        end
        clear('dummy','usedVoxelsInBox');
        d.data = permute(d.data,[1 3 4 2]);
        d.roiPositionInBox = any(whichRoi(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:),4);
        d.marginVoxels = marginVoxels(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:);
        if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None')&& ~isempty(covEstimationBrainMask)
          d.covEstimationBrainMask = covEstimationBrainMask(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:);
        else
          d.covEstimationBrainMask = [];
        end
    end
    clear('dummy');
    
    if iBatch==1
      %----------------------------------------Design matrix: same for all voxels/runs
      % get the stim volumes, if empty then abort
      d = getStimvol(d,scanParams{iScan});
      if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
      % do any call for preprocessing
      if ~isempty(scanParams{iScan}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
      end

      actualStimvol = d.stimvol;
      %precompute permutation vectors
      if computePermutations
        nResamples = params.nResamples;
        %find which event types will be randomized
        %first which EVs are involved in all contrasts/restrictions
        contrastEVs = any([params.contrasts;cell2mat(params.restrictions)],1);
        %then which event types constitute these EVs
        randEvents = any(params.scanParams{iScan}.stimToEVmatrix(:,contrastEVs),2);
        %now compute permutation indices for those events, while keeping other events indices unchanged
        nEventTypes = length(actualStimvol);
        numberEvents = NaN(nEventTypes,1);
        totalNumberEvents = 0;
        stimvolToPermuteIndices = [];
        for iEventType = 1:nEventTypes
           numberEvents(iEventType) = length(actualStimvol{iEventType});
           if randEvents(iEventType)
             stimvolToPermuteIndices = [stimvolToPermuteIndices totalNumberEvents+(1:numberEvents(iEventType))];
           end
           totalNumberEvents = totalNumberEvents+numberEvents(iEventType);
        end
        stimvolToPermute = cell2mat(actualStimvol(randEvents));
        nStimsToPermute = length(stimvolToPermute);
        permutations = repmat(cell2mat(actualStimvol),nResamples,1);
        for iPerm = 1:nResamples
           permutations(iPerm,stimvolToPermuteIndices) = stimvolToPermute(randperm(nStimsToPermute));
        end
      else
        nResamples = 0;
      end

      %create model HRF
      [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,0,1);

      d.volumes = 1:d.dim(4);
      %make a copy of d
      actualD = d;

    else
      thisD=d;
      d = actualD;
      d.data = thisD.data;
      clear('thisD');
    end
    d.dim = size(d.data);
      
    %-------------------------------Permutations------------------------------------------
    for iPerm = 1:nResamples+1
      if iPerm==2
        hWaitBar = mrWaitBar(0,['(glmAnalysis) Computing Permutations for scan ' num2str(iScan)]);
      end
      if iPerm == 1 %if it is the actual data
        d.stimvol = actualStimvol;
        actualData = 1;
        verbose = 1;
      else
        actualData = 0;
        mrWaitBar(iPerm/nResamples,hWaitBar);
        permutedStimvol = permutations(iPerm-1,:);
        for stimnum = 1:length(actualStimvol)%randomize stim events
           d.stimvol{stimnum} = permutedStimvol(1:numberEvents(stimnum));
           permutedStimvol = permutedStimvol(numberEvents(stimnum)+1:end);
        end
        verbose = 0;
      end

      %compute the design matrix for this permutation
      d = makeDesignMatrix(d,params,verbose, iScan);
      % compute estimates and statistics
      [d, out] = getGlmStatistics(d, params, verbose, precision, actualData);%, computeTtests,computeBootstrap);
       
      
    
      if iPerm==1
        %keep values for permutation tests
        if computePermutations
          actualStatistic = out.statistic;
          if params.permutationTests
            permuteStatisticCount = zeros(size(actualStatistic),precision); %will count how many permutation values are above the actual statistic value 
            permuteStatisticCount(isnan(actualStatistic)) = NaN; 
          end
          if params.parametricTests && params.permutationFweAdjustment
            permuteFweStatisticCount = zeros(size(actualStatistic),precision); %same for the max permutation value across space, for FWE adjustment
            permuteFweStatisticCount(isnan(actualStatistic)) = NaN; 
            permuteFweData = initResampleFWE(reshape(actualStatistic,prod(d.dim(1:3)),numberTests),params,out.parametricP);
          end
        end
        
        %reshape data into the subsetbox
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) 
          d.ehdr = reshapeToRoiBox(d.ehdr,d.roiPositionInBox);
          out.r2 = reshapeToRoiBox(out.r2,d.roiPositionInBox);
          d.s2 = reshapeToRoiBox(d.s2,d.roiPositionInBox);
          if params.covCorrection
            d.autoCorrelationParameters = reshapeToRoiBox(d.autoCorrelationParameters,d.roiPositionInBox);
          end
          if params.bootstrapTests && params.bootstrapIntervals
            d.ehdrBootstrapCIs = reshapeToRoiBox(d.ehdrBootstrapCIs,d.roiPositionInBox);
          end
          if numberTests
            if numberContrasts && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              out.contrast = reshapeToRoiBox(out.contrast,d.roiPositionInBox);
              if params.bootstrapTests && params.bootstrapIntervals
                d.contrastBootstrapCIs = reshapeToRoiBox(d.contrastBootstrapCIs,d.roiPositionInBox);
              end
            end
            if params.parametricTests
              out.parametricP = reshapeToRoiBox(out.parametricP,d.roiPositionInBox);
              if params.outputStatistic
                out.statistic = reshapeToRoiBox(out.statistic,d.roiPositionInBox);
              end
              if params.bootstrapFweAdjustment
                out.bootstrapFweParametricP = reshapeToRoiBox(out.bootstrapFweParametricP,d.roiPositionInBox);
                if params.TFCE 
                  out.bootstrapFweTfceP = reshapeToRoiBox(out.bootstrapFweTfceP,d.roiPositionInBox);
                end
              end
            end
            if params.bootstrapTests
              out.bootstrapP = reshapeToRoiBox(out.bootstrapP,d.roiPositionInBox);
              if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                out.bootstrapFweBootstrapP = reshapeToRoiBox(out.bootstrapFweBootstrapP,d.roiPositionInBox);
              end
              if params.TFCE 
                out.tfceBootstrapP = reshapeToRoiBox(out.tfceBootstrapP,d.roiPositionInBox);
                if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                   out.bootstrapFweTfceBootstrapP = reshapeToRoiBox(out.bootstrapFweTfceBootstrapP,d.roiPositionInBox);
                end
              end
            end
          end
        end
        
        %concatenate data
        ehdr = cat(3,ehdr,d.ehdr);
        r2{iScan} = cat(3,r2{iScan},out.r2);
        s2 = cat(3,s2,d.s2);
        if params.covCorrection
          autoCorrelationParameters = cat(3,autoCorrelationParameters,d.autoCorrelationParameters);
        end
        if params.bootstrapTests && params.bootstrapIntervals
          ehdrBootstrapCIs = cat(3,ehdrBootstrapCIs,d.ehdrBootstrapCIs);
        end
        if numberTests
          if numberContrasts
            if length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
              contrast{iScan} = cat(3,contrast{iScan},out.contrast);
            end
            if params.bootstrapTests && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              contrastBootstrapCIs = cat(3,contrastBootstrapCIs,d.contrastBootstrapCIs);
            end
          end
          if params.parametricTests
            parametricP{iScan} = cat(3,parametricP{iScan},out.parametricP);
            if params.outputStatistic
              statistic{iScan} = cat(3,statistic{iScan},out.statistic);
            end
            if params.bootstrapFweAdjustment
              bootstrapFweParametricP{iScan} = cat(3,bootstrapFweParametricP{iScan},out.bootstrapFweParametricP);
              if params.TFCE 
                bootstrapFweTfceP{iScan} = cat(3,bootstrapFweTfceP{iScan},out.bootstrapFweTfceP);
              end
            end
          end
          if params.bootstrapTests
            bootstrapP{iScan} = cat(3,bootstrapP{iScan},out.bootstrapP);
            if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
              bootstrapFweBootstrapP{iScan} = cat(3,bootstrapFweBootstrapP{iScan},out.bootstrapFweBootstrapP);
            end
            if params.TFCE 
              tfceBootstrapP{iScan} = cat(3,tfceBootstrapP{iScan},out.tfceBootstrapP);
              if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                bootstrapFweTfceBootstrapP{iScan} = cat(3,bootstrapFweTfceBootstrapP{iScan},out.bootstrapFweTfceBootstrapP);
              end
            end
          end
        end
        %We compute the TFCE here only in the case we forced loading the whole volume at once
        if params.TFCE && computePermutations
          thisTfce = applyFslTFCE(out.statistic,'',0);  
          %put NaNs back
          thisTfce(isnan(out.statistic)) = NaN;
          actualTfce = thisTfce;           %keep the TFCE transform
          if params.permutationTests
            permuteTfceCount = zeros(size(out.statistic),precision);     %these will count how many permutation values are above the actual TFCE values voxelwise
            permuteTfceCount(isnan(out.statistic)) = NaN;%put NaNs back
          end
          if params.parametricTests && params.permutationFweAdjustment
            permuteFweTfceCount = zeros(size(out.statistic),precision); %same for permutation adjusted FWE
            permuteFweTfceCount(isnan(out.statistic)) = NaN;%put NaNs back
            permuteFweTfceData = initResampleFWE(reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),params);
          end
          clear('out');
        end
      
      else % iPerm>1
        if params.permutationTests
          permuteStatisticCount = permuteStatisticCount+double(out.statistic>actualStatistic);
        end
        if params.parametricTests && params.permutationFweAdjustment
          permuteFweStatisticCount = permuteFweStatisticCount + reshape(resampleFWE(reshape(out.statistic,prod(d.dim(1:3)),numberTests),...
                                                            reshape(actualStatistic,prod(d.dim(1:3)),numberTests),...
                                                            permuteFweData), [d.dim(1:3) numberTests]);
        end
        if params.TFCE
          if isfield(d,'roiPositionInBox') 
            out.statistic = reshapeToRoiBox(out.statistic,d.roiPositionInBox); 
          end
          %We compute the TFCE here only in the case we forced loading the whole volume at once
          %NaNs in the data will be transformed to 0 by FSL, so for ROIs, 
          %TFCE is applied to the smallest box including all the ROIs, but replacing non-computed data by 0
          thisTfce = applyFslTFCE(out.statistic,'',0);  
          %put NaNs back
          thisTfce(isnan(out.statistic)) = NaN;
          clear('out');
          if params.permutationTests
            permuteTfceCount = permuteTfceCount+ double(thisTfce>actualTfce);
          end
          if params.parametricTests && params.permutationFweAdjustment
            %permuteFweTfceCount = permuteFweTfceCount+ double(repmat(max(max(max(thisTfce,[],3),[],2),[],1),[subsetDims 1])>actualTfce);
            permuteFweTfceCount = permuteFweTfceCount + reshape(resampleFWE(reshape(thisTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                permuteFweTfceData),size(permuteFweTfceCount));
          end
          clear('thisTfce');
        end
      end
    end
    if computePermutations
      clear('actualStatistic','actualTfce')
      mrCloseDlg(hWaitBar);
    end
    d = rmfield(d,'data');

    if computePermutations %compute P-values, reshape and concatenate if needed
      if params.permutationTests
        permuteStatisticCount = computeRandP(permuteStatisticCount,nResamples);
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) %reshape data into the subsetbox
          permuteStatisticCount = reshapeToRoiBox(permuteStatisticCount,d.roiPositionInBox);
        end
        permuteP{iScan} = cat(3,permuteP{iScan},permuteStatisticCount);  
        clear('permuteStatisticCount');
        if params.TFCE 
          tfcePermuteP{iScan} = computeRandP(permuteTfceCount,nResamples);
          clear('permuteTfceCount');
          tfcePermuteP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
        end
      end
      if params.parametricTests && params.permutationFweAdjustment
        permuteFweStatisticCount = computeRandP(permuteFweStatisticCount,nResamples);
        permuteFweStatisticCount(isnan(permuteFweStatisticCount)) = NaN; %put NaNs back in place
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) %reshape data into the subsetbox
          permuteFweStatisticCount = reshapeToRoiBox(permuteFweStatisticCount,d.roiPositionInBox);
        end
        permuteFweParametricP{iScan} = cat(3,permuteFweParametricP{iScan},permuteFweStatisticCount);
        clear('permuteFweStatisticCount');
        if params.TFCE 
          permuteFweTfceP{iScan} = computeRandP(permuteFweTfceCount,nResamples);
          clear('permuteFweTfceCount');
          permuteFweTfceP{iScan} = reshape(enforceMonotonicityResampleFWE(...
                                      reshape(permuteFweTfceP{iScan},numel(permuteFweTfceP{iScan}(:,:,:,1)),numberContrasts+numberFtests),...
                                      permuteFweTfceData), size(permuteFweTfceP{iScan}));
          permuteFweTfceP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
        end
      end
      d = rmfield(d,'roiPositionInBox');
    end
  end
  
  if numberTests
    %make parametric probability maps
    if params.parametricTests
      [parametricP{iScan},fdrParametricP{iScan},fweParametricP{iScan}] = transformStatistic(parametricP{iScan},precision,params); 
      if params.permutationFweAdjustment
        permuteFweParametricP{iScan} = transformStatistic(permuteFweParametricP{iScan},precision,params); 
        if params.TFCE 
          permuteFweTfceP{iScan} = transformStatistic(permuteFweTfceP{iScan},precision,params); 
        end
      end
      if params.bootstrapFweAdjustment
        bootstrapFweParametricP{iScan} = transformStatistic(bootstrapFweParametricP{iScan},precision,params);
        if params.TFCE 
          bootstrapFweTfceP{iScan} = transformStatistic(bootstrapFweTfceP{iScan},precision,params);
        end
      end
    end
    if params.permutationTests
      [permuteP{iScan},fdrPermuteP{iScan},fwePermuteP{iScan}] = transformStatistic(permuteP{iScan},precision,params); 
      if params.TFCE 
        [tfcePermuteP{iScan},fdrTfcePermuteP{iScan},fweTfceRandT{iScan}] = transformStatistic(tfcePermuteP{iScan},precision,params); 
      end
    end
    if params.bootstrapTests
      [bootstrapP{iScan},fdrBootstrapP{iScan},fweBootstrapP{iScan}] = transformStatistic(bootstrapP{iScan},precision,params); 
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        bootstrapFweBootstrapP{iScan} = transformStatistic(bootstrapFweBootstrapP{iScan},precision,params);
      end
      if params.TFCE 
        [tfceBootstrapP{iScan},fdrTfceBootstrapP{iScan},fweTfceBootstrapP{iScan}] = transformStatistic(tfceBootstrapP{iScan},precision,params);
        if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
          bootstrapFweTfceBootstrapP{iScan} = transformStatistic(bootstrapFweTfceBootstrapP{iScan},precision,params);
        end
      end
    end
  end
  

  % save eventRelated parameters
  glmAnal.d{iScan} = copyFields(d);
  stimvol = d.stimvol;
  for i=1:length(stimvol)
    stimvol{i} = unique(ceil(stimvol{i}/d.designSupersampling));
  end
  glmAnal.d{iScan}.stimvol = stimvol;
  glmAnal.d{iScan}.hrf = downsample(d.hrf, d.designSupersampling/scanParams{iScan}.estimationSupersampling);%/d.designSupersampling*scanParams{iScan}.estimationSupersampling;
  glmAnal.d{iScan}.actualhrf = d.hrf;
  clear('d');
  glmAnal.d{iScan}.estimationSupersampling = scanParams{iScan}.estimationSupersampling;
  glmAnal.d{iScan}.acquisitionSubsample = scanParams{iScan}.acquisitionSubsample;
  glmAnal.d{iScan}.EVnames = params.EVnames;                %this should be removed if viewGet can get params from GLM analysis
  glmAnal.d{iScan}.dim = [scanDims numVolumes];
  glmAnal.d{iScan}.ehdr = NaN([scanDims size(ehdr,4) size(ehdr,5)],precision);
  glmAnal.d{iScan}.ehdr(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdr;
  glmAnal.d{iScan}.s2 = NaN(scanDims,precision);
  glmAnal.d{iScan}.s2(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = s2; %JB
  clear('ehdr','s2');
  if params.covCorrection
    glmAnal.d{iScan}.autoCorrelationParameters = NaN([scanDims size(autoCorrelationParameters,4)],precision);
    glmAnal.d{iScan}.autoCorrelationParameters(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = autoCorrelationParameters; %JB
    clear('autoCorrelationParameters')
  end
  if params.bootstrapTests && params.bootstrapIntervals
    glmAnal.d{iScan}.ehdrBootstrapCIs = NaN([scanDims size(ehdrBootstrapCIs,4) size(ehdrBootstrapCIs,5)],precision);
    glmAnal.d{iScan}.ehdrBootstrapCIs(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdrBootstrapCIs; %JB
    clear('ehdrBootstrapCIs')
  end
  if numberContrasts
    glmAnal.d{iScan}.contrasts = params.contrasts;
    if params.bootstrapTests && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
      glmAnal.d{iScan}.contrastBootstrapCIs = NaN([scanDims size(contrastBootstrapCIs,4) size(contrastBootstrapCIs,5)],precision);
      glmAnal.d{iScan}.contrastBootstrapCIs(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = contrastBootstrapCIs; %JB
      clear('contrastBootstrapCIs')
    end
  end
  if numberFtests
    glmAnal.d{iScan}.fTestNames = params.fTestNames;        %this should be removed if viewGet can get params from GLM analysis
    glmAnal.d{iScan}.restrictions = params.restrictions;               %this should be removed if viewGet can get params from GLM analysis
  end


end

tic
%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
glmAnal.name = params.saveName;
if strcmp(params.hrfModel,'hrfDeconvolution')
  glmAnal.type = 'deconvAnal';
elseif isempty(params.restrictions) && ~params.computeTtests
  glmAnal.type = 'glmAnal';
else
  glmAnal.type = 'glmAnalStats';
end
glmAnal.groupName = params.groupName;
glmAnal.function = 'glmAnalysis';
glmAnal.reconcileFunction = 'defaultReconcileParams';
glmAnal.mergeFunction = 'defaultMergeParams';
glmAnal.guiFunction = 'glmAnalysisGUI';
glmAnal.params = params;
glmAnal.date = dateString;

%--------------------------------------------------------- Output overlay structures
nScans = viewGet(thisView,'nScans');
% create generic parameters 
defaultOverlay.groupName = params.groupName;
defaultOverlay.function = 'glmAnalysis';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.date = dateString;
defaultOverlay.params = cell(1,nScans);
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = 'glmPlot';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [0 1];
defaultOverlay.clip = [0 1];
defaultOverlay.alphaOverlay='';
defaultOverlay.alphaOverlayExponent=1;
defaultOverlay.data = cell(1,nScans);
defaultOverlay.name = '';
for iScan = params.scanNum
   defaultOverlay.data{iScan} = NaN(scanDims,precision); %to make values outside the box transparent
end


%------------------------------------------------------ save the r2 overlay
overlays = defaultOverlay;
overlays.name = 'r2';
overlays.colormapType = 'setRangeToMax';
for iScan = params.scanNum
   overlays.data{iScan}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) = r2{iScan};
   overlays.params{iScan} = scanParams{iScan};
end
clear('r2');

%--------------------------------------------- save the contrast beta weights overlay(s) (ehdr if no contrast)
contrastNames = makeContrastNames(params.contrasts,params.EVnames,params.tTestSide);
if numberContrasts && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
  %this is to mask the beta values by the probability/Z maps     
  betaAlphaOverlay = cell(numberContrasts,1);

  %find the values for the scale of the beta overlays
  thisOverlay = defaultOverlay;
  ordered_abs_betas = cell2mat(contrast(params.scanNum));
  ordered_abs_betas = ordered_abs_betas(~isnan(ordered_abs_betas));
  min_beta = min(min(min(min(min(ordered_abs_betas)))));
  max_beta = max(max(max(max(max(ordered_abs_betas)))));
  ordered_abs_betas = sort(abs(ordered_abs_betas));
  beta_perc95 = 0; 
  beta_perc95 = max(beta_perc95,ordered_abs_betas(round(numel(ordered_abs_betas)*.95))); %take the 95th percentile for the min/max
  thisOverlay.range = [-beta_perc95 beta_perc95];
  thisOverlay.clip = [min_beta max_beta];
  thisOverlay.colormap = jet(256);
  for iContrast = 1:numberContrasts
    overlays(end+1)=thisOverlay;
    overlays(end).alphaOverlay=betaAlphaOverlay{iContrast};
    overlays(end).name = contrastNames{iContrast};
    for iScan = params.scanNum
      overlays(end).data{iScan} = NaN(scanDims,precision); %to make values outside the box transparent
      overlays(end).data{iScan}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
         contrast{iScan}(:,:,:,iContrast);
      overlays(end).params{iScan} = scanParams{iScan};
    end
  end
  clear('contrast');
end

%--------------------------------------------- save the Test overlay(s)
if numberTests
  
  testNames = [contrastNames params.fTestNames];
  for iTest = 1:numberContrasts+numberFtests
    if iTest<=numberContrasts
      testNames{iTest} = ['T (' testNames{iTest} ')'];
    else
      testNames{iTest} = ['F (' testNames{iTest} ')'];
    end
  end
  
  if params.outputStatistic
    tempStatistic = cell2mat(statistic);
    if numberContrasts && params.computeTtests
      thisOverlay = defaultOverlay;
      thisOverlay.colormap = jet(256);
      max_abs_T = max(max(max(max(max(abs(tempStatistic(:,:,:,1:numberContrasts)))))));
      thisOverlay.range = [-max_abs_T max_abs_T];
      thisOverlay.clip = thisOverlay.range;
      for iTest = 1:numberContrasts
        overlays(end+1)=thisOverlay;
        overlays(end).name = testNames{iTest};
        for iScan = params.scanNum
          overlays(end).data{iScan}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             statistic{iScan}(:,:,:,iTest);
          overlays(end).params{iScan} = scanParams{iScan};
        end
      end
    end
    if numberFtests
      thisOverlay = defaultOverlay;
      thisOverlay.range(1) = min(min(min(min(min(tempStatistic(:,:,:,numberContrasts+(1:numberFtests)))))));
      thisOverlay.range(2) = min(max(max(max(max(tempStatistic(:,:,:,numberContrasts+(1:numberFtests)))))));
      thisOverlay.clip = thisOverlay.range;
      for iTest = numberContrasts+(1:numberFtests)
        overlays(end+1)=thisOverlay;
        overlays(end).name = testNames{iTest};
        for iScan = params.scanNum
          overlays(end).data{iScan}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             statistic{iScan}(:,:,:,iTest);
          overlays(end).params{iScan} = scanParams{iScan};
        end
      end
    end
    clear('statistic','tempStatistic');
  end

  if params.parametricTests
    overlays = [overlays makeOverlay(defaultOverlay, parametricP, subsetBox, params.scanNum, scanParams, ...
                                       '', params.testOutput, testNames)];
    clear('parametricP');
    if params.bootstrapFweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap-FWE-adjusted ', params.testOutput, testNames)];
      clear('bootstrapFweParametricP');
      if params.TFCE  
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-FWE-adjusted TFCE ', params.testOutput, testNames)];
        clear('bootstrapFweTfceP');
      end
    end
    if params.permutationFweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, permuteFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                         'permutation-FWE-adjusted ', params.testOutput, testNames)];
      clear('permuteFweParametricP');
      if params.TFCE
        overlays = [overlays makeOverlay(defaultOverlay, permuteFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                         'permutation-FWE-adjusted TFCE ',params.testOutput, testNames)];
        clear('permuteFweTfceP');
      end
    end
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrParametricP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted ',params.testOutput, testNames)];
      clear('fdrParametricP');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fweParametricP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted ',params.testOutput, testNames)];
      clear('fweParametricP');
    end
  end

  if params.bootstrapTests
    overlays = [overlays makeOverlay(defaultOverlay, bootstrapP, subsetBox, params.scanNum, scanParams, ...
                                       'bootstrap ', params.testOutput, testNames)];
    clear('bootstrapP');
    if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap-FWE-adjusted bootstrap ', params.testOutput, testNames)];
      clear('bootstrapFweBootstrapP');
    end
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted bootstrap ',params.testOutput, testNames)];
      clear('fdrBootstrapTp');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fweBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted bootstrap ',params.testOutput, testNames)];
      clear('fweBootstrapTp');
    end
    if params.TFCE  
      overlays = [overlays makeOverlay(defaultOverlay, tfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap TFCE ', params.testOutput, testNames)];
      clear('tfceBootstrapP');
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-FWE-adjusted bootstrap TFCE ', params.testOutput, testNames)];
        clear('bootstrapFweTfceBootstrapP');
      end
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'FDR-adjusted bootstrap TFCE ', params.testOutput, testNames)];
        clear('fdrTfceBootstrapP');
      end
      if params.fweAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fweTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'FWE-adjusted bootstrap TFCE ', params.testOutput, testNames)];
        clear('fweTfceBootstrapP');
      end
    end
  end

  if params.permutationTests
    overlays = [overlays makeOverlay(defaultOverlay, permuteP, subsetBox, params.scanNum, scanParams, ...
                                       'permutation ', params.testOutput, testNames)];
    clear('permuteP');
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrPermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted permutation ',params.testOutput, testNames)];
      clear('fdrPermuteP');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fwePermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted permutation ',params.testOutput, testNames)];
      clear('fwePermuteP');
    end
    if params.TFCE
      overlays = [overlays makeOverlay(defaultOverlay, tfcePermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'permutation TFCE ',params.testOutput, testNames)];
      clear('tfcePermuteP');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTfcePermuteP, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted permutation TFCE ',params.testOutput, testNames)];
        clear('fdrTfcePermuteP');
      end
      if params.fweAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fweTfceRandT, subsetBox, params.scanNum, scanParams, ...
                                        'FWE-adjusted permutation TFCE ',params.testOutput, testNames)];
        clear('fweTfceRandT');
      end
    end
  end
  
  if params.computeTtests && params.maskContrastOverlay && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
    %set the contrast alpha overlay to the statistical value
    switch(params.testOutput)
      case 'P value'                                                  %statistical maps
        betaAlphaOverlayExponent = -1;      % with inverse masking for p values
      case {'Z value','-log10(P) value}'}
        betaAlphaOverlayExponent = .5;      %or normal masking for Z or log10(p) values
    end
    for iContrast = 2:numberContrasts+1
      overlays(iContrast).alphaOverlayExponent=betaAlphaOverlayExponent;
      overlays(iContrast).alphaOverlay = overlays(end-numberFtests-numberContrasts+iContrast-1).name;
    end
  end
  
end
    
%-------------------------------------------------------- Set the analysis in view
glmAnal.overlays = overlays;
thisView = viewSet(thisView,'newAnalysis',glmAnal);
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
toc

%-------------------------------------------------------- Save the analysis
saveAnalysis(thisView,glmAnal.name);

oneTimeWarning('nonZeroHrfStart',0);
oneTimeWarning('tfceOutputsZeros',0);
set(viewGet(thisView,'figNum'),'Pointer','arrow');
refreshMLRDisplay(viewGet(thisView,'viewNum'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sub routines

function outputData = reshapeToRoiBox(inputData,dataPosition)
%inputData must be size ([A 1 1 B], where A equals the number of non-zero values in dataPosition
%outputData will be size ([size(dataPosition) B])
outputData = NaN([numel(dataPosition) size(inputData,4) size(inputData,5)]);
inputData = permute(inputData,[1 4 5 2 3]); %remove singleton dimensions (y and z)
outputData(dataPosition>0,:,:) = inputData;
outputData = reshape(outputData,[size(dataPosition,1) size(dataPosition,2) size(dataPosition,3) size(outputData,2) size(outputData,3)]);



function [convertedStatistic, fdrAdjustedStatistic, fweAdjustedStatistic] = transformStatistic(p, outputPrecision, params)

convertedStatistic = convertStatistic(p, params.testOutput, outputPrecision);

if params.fdrAdjustment
  fdrAdjustedP = p;
  for iTest = 1:size(p,4)
    fdrAdjustedP(:,:,:,iTest) = fdrAdjust(p(:,:,:,iTest),params);
  end
  fdrAdjustedStatistic = convertStatistic(fdrAdjustedP, params.testOutput, outputPrecision);
else
  fdrAdjustedStatistic = [];
end

if params.fweAdjustment
  fweAdjustedP = p;
  for iTest = 1:size(p,4)
    if ismember(params.fweMethod,{'Adaptive Step-down','Adaptive Single-step'})
      [numberTrueH0,lambda] = estimateNumberTrueH0(p(:,:,:,iTest),params);
    else
      numberTrueH0 = nnz(~isnan(p(:,:,:,iTest)));
      lambda = 0;
    end
    fweAdjustedP(:,:,:,iTest) = fweAdjust(p(:,:,:,iTest),params,numberTrueH0,lambda);
  end
  fweAdjustedStatistic = convertStatistic(fweAdjustedP, params.testOutput, outputPrecision);
else
  fweAdjustedStatistic = [];
end


function p = computeRandP(count,nResamples)
%we do not allow probabilities of 0 and replace them by minP
%(this can occur if no resampling return a value as high as the actual value)
p = max(count/nResamples,1/(nResamples+1));
p(isnan(count)) = NaN; %NaNs must remain NaNs (they became minP when using max)


function p = convertStatistic(p, outputStatistic, outputPrecision)

switch(outputStatistic)
  case 'Z value'
    %convert P value to Z value, 
    p = norminv(1-p);  
    p(p<0) = 0; % we're not interested in negative Z value
    %if there was no round-off error from cdf, we could do the following:
    %Z = max(-norminv(p),0);  %because the left side of norminv seems to be less sensitive to round-off errors,
    %we get -norminv(x) instead of norminv(1-x). also we'renot interested in negative Z value
  case '-log10(P) value'
    %convert P to -log10(P) 
    p = -log10(p);
    
end
p = eval([outputPrecision '(p)']);



function adjustedPdata = fdrAdjust(pData,params)

isNotNan = ~isnan(pData);
sizePdata = size(pData);
pData = pData(isNotNan);
[pData,sortingIndex] = sort(pData);
numberH0 = length(pData); 

switch(params.fdrAssumption)
  case 'Independence/Positive dependence'
    dependenceCorrectionConstant =1;
  case 'None'
    dependenceCorrectionConstant = sum(1./(1:numberH0));
    %approximate value: 
    %dependenceCorrectionConstant = log(numberH0)+0.57721566490153 %Euler constant
end

%estimate number of True Null hypotheses
% ref: Benjamini et al. 2006, Biometrika, 93(3) p491
switch(params.fdrMethod)
  case 'Step-up'
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1); 
      
  case 'Adaptive Step-up' %definition 3 in Benjamini et al. 2006
    %estimate the number of true H0 using one of the methods specified by params.trueNullsEstimationMethod
    numberTrueH0 = min(ceil(estimateNumberTrueH0(pData,params)),numberH0);
    %recompute the adjusted p-value using the estimated ratio of true H0
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,numberTrueH0/numberH0);
    
  case 'Two-stage Step-up'
    %we will estimate the number of true H0 using definition 6 in Benjamini et al. 2006
    % first adjust the p values with q'=(q+1)
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1);
    numberFalseH0 = length(find(adjustedP<(params.trueNullsEstimationThreshold/(params.trueNullsEstimationThreshold+1))));
    if numberFalseH0 %if any voxel is significant at the chosen level usually .05/(.05+1)
      numberH0 = length(adjustedP); 
      %recompute the adjusted p-value using the estimated number of true H0
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,(numberH0 - numberFalseH0)/numberH0);
    end
    
  case 'Multiple-stage Step-up'
    %we will estimate the number of true H0 using definition 7 in Benjamini et al. 2006
    %but they only give the definition for p-value correction with a fixed q
    q = params.trueNullsEstimationThreshold;
    %step-up version with L>=J (I think it's a step-up procedure within a step-down procedure)
    %k=max{i : for all j<=i there exists l>=j so that p(l)<=ql/{m+1?j(1?q)}}.
    i=0;          %
    l = numberH0; %these are just arbitrary values to pass the first while test
    while i<numberH0 && ~isempty(l)
      i = i+1;
      l = find(pData(i:numberH0)<=q*(i:numberH0)'./(numberH0+1-i*(1-q)),1,'last');
    end
    k=i;
    correctedP = pData;
    correctedP(k:end) = 1;
    % adjustement procedure seems a bit complicated... so I'll leave it for now
    % I don't think the results would be much different from the next, simpler, step-down procedure
    adjustedP = correctedP;
    
  case 'Multiple-stage Step-down'
    %this is the step-down version of definition 7 in Benjamini et al. 2006 with L=J
%     q = params.trueNullsEstimationThreshold;
    %k=max{i : for all j<=i there exists l=j so that p(l)<=ql/{m+1?j(1?q)}}.
    %k=max{i : for all j<=i p(j)<=qj/{m+1?j(1?q)}}.
    %this is just the loop version that corresponds to the definition
% % % %     i=1;
% % % %     while  i<numberH0 && all(pData(1:i)<=q*(1:i)'./(numberH0+1-(1:i)'*(1-q)))
% % % %       i = i+1;
% % % %     end
% % % %     k=i;
    %and this a one line version: 
% % % %     k = find(pData>q*(1:numberH0)'./(numberH0+1-(1:numberH0)'*(1-q)),1,'first');
% % % % 
% % % %     correctedP = pData;
% % % %     correctedP(k:end) = 1;
    
    %but I'd rather have an adjustment version:
    %we need to find for each test i, the largest q such that there exists k=max{i : for all j<=i p(j)<=qj/{m+1?j(1?q)}}
    %this is the derivation for q, starting from the previous one-line p-correction:
%     pData>q*(1:numberH0)'./(numberH0+1-(1:numberH0)'*(1-q))
%     pData.*(numberH0+1-(1:numberH0)'*(1-q))>q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'*(1-q) > q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'+ pData.*(1:numberH0)'*q > q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'  > q*(1:numberH0)' - pData.*(1:numberH0)'*q
%     pData*(numberH0+1) - pData.*(1:numberH0)'  > q*((1:numberH0)' - pData.*(1:numberH0)')
%     (pData*(numberH0+1) - pData.*(1:numberH0)')./((1:numberH0)' - pData.*(1:numberH0)')  > q
    adjustedP = (pData*(numberH0+1) - pData.*(1:numberH0)')./((1:numberH0)' - pData.*(1:numberH0)');
    
    %i think we need to enforce monotonicity, but in contrast with the step-up method (see below),
    %it should be done from the smallest to the largest p-value
    for i=2:numberH0
     adjustedP(i) = max(adjustedP(i-1),adjustedP(i));
    end
    adjustedP = min(adjustedP,1); %bound with 1
          
end

adjustedP(sortingIndex) = adjustedP;
adjustedPdata = NaN(sizePdata);
adjustedPdata(isNotNan) = adjustedP;



function qData = linearStepUpFdr(pData,dependenceCorrectionConstant,trueH0Ratio)

numberH0 = length(pData);
%adjustment (does not depend on a pre-chosen threshold)
% it consists in finding, for each p-value, 
% the largest FDR threshold such that the p-value would be considered significant
% it is the inverse operation to p-correction (see commented code below)
% ref: Yekutieli and Benjamini 1999, Journal of Statistical Planning and Inference, 82(1-2) p171
qData = min(pData.*numberH0./(1:numberH0)'.*dependenceCorrectionConstant*trueH0Ratio,1);

%Now, there may be situations where the largest q for a given uncorrected p-value
% is smaller than the largest q for a higher uncorrected p-value 
%(e.g. consecutive identical uncorrected p-values)
% which means that when looking at the map of corrected p-value (=q-values),
% for a given q-value you would reject less hypothesis than in the p-correction method...

%So maybe we should make corrected p-values monotonic in order to have the same number of rejections for the same threshold
%(although i can't find a paper that says that)
%make it monotonic from the end (for each i, q(i) must be smaller than all q between i and n) 
for i=numberH0-1:-1:1
 qData(i) = min(qData(i),qData(i+1));
end

% p-correction method
% % % %correction (p-correction depends on a fixed threshold, here .05)
% % % qAlpha = .05*(1:numberH0)'/numberH0/dependenceCorrectionConstant;
% % % correctedP = pData;
% % % correctedP(find(pData<=qAlpha,1,'last'):end) = 1;
% % % correctedP(sortingIndex) = correctedP;
% % % correctedPdata = NaN(size(pData));
% % % correctedPdata(isNotNan) = correctedP;


function overlays = makeOverlay(overlays,data,subsetBox,scanList,scanParams,testTypeString, testOutput, testNames)
  switch testOutput
    case 'P value'
      overlays.colormap = statsColorMap(256);
      namePrefix = [testTypeString 'P ['];
    case 'Z value'
      overlays.range = [0 8.2096];
      namePrefix = [testTypeString 'Z ['];
    case '-log10(P) value'
      overlays.range = [0 16];
      namePrefix = ['-log10(' testTypeString 'P) ['];
  end
  overlays.clip = overlays.range;
  overlays = repmat(overlays,1,length(testNames));
  for iTest = 1:length(testNames)
    overlays(iTest).name = [namePrefix testNames{iTest} ']'];
    for iScan = scanList
      overlays(iTest).data{iScan}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
         data{iScan}(:,:,:,iTest);
      overlays(iTest).params{iScan} = scanParams{iScan};
    end
  end



