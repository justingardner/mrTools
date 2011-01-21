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

function thisView = glmAnalysis(thisView,params,varargin)

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
if ieNotDefined('scanList')
   scanList = [];
end
if ieNotDefined('params')
   % First get parameters
   params = glmAnalysisGUI('thisView',thisView,'defaultParams',defaultParams,'scanList',scanList);
end

% Abort if params empty
if ieNotDefined('params')
  disp('(glmAnalysis) GLM analysis cancelled');
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

if params.covCorrection   %number of voxels to get around the ROI/subset box in case the covariance matrix is estimated
   voxelsMargin = floor(params.covEstimationAreaSize/2);
else
   voxelsMargin = 0;
end

% % % % if any( sum(contrasts,2)==1 & sum(contrasts~=0,2)==1 )
% % % %    scanNum = params.scanNum(1);
% % % %    fprintf('Getting condition names from scan %d', scanNum); 
% % % %    d = loadScan(thisView, scanNum, [], 0);
% % % %    d = getStimvol(d,scanParams{iScan});
% % % %    % do any call for preprocessing
% % % %    if ~isempty(scanParams{iScan}.preprocess)
% % % %       d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
% % % %    end
% % % % end


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
    if params.TFCE && params.permutationTests
      permuteFweTfceP = r2;
    end
  end
  if params.permutationTests
    permuteP = r2;
    permuteFweParametricP = r2;
    fdrPermuteP=r2;
    fwePermuteP=r2;
    if params.TFCE
      tfcePermuteP = r2;
      fdrTfcePermuteP=r2;
      fweTfceRandT=r2;
    end
  end
  if params.bootstrapStatistics
    bootstrapP = r2;
    fdrBootstrapP=r2;
    fweBootstrapP=r2;
    if ~strcmp(params.resampleFWEadjustment,'None')
      bootstrapFweParametricP = r2;
      bootstrapFweBootstrapP = r2;
    end
    if params.TFCE  
      tfceBootstrapP = r2;
      fdrTfceBootstrapP=r2;
      fweTfceBootstrapP=r2;
      if ~strcmp(params.resampleFWEadjustment,'None')
        bootstrapFweTfceP = r2;
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
   case 'Loaded ROI(s)'
      %get the smallest box containing all the voxels from all the boxes 
      [subsetBox, whichRoi, marginVoxels] = getRoisBox(thisView,iScan,[voxelsMargin voxelsMargin 0]);
      usedVoxelsInBox = marginVoxels | any(whichRoi,4);
      %clear('whichRoi','marginVoxels');
  end
  subsetDims = diff(subsetBox,1,2)'+1;
  %Compute the number of slices to load at once in order to minimize memory usage
  if params.TFCE && (params.permutationTests || params.bootstrapStatistics)
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

      case 'Loaded ROI(s)'
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
      case 'Loaded ROI(s)'
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
    end
    clear('dummy');
    
    if iBatch==1
      %----------------------------------------Design matrix: same for all voxels/runs
      % get the stim volumes, if empty then abort
      d = getStimvol(d,scanParams{iScan});
      if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
      % do any call for preprocessing
      if ~isempty(scanParams{iScan}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess,verbose);
      end

      actualStimvol = d.stimvol;
      %precompute permutation vectors
      if params.permutationTests
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
       
      
      if params.permutationTests
        if isfield(out,'statistic');
          thisStatistic = out.statistic;
        end
      end
      
      if iPerm==1
        %keep values for permutation tests
        if params.permutationTests
          actualStatistic = thisStatistic;
          permuteStatisticCount = zeros(size(actualStatistic),precision); %will count how many permutation values are above the actual statistic value 
          permuteStatisticCount(isnan(thisStatistic)) = NaN; 
          permuteFweStatisticCount = zeros(size(actualStatistic),precision); %same for the max permutation value across space, for FWE adjustment
          %this is for resampled-based FWE control
          permuteFweData = initResampleFWE(reshape(actualStatistic,prod(d.dim(1:3)),numberTests),params);
        end
        
        %reshape data into the subsetbox
        if strcmp(params.analysisVolume,'Loaded ROI(s)') 
          d.ehdr = reshapeToRoiBox(d.ehdr,d.roiPositionInBox);
          out.r2 = reshapeToRoiBox(out.r2,d.roiPositionInBox);
          d.s2 = reshapeToRoiBox(d.s2,d.roiPositionInBox);
          if params.covCorrection
            d.autoCorrelationParameters = reshapeToRoiBox(d.autoCorrelationParameters,d.roiPositionInBox);
          end
          if params.bootstrapStatistics && params.bootstrapIntervals
            d.ehdrBootstrapCIs = reshapeToRoiBox(d.ehdrBootstrapCIs,d.roiPositionInBox);
          end
          if numberTests
            if numberContrasts && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              out.contrast = reshapeToRoiBox(out.contrast,d.roiPositionInBox);
              if params.bootstrapStatistics && params.bootstrapIntervals
                d.contrastBootstrapCIs = reshapeToRoiBox(d.contrastBootstrapCIs,d.roiPositionInBox);
              end
            end
            if params.parametricTests
              out.parametricP = reshapeToRoiBox(out.parametricP,d.roiPositionInBox);
              if params.outputStatistic
                out.statistic = reshapeToRoiBox(out.statistic,d.roiPositionInBox);
              end
            end
            if params.bootstrapStatistics
              out.bootstrapP = reshapeToRoiBox(out.bootstrapP,d.roiPositionInBox);
              if ~strcmp(params.resampleFWEadjustment,'None')
                out.bootstrapFweParametricP = reshapeToRoiBox(out.bootstrapFweParametricP,d.roiPositionInBox);
                out.bootstrapFweBootstrapP = reshapeToRoiBox(out.bootstrapFweBootstrapP,d.roiPositionInBox);
              end
              if params.TFCE 
                out.tfceBootstrapP = reshapeToRoiBox(out.tfceBootstrapP,d.roiPositionInBox);
                if ~strcmp(params.resampleFWEadjustment,'None')
                  out.bootstrapFweTfceP = reshapeToRoiBox(out.bootstrapFweTfceP,d.roiPositionInBox);
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
        if params.bootstrapStatistics && params.bootstrapIntervals
          ehdrBootstrapCIs = cat(3,ehdrBootstrapCIs,d.ehdrBootstrapCIs);
        end
        if numberTests
          if numberContrasts
            if length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
              contrast{iScan} = cat(3,contrast{iScan},out.contrast);
            end
            if params.bootstrapStatistics && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              contrastBootstrapCIs = cat(3,contrastBootstrapCIs,d.contrastBootstrapCIs);
            end
          end
          if params.parametricTests
            parametricP{iScan} = cat(3,parametricP{iScan},out.parametricP);
            if params.outputStatistic
              statistic{iScan} = cat(3,statistic{iScan},out.statistic);
            end
          end
          if params.bootstrapStatistics
            bootstrapP{iScan} = cat(3,bootstrapP{iScan},out.bootstrapP);
            if ~strcmp(params.resampleFWEadjustment,'None')
              bootstrapFweParametricP{iScan} = cat(3,bootstrapFweParametricP{iScan},out.bootstrapFweParametricP);
              bootstrapFweBootstrapP{iScan} = cat(3,bootstrapFweBootstrapP{iScan},out.bootstrapFweBootstrapP);
            end
            if params.TFCE 
              tfceBootstrapP{iScan} = cat(3,tfceBootstrapP{iScan},out.tfceBootstrapP);
              if ~strcmp(params.resampleFWEadjustment,'None')
                bootstrapFweTfceP{iScan} = cat(3,bootstrapFweTfceP{iScan},out.bootstrapFweTfceP);
                bootstrapFweTfceBootstrapP{iScan} = cat(3,bootstrapFweTfceBootstrapP{iScan},out.bootstrapFweTfceBootstrapP);
              end
            end
          end
        end
        
        clear('out')
      
      % iPerm>1  
      else
        if numberTests
          permuteStatisticCount = permuteStatisticCount+double(thisStatistic>actualStatistic);
          %permuteFweStatisticCount = permuteFweStatisticCount+double(repmat(max(max(max(thisStatistic,[],3),[],2),[],1),[d.dim(1:3) 1])>actualStatistic);
          if ~strcmp(params.resampleFWEadjustment,'None')
            permuteFweStatisticCount = permuteFweStatisticCount + reshape(resampleFWE(reshape(thisStatistic,prod(d.dim(1:3)),numberTests),...
                                                              reshape(actualStatistic,prod(d.dim(1:3)),numberTests),...
                                                              permuteFweData,params),...
                                                  [d.dim(1:3) numberTests]);
          end
        end
      end

      %We compute the TFCE here only in the case we forced loading the whole volume at once
      if params.TFCE && params.permutationTests
        if numberTests
          if isfield(d,'roiPositionInBox') 
            thisStatistic = reshapeToRoiBox(thisStatistic,d.roiPositionInBox); 
          end
          %NaNs in the data will be transformed to 0 by FSL, so for ROIs, 
          %TFCE is applied to the smallest box including all the ROIs, but replacing non-computed data by 0
          thisTfce = applyFslTFCE(thisStatistic,'',0);  
          %put NaNs back
          thisTfce(isnan(thisStatistic)) = NaN;
          if iPerm == 1 %if it is the actual data
            permuteTfceCount = zeros(size(thisStatistic),precision);     %these will count how many permutation values are above the actual TFCE values voxelwise
            %put NaNs back
            permuteTfceCount(isnan(thisStatistic)) = NaN;
            actualTfce = thisTfce;           %keep the TFCE transform
            if ~strcmp(params.resampleFWEadjustment,'None')
              permuteFweTfceCount = permuteTfceCount; %same for permutation adjusted FWE
              permuteFweTfceData = initResampleFWE(reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),params);
            end
          else
            permuteTfceCount = permuteTfceCount+ double(thisTfce>actualTfce);
            if ~strcmp(params.resampleFWEadjustment,'None')
              %permuteFweTfceCount = permuteFweTfceCount+ double(repmat(max(max(max(thisTfce,[],3),[],2),[],1),[subsetDims 1])>actualTfce);
              permuteFweTfceCount = permuteFweTfceCount + reshape(resampleFWE(reshape(thisTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                  reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                  permuteFweTfceData,params),...
                                                      size(permuteFweTfceCount));
            end
          end
        end
        
        clear('thisTfce');
      end
    end
    if params.permutationTests 
       mrCloseDlg(hWaitBar);
    end
    d = rmfield(d,'data');
    clear('actualTfce')

    if params.permutationTests %compute P-values, reshape and concatenate if needed
      if numberTests 
        permuteStatisticCount = computeRandP(permuteStatisticCount,nResamples);
        if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
            permuteStatisticCount = reshapeToRoiBox(permuteStatisticCount,d.roiPositionInBox);
        end
        permuteP{iScan} = cat(3,permuteP{iScan},permuteStatisticCount);          
        if ~strcmp(params.resampleFWEadjustment,'None')
          permuteFweStatisticCount = computeRandP(permuteFweStatisticCount,nResamples);
          permuteFweStatisticCount(isnan(permuteFweStatisticCount)) = NaN; %put NaNs back in place
          if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
            permuteFweStatisticCount = reshapeToRoiBox(permuteFweStatisticCount,d.roiPositionInBox);
          end
          permuteFweParametricP{iScan} = cat(3,permuteFweParametricP{iScan},permuteFweStatisticCount);
        end
        if params.TFCE 
          tfcePermuteP{iScan} = computeRandP(permuteTfceCount,nResamples);
          clear('permuteTfceCount');
          tfcePermuteP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
          if ~strcmp(params.resampleFWEadjustment,'None')
            permuteFweTfceP{iScan} = computeRandP(permuteFweTfceCount,nResamples);
            clear('permuteFweTfceCount');
            permuteFweTfceP{iScan} = reshape(enforceMonotonicityResampleFWE(...
                                        reshape(permuteFweTfceP{iScan},numel(permuteFweTfceP{iScan}(:,:,:,1)),numberContrasts+numberFtests),...
                                        permuteFweTfceData,params),...
                                        size(permuteFweTfceP{iScan}));
            permuteFweTfceP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
          end
        end
      end
      d = rmfield(d,'roiPositionInBox');
    end
    

  end
  clear('permuteStatisticCount','actualStatistic')
  
  if numberTests
    %make parametric probability maps
    if params.parametricTests
      [parametricP{iScan},fdrParametricP{iScan},fweParametricP{iScan}] = transformStatistic(parametricP{iScan},precision,params); 
    end
    if params.permutationTests
      [permuteP{iScan},fdrPermuteP{iScan},fwePermuteP{iScan}] = transformStatistic(permuteP{iScan},precision,params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        permuteFweParametricP{iScan} = transformStatistic(permuteFweParametricP{iScan},precision,params); 
      end
      if params.TFCE 
        [tfcePermuteP{iScan},fdrTfcePermuteP{iScan},fweTfceRandT{iScan}] = transformStatistic(tfcePermuteP{iScan},precision,params); 
        if ~strcmp(params.resampleFWEadjustment,'None')
          permuteFweTfceP{iScan} = transformStatistic(permuteFweTfceP{iScan},precision,params); 
        end
      end
    end
    if params.bootstrapStatistics
      [bootstrapP{iScan},fdrBootstrapP{iScan},fweBootstrapP{iScan}] = transformStatistic(bootstrapP{iScan},precision,params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        bootstrapFweParametricP{iScan} = transformStatistic(bootstrapFweParametricP{iScan},precision,params);
        bootstrapFweBootstrapP{iScan} = transformStatistic(bootstrapFweBootstrapP{iScan},precision,params);
      end
      if params.TFCE 
        [tfceBootstrapP{iScan},fdrTfceBootstrapP{iScan},fweTfceBootstrapP{iScan}] = transformStatistic(tfceBootstrapP{iScan},precision,params);
        if ~strcmp(params.resampleFWEadjustment,'None')
          bootstrapFweTfceP{iScan} = transformStatistic(bootstrapFweTfceP{iScan},precision,params);
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
  glmAnal.d{iScan}.hrf = downsample(d.hrf, d.designSupersampling/scanParams{iScan}.estimationSupersampling)/d.designSupersampling*scanParams{iScan}.estimationSupersampling;
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
  if params.bootstrapStatistics && params.bootstrapIntervals
    glmAnal.d{iScan}.ehdrBootstrapCIs = NaN([scanDims size(ehdrBootstrapCIs,4) size(ehdrBootstrapCIs,5)],precision);
    glmAnal.d{iScan}.ehdrBootstrapCIs(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdrBootstrapCIs; %JB
    clear('ehdrBootstrapCIs')
  end
  if numberContrasts
    glmAnal.d{iScan}.contrasts = params.contrasts;
    if params.bootstrapStatistics && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
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

%--------------------------------------------- save the Test overlay(s)

if numberTests
  
  contrastNames = makeContrastNames(params.contrasts,params.EVnames,params.tTestSide);
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

    if ~strcmp(params.resampleFWEadjustment,'None')
      if params.bootstrapStatistics
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-adjusted ', params.testOutput, testNames)];
        clear('bootstrapFweParametricP');
        if params.TFCE  
          if ~strcmp(params.resampleFWEadjustment,'None')
            overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                               'bootstrap-adjusted TFCE ', params.testOutput, testNames)];
            clear('bootstrapFweTfceP');
          end
        end
      end 
      if params.permutationTests
        overlays = [overlays makeOverlay(defaultOverlay, permuteFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                           'permutation-adjusted ', params.testOutput, testNames)];
        clear('permuteFweParametricP');
        if params.TFCE
          overlays = [overlays makeOverlay(defaultOverlay, permuteFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                           'permutation-adjusted TFCE ',params.testOutput, testNames)];
          clear('permuteFweTfceP');
        end
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

  if params.bootstrapStatistics
    overlays = [overlays makeOverlay(defaultOverlay, bootstrapP, subsetBox, params.scanNum, scanParams, ...
                                       'bootstrap ', params.testOutput, testNames)];
    clear('bootstrapP');
    if ~strcmp(params.resampleFWEadjustment,'None')
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap-adjusted bootstrap ', params.testOutput, testNames)];
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
      if ~strcmp(params.resampleFWEadjustment,'None')
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-adjusted bootstrap TFCE ', params.testOutput, testNames)];
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
end
    
%--------------------------------------------- save the contrast beta weights overlay(s) (ehdr if no contrast)
if numberContrasts && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
  %this is to mask the beta values by the probability/Z maps     
  betaAlphaOverlay = cell(numberContrasts,1);
  if params.computeTtests 
    %set the contrast alpha overlay to the statistical value
    switch(params.testOutput)
      case 'P value'                                                  %statistical maps
        betaAlphaOverlayExponent = -1;      % with inverse masking for p values
      case {'Z value','-log10(P) value}'}
        betaAlphaOverlayExponent = .5;      %or normal masking for Z or log10(p) values
    end
    for iContrast = 1:numberContrasts
      betaAlphaOverlay{iContrast} = overlays(end-numberFtests-numberContrasts+iContrast).name;
    end
  else
    betaAlphaOverlayExponent = 1;
  end

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
    overlays(end).alphaOverlayExponent=betaAlphaOverlayExponent;
  end
  clear('contrast');
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

[convertedStatistic,p] = convertStatistic(p, params.testOutput, outputPrecision);

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
    if ~strcmp(params.adaptiveFweMethod,'None')
      [trueH0ratio,lambda] = estimateTrueH0Ratio(p(:,:,:,iTest),params);
    else
      trueH0ratio = 1;
      lambda = 0;
    end
    fweAdjustedP(:,:,:,iTest) = fweAdjust(p(:,:,:,iTest),params,trueH0ratio,lambda);
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


function [convertedP,p] = convertStatistic(p, outputStatistic, outputPrecision)

switch(outputStatistic)
  case 'Z value'
    %convert P value to Z value, 
    convertedP = norminv(1-p);  
    convertedP(convertedP<0) = 0; % we're not interested in negative Z value
    %if there was no round-off error from cdf, we could do the following:
    %Z = max(-norminv(p),0);  %because the left side of norminv seems to be less sensitive to round-off errors,
    %we get -norminv(x) instead of norminv(1-x). also we'renot interested in negative Z value
  case '-log10(P) value'
    %convert P to -log10(P) 
    convertedP = -log10(p);
    
end
convertedP = eval([outputPrecision '(convertedP)']);



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
  case 'Adaptive Step-up'
    %we will estimate the number of true H0 using definition 3 in Benjamini et al. 2006
    % first adjust the p values 
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1);
    if any(adjustedP<params.adaptiveFdrThreshold) %if any voxel is significant at the chosen level (usually .05) (most likely)
      %compute the estimated number of true H0 when considering the number of rejections at all possible levels
      numberH0 = length(pData);
      numberTrueH0 = (numberH0+1-(1:numberH0)')./(1-pData);
      %find the first estimated number that is greater than the previous one
      numberTrueH0 = ceil(min(numberTrueH0(find(diff(numberTrueH0)>0,1,'first')+1),numberH0));
      %recompute the adjusted p-value using the estimated ratio of true H0
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,numberTrueH0/numberH0);
    end
  case 'Two-stage Adaptive Step-up'
    %we will estimate the number of true H0 using definition 6 in Benjamini et al. 2006
    % first adjust the p values with q'=(q+1)
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1);
    numberFalseH0 = length(find(adjustedP<(params.adaptiveFdrThreshold/(params.adaptiveFdrThreshold+1))));
    if numberFalseH0 %if any voxel is significant at the chosen level usually .05/(.05+1)
      numberH0 = length(adjustedP); 
      %recompute the adjusted p-value using the estimated number of true H0
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,(numberH0 - numberFalseH0)/numberH0);
    end
    
  case 'Multiple-stage Adaptive Step-up'
    %we will estimate the number of true H0 using definition 7 in Benjamini et al. 2006
    %but they only give the definition for p-value correction with a fixed q
    q = params.adaptiveFdrThreshold;
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
    
  case 'Adaptive Step-down'
    %this is the step-down version of definition 7 in Benjamini et al. 2006 with L=J
%     q = params.adaptiveFdrThreshold;
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



