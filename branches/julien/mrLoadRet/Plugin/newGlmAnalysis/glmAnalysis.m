% glmAnalysis.m
%
%      usage: thisView = glmAnalysis(thisView,params)
%         by: farshad moradi,  modified by julien besle 12/01/2010 to 25/10/2010 to perform an ANOVAs(F-tests) and T-tests
%       date: 06/14/07
%    purpose: GLM analysis using design matrix convolved with any HRF model (including deconvolution)
%              $Id$
%
% 
%     Fits a GLM to the data using ordinary or generalized Least Squares
%     GLM is composed of EVs which can be stimulus times or combinations thereof
%     outputs r2 overlay
%     computes any number of contrasts
%     T-tests are performed on any linear combinations of Explanatory Variables (EVs) against zero (H0 = contrast'*beta
%     F-tests tests any set of contrasts, not necessarily identical to above contrasts to subsets of  (H0 = no contrast differs from zero)

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
restrictions = params.restrictions;
numberFtests = length(restrictions);
contrasts = params.contrasts;
numberContrasts = size(contrasts,1);

if params.covCorrection   %number of voxels to get around the ROI/subset box in case the covariance matrix is estimated
   voxelsMargin = floor(params.covEstimationAreaSize/2);
else
   voxelsMargin = 0;
end

if any( sum(contrasts,2)==1 & sum(contrasts~=0,2)==1 )
   scanNum = params.scanNum(1);
   fprintf('Getting condition names from scan %d', scanNum); 
   d = loadScan(thisView, scanNum, [], 0);
   d = getStimvol(d,scanParams{scanNum});
   % do any call for preprocessing
   if ~isempty(scanParams{scanNum}.preprocess)
      d = eventRelatedPreProcess(d,scanParams{scanNum}.preprocess);
   end
end


%--------------------------------------------------------- Main loop over scans ---------------------------------------------------
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
%initialize the data we're keeping for output overlays
precision = mrGetPref('defaultPrecision');
r2 = cell(1,params.scanNum(end));
if ~isempty(contrasts)
  C = r2;
  if params.computeTtests
    if params.parametricTests
      T = r2;
      Tp = r2;  
      if params.fdrAdjustment
        fdrTp=r2;
      end
      if params.TFCE
        %if params.randomizationTests
          randMaxTfceT = r2;
        %end
      end
    end
    if params.randomizationTests
      randT = r2;
      randMaxT = r2;
      if params.fdrAdjustment
        fdrRandC=r2;
      end
      if params.TFCE
        randTfceT = r2;
        if params.fdrAdjustment
          fdrTfceRandT=r2;
        end
      end
    end
    if params.bootstrapStatistics
      bootstrapT = r2;
      if params.fdrAdjustment
        fdrBootstrapT=r2;
      end
      if ~strcmp(params.resampleFWEadjustment,'None')
        bootstrapMaxT = r2;
        if params.TFCE  
          bootstrapTfceT = r2;
          if params.fdrAdjustment
            fdrBootstrapTfceT=r2;
          end
          if ~strcmp(params.resampleFWEadjustment,'None')
            bootstrapMaxTfceT = r2;
          end
        end
      end
    end
  end
end
if ~isempty(restrictions)
  F = r2;  
  if params.parametricTests
    Fp = r2;
    if params.fdrAdjustment
      fdrFp=r2;
    end
    if params.TFCE
      %if params.randomizationTests
        randMaxTfceF = r2;
      %end
    end
  end
  if params.randomizationTests
    randF = r2;
    randMaxF = r2;
    if params.fdrAdjustment
      fdrRandF=r2;
    end
    if params.TFCE
      randTfceF = r2;
      if params.fdrAdjustment
        fdrTfceRandF=r2;
      end
    end
  end
  if params.bootstrapStatistics
    bootstrapF = r2;
    if params.fdrAdjustment
      fdrBootstrapF=r2;
    end
    if ~strcmp(params.resampleFWEadjustment,'None')
      bootstrapMaxF = r2;
      if params.TFCE  
        bootstrapTfceF = r2;
        if params.fdrAdjustment
          fdrBootstrapTfceF=r2;
        end
        if ~strcmp(params.resampleFWEadjustment,'None')
          bootstrapMaxTfceF = r2;
        end
      end
    end
  end
end

for scanNum = params.scanNum
  numVolumes = viewGet(thisView,'nFrames',scanNum);
  scanDims = viewGet(thisView,'dims',scanNum);

  %compute the dimensions of the subset of voxels to load
  switch(params.analysisVolume)
   case {'Whole volume','Subset box'}
       subsetBox = eval(scanParams{scanNum}.subsetBox);
   case 'Loaded ROI(s)'
      %get the smallest box containing all the voxels from all the boxes 
      [subsetBox, whichRoi, marginVoxels] = getRoisBox(thisView,scanNum,[voxelsMargin voxelsMargin 0]);
      usedVoxelsInBox = marginVoxels | any(whichRoi,4);
      %clear('whichRoi','marginVoxels');
  end
  subsetDims = diff(subsetBox,1,2)'+1;
  %Compute the number of slices to load at once in order to minimize memory usage
  if params.TFCE && (params.randomizationTests || params.bootstrapStatistics)
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

  ehdr = []; rdf = []; mdf = [];s2 = []; autoCorrelationParameters = [];
  ehdrBootstrapCIs = [];  contrastBootstrapCIs = [];


  % calculate which slices we will be working on
  lastSlices = cumsum(slicesToLoad); 
  firstSlices = lastSlices - slicesToLoad + 1; 
  %loop
  for iBatch = 1:length(loadCallsPerBatch)
    switch(params.analysisVolume)
      case {'Whole volume','Subset box'}
        %Load data
        d = loadScan(thisView,scanNum,[],subsetBox(3,1) + [firstSlices(iBatch) lastSlices(iBatch)] -1, precision,subsetBox(1,:),subsetBox(2,:));
      case 'Loaded ROI(s)'
        for iLoad= loadCallsPerBatch{iBatch}
          dummy = loadScan(thisView,scanNum,[],subsetBox(3,1) + [firstSlices(iLoad) lastSlices(iLoad)] -1, precision,subsetBox(1,:),subsetBox(2,:));
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
      d = getStimvol(d,scanParams{scanNum});
      if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
      % do any call for preprocessing
      if ~isempty(scanParams{scanNum}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{scanNum}.preprocess,verbose);
      end

      actualStimvol = d.stimvol;
      %precompute permutation vectors
      if params.randomizationTests
        nRand = params.nRand;
        %find which event types will be randomized
        %first which EVs are involved in all contrasts/restrictions
        contrastEVs = any([contrasts;cell2mat(restrictions)],1);
        %then which event types constitute these EVs
        randEvents = any(params.scanParams{scanNum}.stimToEVmatrix(:,contrastEVs),2);
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
        randomizations = repmat(cell2mat(actualStimvol),nRand,1);
        for iRand = 1:nRand
           randomizations(iRand,stimvolToPermuteIndices) = stimvolToPermute(randperm(nStimsToPermute));
        end
      else
        nRand = 0;
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
    %we only need to compute Ttests
    computeTtests = params.computeTtests & (params.parametricTests | ... %if parametric tests are asked for
       (length(params.componentsToTest)>1 & strcmp(params.componentsCombination,'Or')) );  % or if contrasts are tested on several components using 'Or' combination
    for iRand = 1:nRand+1
      if iRand==2
        hWaitBar = mrWaitBar(0,['(glmAnalysis) Computing Permutations for scan ' num2str(scanNum)]);
      end
      if iRand == 1 %if it is the actual data
        d.stimvol = actualStimvol;
        computeEstimates = 1;
        computeBootstrap = params.bootstrapStatistics;
        computeTtests = computeTtests || computeBootstrap;
        verbose = 1;
      else
        computeEstimates = 0;
        computeBootstrap = 0;
        mrWaitBar(iRand/nRand,hWaitBar);
        randStimvol = randomizations(iRand-1,:);
        for stimnum = 1:length(actualStimvol)%randomize stim events
           d.stimvol{stimnum} = randStimvol(1:numberEvents(stimnum));
           randStimvol = randStimvol(numberEvents(stimnum)+1:end);
        end
        verbose = 0;
      end

      %compute the design matrix for this randomization
      d = makeDesignMatrix(d,params,verbose, scanNum);
      % compute estimates and statistics
      [d, tempR2, tempC, tempT, tempF, bootstrap] = getGlmStatistics(d, params, verbose, precision, computeEstimates, computeTtests, computeBootstrap);
       
      
      if params.randomizationTests
        thisF = tempF;
        %we'll use the T values for randomizations:
        if ~isempty(contrasts) && params.computeTtests && ...
          (params.parametricTests ... %if parametric tests are asked for (in order to compute resampled FWE adjustments)
        || (length(params.componentsToTest)>1 && strcmp(params.componentsCombination,'Or'))) 
        %or if contrasts are tested on several components using 'Or' combination, since we can't compute contrasts in this case
          thisT = tempT; 
        else %otherwise, we can use the contrast values instead
          thisT = tempC;
        end
      end
      
      if iRand==1
        %keep values for randomization tests
        if params.randomizationTests
          if ~isempty(contrasts) && params.computeTtests
            tempActualT = thisT;
            tempRandT = zeros(size(tempActualT),precision); %will count how many randomization values are above (below) the actual contrast value 
            tempRandT(isnan(thisT)) = NaN; 
            tempRandMaxT = zeros(size(tempActualT),precision); %same for the max (min) randomization value across space, for FWE adjustment
            %this is for resampled-based FWE control
            [indexSortActualT,indexReorderActualT,actualTisNotNaN] = initResampleFWE(reshape(tempActualT,prod(d.dim(1:3)),numberContrasts),params);
          end
          if ~isempty(restrictions)
            tempActualF = thisF;
            tempRandF = zeros(size(tempActualF),precision); %will count how many randomization values are above the actual F value 
            tempRandF(isnan(thisF)) = NaN; 
            tempRandMaxF = zeros(size(tempActualF),precision); %same for the max randomization value across space, for FWE adjustment
            %this is for resampled-based FWE control
            [indexSortActualF,indexReorderActualF,actualFisNotNaN] = initResampleFWE(reshape(tempActualF,prod(d.dim(1:3)),numberFtests),params);
          end
        end
        
        %reshape data into the subsetbox
        if strcmp(params.analysisVolume,'Loaded ROI(s)') 
          d.ehdr = reshapeToRoiBox(d.ehdr,d.roiPositionInBox);
          tempR2 = reshapeToRoiBox(tempR2,d.roiPositionInBox);
          d.s2 = reshapeToRoiBox(d.s2,d.roiPositionInBox);
          if params.covCorrection
            d.autoCorrelationParameters = reshapeToRoiBox(d.autoCorrelationParameters,d.roiPositionInBox);
          end
          if params.bootstrapStatistics && params.bootstrapIntervals
            d.ehdrBootstrapCIs = reshapeToRoiBox(d.ehdrBootstrapCIs,d.roiPositionInBox);
          end
          if ~isempty(contrasts)
            if length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
              tempC = reshapeToRoiBox(tempC,d.roiPositionInBox);
            end
            if params.computeTtests
              if params.parametricTests
                tempT = reshapeToRoiBox(tempT,d.roiPositionInBox);
              end
              if params.bootstrapStatistics
                bootstrap.Tp = reshapeToRoiBox(bootstrap.Tp,d.roiPositionInBox);
                bootstrap.maxTp = reshapeToRoiBox(bootstrap.maxTp,d.roiPositionInBox);
                if params.TFCE 
                  bootstrap.tfceTp = reshapeToRoiBox(bootstrap.tfceTp,d.roiPositionInBox);
                  bootstrap.maxTfceTp = reshapeToRoiBox(bootstrap.maxTfceTp,d.roiPositionInBox);
                end
              end
            end
            if params.bootstrapStatistics && params.bootstrapIntervals && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              d.contrastBootstrapCIs = reshapeToRoiBox(d.contrastBootstrapCIs,d.roiPositionInBox);
            end
          end
          if ~isempty(restrictions) 
            if params.parametricTests
              tempF = reshapeToRoiBox(tempF,d.roiPositionInBox);
              if numel(d.rdf)>1
                d.rdf = reshapeToRoiBox(d.rdf,d.roiPositionInBox);
                d.mdf = reshapeToRoiBox(d.mdf,d.roiPositionInBox);
              end
            end
            if params.bootstrapStatistics
              bootstrap.Fp = reshapeToRoiBox(bootstrap.Fp,d.roiPositionInBox);
              bootstrap.maxFp = reshapeToRoiBox(bootstrap.maxFp,d.roiPositionInBox);
              if params.TFCE 
                bootstrap.tfceFp = reshapeToRoiBox(bootstrap.tfceFp,d.roiPositionInBox);
                bootstrap.maxTfceFp = reshapeToRoiBox(bootstrap.maxTfceFp,d.roiPositionInBox);
              end
            end
          end
        end
        
        %concatenate data
        ehdr = cat(3,ehdr,d.ehdr);
        r2{scanNum} = cat(3,r2{scanNum},tempR2);
        s2 = cat(3,s2,d.s2);
        if params.covCorrection
          autoCorrelationParameters = cat(3,autoCorrelationParameters,d.autoCorrelationParameters);
        end
        if params.bootstrapStatistics && params.bootstrapIntervals
          ehdrBootstrapCIs = cat(3,ehdrBootstrapCIs,d.ehdrBootstrapCIs);
        end
        if ~isempty(contrasts)
          if length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
            C{scanNum} = cat(3,C{scanNum},tempC);
          end
          if params.computeTtests
            if params.parametricTests
              T{scanNum} = cat(3,T{scanNum},tempT);
            end
            if params.bootstrapStatistics
              bootstrapT{scanNum} = cat(3,bootstrapT{scanNum},bootstrap.Tp);
              if ~strcmp(params.resampleFWEadjustment,'None')
                bootstrapMaxT{scanNum} = cat(3,bootstrapMaxT{scanNum},bootstrap.maxTp);
                if params.TFCE 
                  bootstrapTfceT{scanNum} = cat(3,bootstrapTfceT{scanNum},bootstrap.tfceTp);
                  bootstrapMaxTfceT{scanNum} = cat(3,bootstrapMaxTfceT{scanNum},bootstrap.maxTfceTp);
                end
              end
            end
          end
          if params.bootstrapStatistics && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
            contrastBootstrapCIs = cat(3,contrastBootstrapCIs,d.contrastBootstrapCIs);
          end
        end
        if ~isempty(restrictions)
          if params.parametricTests
            F{scanNum} = cat(3,F{scanNum},tempF);
            if numel(d.rdf)>1
               rdf = cat(3,rdf,d.rdf);
               mdf = cat(3,mdf,d.mdf);
            end
          end
          if params.bootstrapStatistics
            bootstrapF{scanNum} = cat(3,bootstrapF{scanNum},bootstrap.Fp);
            if ~strcmp(params.resampleFWEadjustment,'None')
              bootstrapMaxF{scanNum} = cat(3,bootstrapMaxF{scanNum},bootstrap.maxFp);
              if params.TFCE 
                bootstrapTfceF{scanNum} = cat(3,bootstrapTfceF{scanNum},bootstrap.tfceFp);
                bootstrapMaxTfceF{scanNum} = cat(3,bootstrapMaxTfceF{scanNum},bootstrap.maxTfceFp);
              end
            end
          end
        end
        
        clear('tempR2','tempC', 'tempT', 'tempF','bootstrap')
      
      % iRand>1  
      else
        if ~isempty(contrasts)
          tempRandT = tempRandT+double(thisT>tempActualT);
          if ~strcmp(params.resampleFWEadjustment,'None')
            %tempRandMaxT = tempRandMaxT+double(repmat(max(max(max(thisT,[],3),[],2),[],1),[d.dim(1:3) 1])>tempActualT);
            tempRandMaxT = tempRandMaxT + reshape(resampleFWE(reshape(thisT,prod(d.dim(1:3)),numberContrasts),...
                                                              reshape(tempActualT,prod(d.dim(1:3)),numberContrasts),...
                                                              indexSortActualT,indexReorderActualT,actualTisNotNaN,params),...
                                                  [d.dim(1:3) numberContrasts]);
          end
        end
        if ~isempty(restrictions)
          tempRandF = tempRandF+double(thisF>tempActualF);
          %tempRandMaxF = tempRandMaxF+double(repmat(max(max(max(thisF,[],3),[],2),[],1),[d.dim(1:3) 1])>tempActualF);
          if ~strcmp(params.resampleFWEadjustment,'None')
            tempRandMaxF = tempRandMaxF + reshape(resampleFWE(reshape(thisF,prod(d.dim(1:3)),numberFtests),...
                                                              reshape(tempActualF,prod(d.dim(1:3)),numberFtests),...
                                                              indexSortActualF,indexReorderActualF,actualFisNotNaN,params),...
                                                  [d.dim(1:3) numberFtests]);
          end
        end
      end

      %We compute the TFCE here only in the case we forced loading the whole volume at once
      if params.TFCE && params.randomizationTests
        
        if computeTtests && ~isempty(contrasts)
          if isfield(d,'roiPositionInBox') 
            thisT = reshapeToRoiBox(thisT,d.roiPositionInBox);
          end
          tempTfce = applyFslTFCE(thisT,'',0);
          %put NaNs back
          tempTfce(isnan(thisT)) = NaN;
          if iRand == 1 %if it is the actual data
            countTfceT = zeros(size(thisT),precision);     %
            %put NaNs back
            countTfceT(isnan(thisT)) = NaN;
            actualTfceT = tempTfce;
            if ~strcmp(params.resampleFWEadjustment,'None')
              countMaxTfceT = countTfceT;     %
              [indexSortActualTfceT,indexReorderActualTfceT,tfceActualTisNotNaN] = ...
                  initResampleFWE(reshape(actualTfceT,numel(tempTfce(:,:,:,1)),numberContrasts),params);
            end
          else
            countTfceT = countTfceT + double(tempTfce>actualTfceT);
            if ~strcmp(params.resampleFWEadjustment,'None')
              %countMaxTfceT = countMaxTfceT+ double(repmat(max(max(max(tempTfce,[],3),[],2),[],1),[subsetDims 1])>actualTfceT);
              countMaxTfceT = countMaxTfceT + reshape(resampleFWE(reshape(tempTfce,numel(tempTfce(:,:,:,1)),numberContrasts),...
                                                                  reshape(actualTfceT,numel(tempTfce(:,:,:,1)),numberContrasts),...
                                                                  indexSortActualTfceT,indexReorderActualTfceT,tfceActualTisNotNaN,params),...
                                                      size(countMaxTfceT));
            end
          end
        end
        
        %same for F-tests
        if ~isempty(restrictions)
          if isfield(d,'roiPositionInBox') 
            thisF = reshapeToRoiBox(thisF,d.roiPositionInBox); 
          end
          %NaNs in the data will be transformed to 0 by FSL, so for ROIs, 
          %TFCE is applied to the smallest box including all the ROIs, but replacing non-computed data by 0
          tempTfce = applyFslTFCE(thisF,'',0);  
          %put NaNs back
          tempTfce(isnan(thisF)) = NaN;
          if iRand == 1 %if it is the actual data
            countTfceF = zeros(size(thisF),precision);     %these will count how many randomization values are above the actual TFCE values voxelwise
            %put NaNs back
            countTfceF(isnan(thisF)) = NaN;
            actualTfceF = tempTfce;           %keep the TFCE transform
            if ~strcmp(params.resampleFWEadjustment,'None')
              countMaxTfceF = countTfceF; %same for randomization adjusted FWE
              [indexSortActualTfceF,indexReorderActualTfceF,tfceActualFisNotNaN] = ...
                  initResampleFWE(reshape(actualTfceF,numel(tempTfce(:,:,:,1)),numberFtests),params);
            end
          else
            countTfceF = countTfceF+ double(tempTfce>actualTfceF);
            if ~strcmp(params.resampleFWEadjustment,'None')
              %countMaxTfceF = countMaxTfceF+ double(repmat(max(max(max(tempTfce,[],3),[],2),[],1),[subsetDims 1])>actualTfceF);
              countMaxTfceF = countMaxTfceF + reshape(resampleFWE(reshape(tempTfce,numel(tempTfce(:,:,:,1)),numberFtests),...
                                                                  reshape(actualTfceF,numel(tempTfce(:,:,:,1)),numberFtests),...
                                                                  indexSortActualTfceF,indexReorderActualTfceF,tfceActualFisNotNaN,params),...
                                                      size(countMaxTfceF));
            end
          end
        end
        
        clear('tempTfce');
      end
    end
    if params.randomizationTests 
       mrCloseDlg(hWaitBar);
    end
    d = rmfield(d,'data');
    clear('actualTfceT','actualTfceF')

    if params.randomizationTests %compute P-values, reshape and concatenate if needed
      if ~isempty(contrasts)
        tempRandT = tempRandT/nRand;
        if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
            tempRandT = reshapeToRoiBox(tempRandT,d.roiPositionInBox);
        end
        randT{scanNum} = cat(3,randT{scanNum},tempRandT);
        if ~strcmp(params.resampleFWEadjustment,'None')
          tempRandMaxT = tempRandMaxT/nRand;
          tempRandMaxT(isnan(tempRandMaxT)) = NaN; %put NaNs back in place
          tempRandMaxT = reshape(enforceMonotonicityResampleFWE(...
                                    reshape(tempRandMaxT,prod(d.dim(1:3)),numberContrasts),...
                                    indexSortActualT,indexReorderActualT,actualTisNotNaN,params),...
                                    size(tempRandMaxT));
          if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
             tempRandMaxT = reshapeToRoiBox(tempRandMaxT,d.roiPositionInBox);
          end
          randMaxT{scanNum} = cat(3,randMaxT{scanNum},tempRandMaxT);
        end
        if params.TFCE 
          randTfceT{scanNum} = countTfceT/nRand;
          clear('countTfceT');
          randTfceT{scanNum}(isnan(randT{scanNum})) = NaN; %put NaNs back in place
          if ~strcmp(params.resampleFWEadjustment,'None')
            randMaxTfceT{scanNum} = countMaxTfceT/nRand;
            clear('countMaxTfceT');
            randMaxTfceT{scanNum} = reshape(enforceMonotonicityResampleFWE(...
                                        reshape(randMaxTfceT{scanNum},numel(randMaxTfceT{scanNum}(:,:,:,1)),numberContrasts),...
                                        indexSortActualTfceT,indexReorderActualTfceT,tfceActualTisNotNaN,params),...
                                        size(randMaxTfceT{scanNum}));
            randMaxTfceT{scanNum}(isnan(randT{scanNum})) = NaN; %put NaNs back in place
          end
        end
      end
      if ~isempty(restrictions) 
        tempRandF = tempRandF/nRand;
        if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
            tempRandF = reshapeToRoiBox(tempRandF,d.roiPositionInBox);
        end
        randF{scanNum} = cat(3,randF{scanNum},tempRandF);          
        if ~strcmp(params.resampleFWEadjustment,'None')
          tempRandMaxF = tempRandMaxF/nRand;
          tempRandMaxF(isnan(tempRandMaxF)) = NaN; %put NaNs back in place
          if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
            tempRandMaxF = reshapeToRoiBox(tempRandMaxF,d.roiPositionInBox);
          end
          randMaxF{scanNum} = cat(3,randMaxF{scanNum},tempRandMaxF);
        end
        if params.TFCE 
          randTfceF{scanNum} = countTfceF/nRand;
          clear('countTfceF');
          randTfceF{scanNum}(isnan(randF{scanNum})) = NaN; %put NaNs back in place
          if ~strcmp(params.resampleFWEadjustment,'None')
            randMaxTfceF{scanNum} = countMaxTfceF/nRand;
            clear('countMaxTfceF');
            randMaxTfceF{scanNum} = reshape(enforceMonotonicityResampleFWE(...
                                        reshape(randMaxTfceF{scanNum},numel(randMaxTfceF{scanNum}(:,:,:,1)),numberFtests),...
                                        indexSortActualTfceF,indexReorderActualTfceF,tfceActualFisNotNaN,params),...
                                        size(randMaxTfceF{scanNum}));
            randMaxTfceF{scanNum}(isnan(randF{scanNum})) = NaN; %put NaNs back in place
          end
        end
      end
      d = rmfield(d,'roiPositionInBox');
    end
    

  end
  clear('tempRandT','tempRandF','tempActualT','tempActualF')
  
  if params.computeTtests && ~isempty(contrasts)
    if params.parametricTests
        %convert to P value
        Tp{scanNum} = 1 - cdf('t', double(T{scanNum}), d.rdf); %here use doubles to deal with small Ps
        if ~params.outputParametricStatistic
          T{scanNum} = [];
        end
        if strcmp(params.tTestSide,'Both')
          Tp{scanNum} = 2*Tp{scanNum};
        end
        [Tp{scanNum},fdrTp{scanNum}] = transformStatistic(Tp{scanNum},precision,1e-16,params); 
    end
    if params.randomizationTests
      [randT{scanNum},fdrRandC{scanNum}] = transformStatistic(randT{scanNum},precision,1/(nRand+1),params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        randMaxT{scanNum} = transformStatistic(randMaxT{scanNum},precision,1/(nRand+1),params); 
      end
      if params.TFCE 
        [randTfceT{scanNum},fdrTfceRandT{scanNum}] = transformStatistic(randTfceT{scanNum},precision,1/(nRand+1),params); 
        if ~strcmp(params.resampleFWEadjustment,'None')
          randMaxTfceT{scanNum} = transformStatistic(randMaxTfceT{scanNum},precision,1/(nRand+1),params); 
        end
      end
    end
    if params.bootstrapStatistics
      [bootstrapT{scanNum},fdrBootstrapT{scanNum}] = transformStatistic(bootstrapT{scanNum},precision,1/(params.nBootstrap+1),params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        bootstrapMaxT{scanNum} = transformStatistic(bootstrapMaxT{scanNum},precision,1/(params.nBootstrap+1),params);
        if params.TFCE 
          [bootstrapTfceT{scanNum},fdrBootstrapTfceT{scanNum}] = transformStatistic(bootstrapTfceT{scanNum},precision,1/(params.nBootstrap+1),params);
          if ~strcmp(params.resampleFWEadjustment,'None')
            bootstrapMaxTfceT{scanNum} = transformStatistic(bootstrapMaxTfceT{scanNum},precision,1/(params.nBootstrap+1),params);
          end
        end
      end
    end
  end
  
  if ~isempty(restrictions)
    %make parametric probability maps
    if params.parametricTests 
      if ~strcmp(params.testOutput,'T/F value') || params.fdrAdjustment
        %convert to P value
        if numel(d.rdf)>1
          d.rdf = rdf;
          d.mdf = mdf;
        else
          rdf = d.rdf*ones(size(s2));
          mdf = repmat(permute(d.mdf,[1 3 4 2]),[size(s2) 1]);
        end   
        Fp{scanNum} = 1 - cdf('f', double(F{scanNum}), mdf, repmat(rdf,[1 1 1 numberFtests]));  
        clear('mdf','rdf')
        if ~params.outputParametricStatistic
          F{scanNum} = [];
        end
        [Fp{scanNum},fdrFp{scanNum}] = transformStatistic(Fp{scanNum},precision,1e-16,params); 
      end
    end
    if params.randomizationTests 
      [randF{scanNum},fdrRandF{scanNum}] = transformStatistic(randF{scanNum},precision,1/(nRand+1),params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        randMaxF{scanNum} = transformStatistic(randMaxF{scanNum},precision,1/(nRand+1),params); 
      end
      if params.TFCE 
        [randTfceF{scanNum},fdrTfceRandF{scanNum}] = transformStatistic(randTfceF{scanNum},precision,1/(nRand+1),params); 
        if ~strcmp(params.resampleFWEadjustment,'None')
          randMaxTfceF{scanNum} = transformStatistic(randMaxTfceF{scanNum},precision,1/(nRand+1),params); 
        end
      end
    end
    if params.bootstrapStatistics
      [bootstrapF{scanNum},fdrBootstrapF{scanNum}] = transformStatistic(bootstrapF{scanNum},precision,1/(params.nBootstrap+1),params); 
      if ~strcmp(params.resampleFWEadjustment,'None')
        bootstrapMaxF{scanNum} = transformStatistic(bootstrapMaxF{scanNum},precision,1/(params.nBootstrap+1),params);
      end
      if params.TFCE 
        [bootstrapTfceF{scanNum},fdrBootstrapTfceF{scanNum}] = transformStatistic(bootstrapTfceF{scanNum},precision,1/(params.nBootstrap+1),params);
        if ~strcmp(params.resampleFWEadjustment,'None')
          bootstrapMaxTfceF{scanNum} = transformStatistic(bootstrapMaxTfceF{scanNum},precision,1/(params.nBootstrap+1),params);
        end
      end
    end
  end

  % save eventRelated parameters
  glmAnal.d{scanNum} = copyFields(d);
  stimvol = d.stimvol;
  for i=1:length(stimvol)
    stimvol{i} = unique(ceil(stimvol{i}/d.designSupersampling));
  end
  glmAnal.d{scanNum}.stimvol = stimvol;
  glmAnal.d{scanNum}.hrf = downsample(d.hrf, d.designSupersampling/scanParams{scanNum}.estimationSupersampling)/d.designSupersampling*scanParams{scanNum}.estimationSupersampling;
  glmAnal.d{scanNum}.actualhrf = d.hrf;
  clear('d');
  glmAnal.d{scanNum}.estimationSupersampling = scanParams{scanNum}.estimationSupersampling;
  glmAnal.d{scanNum}.acquisitionSubsample = scanParams{scanNum}.acquisitionSubsample;
  glmAnal.d{scanNum}.EVnames = params.EVnames;                %this should be removed if viewGet can get params from GLM analysis
  glmAnal.d{scanNum}.dim = [scanDims numVolumes];
  glmAnal.d{scanNum}.ehdr = NaN([scanDims size(ehdr,4) size(ehdr,5)],precision);
  glmAnal.d{scanNum}.ehdr(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdr;
  glmAnal.d{scanNum}.s2 = NaN(scanDims,precision);
  glmAnal.d{scanNum}.s2(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = s2; %JB
  clear('ehdr','s2');
  if params.covCorrection
    glmAnal.d{scanNum}.autoCorrelationParameters = NaN([scanDims size(autoCorrelationParameters,4)],precision);
    glmAnal.d{scanNum}.autoCorrelationParameters(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = autoCorrelationParameters; %JB
    clear('autoCorrelationParameters')
  end
  if params.bootstrapStatistics && params.bootstrapIntervals
    glmAnal.d{scanNum}.ehdrBootstrapCIs = NaN([scanDims size(ehdrBootstrapCIs,4) size(ehdrBootstrapCIs,5)],precision);
    glmAnal.d{scanNum}.ehdrBootstrapCIs(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdrBootstrapCIs; %JB
    clear('ehdrBootstrapCIs')
  end
  if ~isempty(contrasts)
    glmAnal.d{scanNum}.contrasts = contrasts;
    if params.bootstrapStatistics && params.bootstrapIntervals  && (length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
      glmAnal.d{scanNum}.contrastBootstrapCIs = NaN([scanDims size(contrastBootstrapCIs,4) size(contrastBootstrapCIs,5)],precision);
      glmAnal.d{scanNum}.contrastBootstrapCIs(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = contrastBootstrapCIs; %JB
      clear('contrastBootstrapCIs')
    end
  end
  if ~isempty(restrictions)
    glmAnal.d{scanNum}.fTestNames = params.fTestNames;        %this should be removed if viewGet can get params from GLM analysis
    glmAnal.d{scanNum}.restrictions = restrictions;               %this should be removed if viewGet can get params from GLM analysis
  end


end

tic
%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
glmAnal.name = params.saveName;
if strcmp(params.hrfModel,'hrfDeconvolution')
  glmAnal.type = 'deconvAnal';
elseif isempty(restrictions) && ~params.computeTtests
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
for scanNum = params.scanNum
   defaultOverlay.data{scanNum} = NaN(scanDims,precision); %to make values outside the box transparent
end


%------------------------------------------------------ save the r2 overlay
overlays = defaultOverlay;
overlays.name = 'r2';
overlays.colormapType = 'setRangeToMax';
for scanNum = params.scanNum
   overlays.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) = r2{scanNum};
   overlays.params{scanNum} = scanParams{scanNum};
end
clear('r2');

%--------------------------------------------- save the contrast overlay(s)

if ~isempty(contrasts)
  
  contrastNames = makeContrastNames(contrasts,params.EVnames,params.tTestSide);
  
  if params.computeTtests
    if params.parametricTests
      
      if params.outputParametricStatistic
        thisOverlay = defaultOverlay;
        thisOverlay.colormap = jet(256);
        max_abs_T = max(max(max(max(max(abs(cell2mat(T)))))));
        thisOverlay.range = [-max_abs_T max_abs_T];
        namePrefix = 'T ( ';
        betaAlphaOverlayExponent = 0;

        thisOverlay.clip = thisOverlay.range;
        for iContrast = 1:numberContrasts
          overlays(end+1)=thisOverlay;
          overlays(end).name = [namePrefix contrastNames{iContrast} ')'];
          for scanNum = params.scanNum
            overlays(end).data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
               T{scanNum}(:,:,:,iContrast);
            overlays(end).params{scanNum} = scanParams{scanNum};
          end
        end
      end
      clear('T');
      
      overlays = [overlays makeOverlay(defaultOverlay, Tp, subsetBox, params.scanNum, scanParams, ...
                                         '', params.testOutput, 'T', contrastNames)];
      clear('Tp');

      if ~strcmp(params.resampleFWEadjustment,'None')
        if params.bootstrapStatistics
          overlays = [overlays makeOverlay(defaultOverlay, bootstrapMaxT, subsetBox, params.scanNum, scanParams, ...
                                             'bootstrap-adjusted ', params.testOutput, 'T', contrastNames)];
          clear('bootstrapMaxT');
          if params.TFCE  
            if ~strcmp(params.resampleFWEadjustment,'None')
              overlays = [overlays makeOverlay(defaultOverlay, bootstrapMaxTfceT, subsetBox, params.scanNum, scanParams, ...
                                                 'bootstrap-adjusted ', params.testOutput, 'TFCE T', contrastNames)];
              clear('bootstrapMaxTfceT');
            end
          end
        end 
        if params.randomizationTests
          overlays = [overlays makeOverlay(defaultOverlay, randMaxT, subsetBox, params.scanNum, scanParams, ...
                                             'randomization-adjusted ', params.testOutput, 'T', contrastNames)];
          clear('randMaxT');
          if params.TFCE
            overlays = [overlays makeOverlay(defaultOverlay, randMaxTfceT, subsetBox, params.scanNum, scanParams, ...
                                             'randomization-adjusted ',params.testOutput, 'TFCE T', contrastNames)];
            clear('randMaxTfceT');
          end
        end
      end
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTp, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted ',params.testOutput, 'T', contrastNames)];
        clear('fdrTp');
      end
      
    end
    
    if params.bootstrapStatistics
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapT, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap ', params.testOutput, 'T', contrastNames)];
      clear('bootstrapT');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapT, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted bootstrap ',params.testOutput, 'T', contrastNames)];
        clear('fdrBootstrapTp');
      end
      if params.TFCE  
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapTfceT, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap ', params.testOutput, 'TFCE T', contrastNames)];
        clear('bootstrapTfceT');
        if params.fdrAdjustment
          overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapTfceT, subsetBox, params.scanNum, scanParams, ...
                                             'FDR-adjusted bootstrap ', params.testOutput, 'TFCE T', contrastNames)];
          clear('fdrBootstrapTfceT');
        end
      end
    end
    
    if params.randomizationTests
      overlays = [overlays makeOverlay(defaultOverlay, randT, subsetBox, params.scanNum, scanParams, ...
                                         'randomization ', params.testOutput, 'T', contrastNames)];
      clear('randT');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrRandC, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted randomization ',params.testOutput, 'T', contrastNames)];
        clear('fdrRandC');
      end
      if params.TFCE
        overlays = [overlays makeOverlay(defaultOverlay, randTfceT, subsetBox, params.scanNum, scanParams, ...
                                        'randomization ',params.testOutput, 'TFCE T', contrastNames)];
        clear('randTfceT');
        if params.fdrAdjustment
          overlays = [overlays makeOverlay(defaultOverlay, fdrTfceRandT, subsetBox, params.scanNum, scanParams, ...
                                          'FDR-adjusted randomization ',params.testOutput, 'TFCE T', contrastNames)];
          clear('fdrTfceRandT');
        end
      end
    end
    
    %this is to mask the beta values by the probability/Z maps     
    betaAlphaOverlay = cell(length(contrastNames),1);
    if params.parametricTests || params.bootstrapStatistics || params.randomizationsTests
      %set the contrast alpha overlay to the statistical value
      switch(params.testOutput)
        case 'P value'                                                  %statistical maps
          betaAlphaOverlayExponent = -1;      % with inverse masking for p values
        case {'Z value','-log10(P) value}'}
          betaAlphaOverlayExponent = .5;      %or normal masking for Z or log10(p) values
      end
      for iContrast = 1:numberContrasts
        betaAlphaOverlay{iContrast} = overlays(end-numberContrasts+iContrast).name;
      end
    else
      betaAlphaOverlayExponent = 1;
    end
  end
  
%--------------------------------------------- save the contrast beta weights overlay(s) (ehdr if no contrast)

  if length(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
    %find the values for the scale of the beta overlays
    thisOverlay = defaultOverlay;
    ordered_abs_betas = cell2mat(C(params.scanNum));
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
      for scanNum = params.scanNum
        overlays(end).data{scanNum} = NaN(scanDims,precision); %to make values outside the box transparent
        overlays(end).data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
           C{scanNum}(:,:,:,iContrast);
        overlays(end).params{scanNum} = scanParams{scanNum};
      end
      overlays(end).alphaOverlayExponent=betaAlphaOverlayExponent;
    end
    clear('C');
  end
end

%----------------------------------------------- save the F-test overlay(s)
if ~isempty(restrictions)
  if params.parametricTests
    if params.outputParametricStatistic
      thisOverlay = defaultOverlay;
      thisOverlay.range(1) = min(min(min(min(min(cell2mat(F))))));
      thisOverlay.range(2) = min(max(max(max(max(cell2mat(F))))));
      namePrefix = 'F ( ';
      thisOverlay.clip = thisOverlay.range;
      for iFtest = 1:numberFtests
        %probability maps
        overlays(end+1)=thisOverlay;
        overlays(end).name = [namePrefix params.fTestNames{iFtest} ')'];
        for scanNum = params.scanNum
          overlays(end).data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             F{scanNum}(:,:,:,iFtest);
          overlays(end).params{scanNum} = scanParams{scanNum};
        end
      end
      clear('F');
    end

    overlays = [overlays makeOverlay(defaultOverlay, Fp, subsetBox, params.scanNum, scanParams, ...
                                       '', params.testOutput, 'F', params.fTestNames)];
    clear('Fp');

    if ~strcmp(params.resampleFWEadjustment,'None')
      if params.bootstrapStatistics
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapMaxF, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-adjusted ', params.testOutput, 'F', params.fTestNames)];
        clear('bootstrapMaxF');
         if params.TFCE  
          if ~strcmp(params.resampleFWEadjustment,'None')
            overlays = [overlays makeOverlay(defaultOverlay, bootstrapMaxTfceF, subsetBox, params.scanNum, scanParams, ...
                                               'bootstrap-adjusted ', params.testOutput, 'TFCE F', params.fTestNames)];
            clear('bootstrapMaxTfceF');
          end
         end
      end
      if params.randomizationTests
        overlays = [overlays makeOverlay(defaultOverlay, randMaxF, subsetBox, params.scanNum, scanParams, ...
                                           'randomization-adjusted ', params.testOutput, 'F', params.fTestNames)];
        clear('randMaxF');
        if params.TFCE
          overlays = [overlays makeOverlay(defaultOverlay, randMaxTfceF, subsetBox, params.scanNum, scanParams, ...
                                           'randomization-adjusted ',params.testOutput, 'TFCE F', params.fTestNames)];
          clear('randMaxTfceF');
        end
      end
    end
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrFp, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted ',params.testOutput, 'F', params.fTestNames)];
      clear('fdrFp');
    end

  end
  
  if params.bootstrapStatistics
    overlays = [overlays makeOverlay(defaultOverlay, bootstrapF, subsetBox, params.scanNum, scanParams, ...
                                     'bootstrap ',params.testOutput, 'F', params.fTestNames)];
    clear('bootstrapF');
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapF, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted bootstrap ',params.testOutput, 'F', params.fTestNames)];
      clear('fdrBootstrapF');
    end
    if params.TFCE  
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapTfceF, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap ', params.testOutput, 'TFCE F', params.fTestNames)];
      clear('bootstrapTfceF');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapTfceF, subsetBox, params.scanNum, scanParams, ...
                                           'FDR-adjusted bootstrap ', params.testOutput, 'TFCE F', params.fTestNames)];
        clear('fdrBootstrapTfceF');
      end
    end
  end
  
  if params.randomizationTests
    overlays = [overlays makeOverlay(defaultOverlay, randF, subsetBox, params.scanNum, scanParams, ...
                                     'randomization ',params.testOutput, 'F', params.fTestNames)];
    clear('randF');
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrRandF, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted randomization ',params.testOutput, 'F', params.fTestNames)];
      clear('fdrRandF');
    end
    
    if params.TFCE
      overlays = [overlays makeOverlay(defaultOverlay, randTfceF, subsetBox, params.scanNum, scanParams, ...
                                       'randomization ',params.testOutput, 'TFCE F', params.fTestNames)];
      clear('randTfceF');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTfceRandF, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted randomization ',params.testOutput, 'TFCE F', params.fTestNames)];
        clear('fdrTfceRandF');
      end
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



function [convertedStatistic, fdrAdjustedStatistic] = transformStatistic(p, outputPrecision, minP, params)

convertedStatistic = convertStatistic(p, params.testOutput, outputPrecision, minP);

if params.fdrAdjustment
  fdrAdjustedP = p;
  for iTest = 1:size(p,4)
    fdrAdjustedP(:,:,:,iTest) = fdrAdjust(p(:,:,:,iTest),params);
  end
  fdrAdjustedStatistic = convertStatistic(fdrAdjustedP, params.testOutput, outputPrecision, minP);
else
  fdrAdjustedStatistic = [];
end


function convertedP = convertStatistic(p, outputStatistic, outputPrecision, minP)
%we do not allow probabilities of 0 and replace them by minP
%(this can occur because cdf cannot return values less than 1e-16, 
%or because no bootstrap resampling return a value as high as the actual value)
convertedP = max(p,minP);
convertedP(isnan(p)) = NaN; %NaNs must remain NaNs (they became minP when using max)

switch(outputStatistic)
  case 'Z value'
    %convert P value to Z value, 
    convertedP = norminv(1-convertedP);  
    convertedP(convertedP<0) = 0; % we're not interested in negative Z value
    %if there was no round-off error from cdf, we could do the following:
    %Z = max(-norminv(p),0);  %because the left side of norminv seems to be less sensitive to round-off errors,
    %we get -norminv(x) instead of norminv(1-x). also we'renot interested in negative Z value
  case '-log10(P) value'
    %convert P to -log10(P) 
    convertedP = -log10(convertedP);
    
end
convertedP = eval([outputPrecision '(convertedP)']);



function adjustedPdata = fdrAdjust(pData,params)

isNotNan = ~isnan(pData);
sizePdata = size(pData);
pData = pData(isNotNan);
[pData,sortingIndex] = sort(pData);

switch(params.fdrAssumption)
  case 'Independence/Positive dependence'
    constant =1;
  case 'None'
    constant = sum(1./(1:nTests));
    %approximate value: 
    %constant = log(nTests)+0.57721566490153 %Euler constant
end

%estimate number of True Null hypotheses
% ref: Benjamini et al. 2006, Biometrika, 93(3) p491
switch(params.fdrMethod)
  case 'Linear Step-up'
    trueHoRatio = 1; %no estimation of the number of true Null hypotheses over the total number of rejected test
  case 'Two-stage Linear Step-up'
  case 'Multiple-stage Linear Step-up'
  case 'Adaptive Linear Step-up'
    %we will estimate the number of true H0 using definition 3 in Benjamini et al. 2006
    % first adjust the p values 
    adjustedP = linearStepUpFdr(pData,constant,1);
    if ~any(adjustedP>params.fdrMultipleStepQ) %if no voxel is significant at a chosen level (usually .05) (very unlikely)
      trueHoRatio = 0; % No H0 will be rejected
    else
      %compute the estimated number of true H0 when considering the number of rejectionsa t all possible levels
      numberH0 = length(adjustedP);
      numberTrueHo = (numberH0+1-(1:numberH0)')./(1-pData);
      %find the first estimated number that is greater than the previous one
      numberTrueHo = ceil(min(numberTrueHo(find(diff(numberTrueHo)>0,1,'first')+1),numberH0));
      trueHoRatio = numberH0/numberTrueHo; 
    end
      
end

if trueHoRatio ==0
  adjustedPdata = ones(sizePdata); 
else
  adjustedP = linearStepUpFdr(pData,constant,trueHoRatio);
  adjustedP(sortingIndex) = adjustedP;
  adjustedPdata = NaN(sizePdata);
  adjustedPdata(isNotNan) = adjustedP;
end



function qData = linearStepUpFdr(pData,constant,trueHoRatio)

nTests = length(pData);
%adjustment (does not depend on a chosen threshold)
% it consists in finding, for each p-value, 
% the largest FDR threshold such that the p-value would be considered significant
% it is the inverse operation to p-correction (see commented code below)
% ref: Yekutieli and Benjamini 1999, Journal of Statistical Planning and Inference, 82(1-2) p171
qData = min(pData.*nTests./(1:nTests)'.*constant/trueHoRatio,1);
%make it monotonic from the end (for each i, q(i) must be smaller than all q between i and n)
minQ = qData(end);
for i=nTests-1:-1:1
 qData(i) = min(minQ,qData(i));
 minQ = qData(i);
end

% % % %correction (p-correction depends on a fixed threshold, here .05)
% % % qAlpha = .05*(1:nTests)'/nTests/constant;
% % % correctedP = pData;
% % % correctedP(find(pData<=qAlpha,1,'last'):end) = 1;
% % % correctedP(sortingIndex) = correctedP;
% % % correctedPdata = NaN(size(pData));
% % % correctedPdata(isNotNan) = correctedP;

function overlays = makeOverlay(overlays,data,subsetBox,scanList,scanParams,testString1, testOutput, testString2, testNames)
  switch testOutput
    case 'P value'
      overlays.colormap = statsColorMap(256);
      namePrefix = [testString1 'P (' testString2 ' '];
    case 'Z value'
      overlays.range = [0 8.2096];
      namePrefix = [testString1 'Z (' testString2 ' '];
    case '-log10(P) value'
      overlays.range = [0 16];
      namePrefix = ['-log10(' testString1 'P) (' testString2 ' '];
  end
  overlays.clip = overlays.range;
  overlays = repmat(overlays,1,length(testNames));
  for iTest = 1:length(testNames)
    overlays(iTest).name = [namePrefix testNames{iTest} ')'];
    for scanNum = scanList
      overlays(iTest).data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
         data{scanNum}(:,:,:,iTest);
      overlays(iTest).params{scanNum} = scanParams{scanNum};
    end
  end



