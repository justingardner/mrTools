% eventRelatedGlm.m
%
%      usage: thisView = eventRelatedGlm(thisView,params)
%         by: farshad moradi
%       date: 06/14/07
%    purpose: same as eventRelated, but uses canonical hrf instead of
%             deconvolution
%              $Id$
%
%  modified by julien besle on 12/01/2010 to perform an ANOVAs and T-tests
%     Anovas are limited to subsets of Explanatory Variables (EVs) (H0 = no EV differs from zero)
%     T-tests are performed on any linear combinations of EVs against zero

function [thisView d] = eventRelatedGlm(thisView,params,varargin)

d = [];

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help eventRelatedGlm
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
   params = eventRelatedGlmGUI('thisView',thisView,'defaultParams',defaultParams,'scanList',scanList);
end

if justGetParams
  d = params;
  return
end

% Abort if params empty
if ieNotDefined('params'),return,end


% set the group
thisView = viewSet(thisView,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

%just to have shorter variable names
scanParams = params.scanParams;
testParams = params.testParams;
fTests = testParams.fTests;
contrasts = testParams.contrasts;

if params.covCorrection   %number of voxels to get around the ROI/subset box in case the covariance matrix is estimated
   voxelsMargin = floor(params.covEstimationAreaSize/2);
else
   voxelsMargin = 0;
end

if any( sum(contrasts,2)==1 & sum(contrasts~=0,2)==1 )
   scanNum = params.scanNum(1);
   disp(sprintf('Getting condition names from scan %d', scanNum)); 
   d = loadScan(thisView, scanNum, [], 0);
   d = getStimvol(d,scanParams{scanNum});
   % do any call for preprocessing
   if ~isempty(scanParams{scanNum}.preprocess)
      d = eventRelatedPreProcess(d,scanParams{scanNum}.preprocess);
   end
end


%--------------------------------------------------------- Main loop over scans ---------------------------------------------------
tic
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
%initialize the data we're keeping for output overlays
precision = mrGetPref('defaultPrecision');
r2 = cell(1,params.scanNum(end));
contrastData = r2;
tData = r2;
contrastRandData = r2;  
tDataTFCE = r2;
tRandTfce = r2;
thresholdTfceT = NaN(size(contrasts,1),params.scanNum(end));
fRandData = r2;
fData = r2;  
fDataTFCE = r2;
fRandTfce = r2;
thresholdTfceF = NaN(size(fTests,1),params.scanNum(end));

for scanNum = params.scanNum
  numVolumes = viewGet(thisView,'nFrames',scanNum);
  scanDims = viewGet(thisView,'dims',scanNum);

  %compute the dimensions of the subset of voxels to load
  switch(params.analysisVolume)
   case {'Whole volume','Subset box'}
       subsetBox = eval(scanParams{scanNum}.subsetBox);
   case 'Loaded ROI(s)'
      %get the smallest box containing all the voxels from all the boxes 
      [subsetBox, whichRoi, marginVoxels] = getRoisBox(thisView,scanNum,voxelsMargin);
      roiPositionInBox = marginVoxels | any(whichRoi,4);
      clear('whichRoi','marginVoxels');
  end
  subsetDims = diff(subsetBox,1,2)+1;
  %Compute the number of slices to load at once in order to minimize memory usage
  if testParams.TFCE && testParams.randomizationTests
   %in this particular case, we force loading the whole dataset at once
   slicesToLoad = subsetDims(3);
   loadCallsPerBatch = {1};
   rawNumSlices = subsetDims(3);
   %compare the size of the array to load to the memory preference and ask for confirmation if 
   switch (precision)
     case {'double'}
      bytesPerNum = 8;
     case{'single'}
      bytesPerNum = 4;
   end
   sizeArray = bytesPerNum*numVolumes*prod(subsetDims);
   if sizeArray> mrGetPref('maxBlocksize')
      if ~askuser(sprintf('(eventRelatedGlm) This will load an array of %.2f Gb in memory. Are you sure you want to proceed ?', sizeArray/1024^3));
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
          nVoxels =  nVoxels+nnz(roiPositionInBox(:,:,currentFirstSlice+(1:slicesToLoad(iLoad))) > 0);
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
      mrWarnDlg(['Too many data points to perform the analysis on whole slices. Implement row analysis support, increase memory block size (by a factor ' num2str(subsetDims(3)/rawNumSlices) ') or reduce subset box']);
      return;
  end

  rdf = []; mdf = [];rss = []; 
  ehdr = []; ehdrste = [];

  %----------------------------------------Design matrix: same for all voxels
  % get the stim volumes, if empty then abort
  d = loadScan(thisView, scanNum, [], 0); %load scan info
  d = getStimvol(d,scanParams{scanNum});
  if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end

  actualStimvol = d.stimvol;
  %precompute randomizations
  if testParams.randomizationTests
    nRand = testParams.nRand;
    randAlpha = .05;
    nEventTypes = length(actualStimvol);
    nStimPerEventType = NaN(nEventTypes,1);
    for iEventType = 1:nEventTypes
       nStimPerEventType(iEventType) = length(actualStimvol{iEventType});
    end
    all_stimvol = cell2mat(actualStimvol);
    nStims = length(all_stimvol);
    randomizations = NaN(nRand-1,nStims);
    for iRand = 1:nRand-1
       randomizations(iRand,:) = randperm(nStims);
    end
  else
    nRand = 1;
  end

  %create model HRF
  [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.supersampling,0,1);

  d.volumes = 1:d.dim(4);

  %make a copy of d
  actual_d = d;

  % calculate which slices we will be working on
  lastSlices = cumsum(slicesToLoad); 
  firstSlices = lastSlices - slicesToLoad + 1; 
  %loop
  for iBatch = 1:length(loadCallsPerBatch)
    d = actual_d;
    switch(params.analysisVolume)
       case {'Whole volume','Subset box'}
          %Load data
          dummy = loadScan(thisView,scanNum,[],subsetBox(3,1) + [firstSlices(iBatch) lastSlices(iBatch)] -1, precision,subsetBox(1,:),subsetBox(2,:));
          d.data = dummy.data;
          d.dim = dummy.dim;
       case 'Loaded ROI(s)'
          d.data = [];
          for iLoad= loadCallsPerBatch{iBatch}
             dummy = loadScan(thisView,scanNum,[],subsetBox(3,1) + [firstSlices(iLoad) lastSlices(iLoad)] -1, precision,subsetBox(1,:),subsetBox(2,:));
             dummy.data = reshape(dummy.data,[prod(dummy.dim(1:3)) dummy.dim(4)]);
             d.data = cat(1,d.data, dummy.data(roiPositionInBox(:,:,firstSlices(iLoad):lastSlices(iLoad))>0,:) );
          end
          d.data = permute(d.data,[1 3 4 2]);
          d.roiPositionInBox = roiPositionInBox(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)));
          d.dim = size(d.data);
    end
    clear dummy;

    %-------------------------------Permutations------------------------------------------
    for iRand = 1:nRand
      if iRand==2
        hWaitBar = mrWaitbar(0,['Computing Permutations for scan ' num2str(scanNum)]);
      end
      if iRand == 1 %if it is the actual data
        d.stimvol = actualStimvol;
        computeTtests = testParams.parametricTests & params.computeTtests;
        computeEstimates = 1;
        verbose = 1;
      else
        %we only need to compute Ttests
        computeTtests = params.computeTtests & (testParams.TFCE  |... %in the case we compute TFCE
          (length(testParams.componentsToTest)>1 & strcmp(testParams.componentsCombination,'Or') & testParams.parametricTests));  % or if contrasts are tested on several components using 'Or' combination

        computeEstimates = 0;
        mrWaitbar(iRand/nRand,hWaitBar);
        randStimvol = all_stimvol(randomizations(iRand-1,:));
        for stimnum = 1:length(actualStimvol)%randomize stim events
           d.stimvol{stimnum} = randStimvol(1:nStimPerEventType(stimnum));
           randStimvol = randStimvol(nStimPerEventType(stimnum)+1:end);
        end
        verbose = 0;
      end

      % do any call for preprocessing
      if ~isempty(scanParams{scanNum}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{scanNum}.preprocess,verbose);
      end
      %compute the design matrix for this randomization
      d = makeglm(d,params,verbose);
      % compute estimates and statistics
      [d, contrastDataTemp, tDataTemp, fDataTemp] = getGlmStatistics(d, params, verbose, precision, computeEstimates, computeTtests);

      if ~isempty(contrasts)
        %in the particular case where contrasts are tested on several components using 'Or' combination, we don't compute the contrasts
        %so we'll use the T-test (in fact F-test) values instead of the contrast values
        if length(testParams.componentsToTest)>1 && strcmp(testParams.componentsCombination,'Or')
          contrastDataTemp = tDataTemp;
        end
          
        if iRand==1
          actualTdataTemp = tDataTemp;
          actualCdataTemp = contrastDataTemp;
          if testParams.randomizationTests
            contrastRandDataTemp = zeros(size(actualCdataTemp),precision); %will count how many randomization values are above the actual contrast value 
          end
        else
          switch testParams.tTestSide
            case 'Right'
              contrastRandDataTemp = contrastRandDataTemp+double(contrastDataTemp>actualCdataTemp);
            case 'Left'
              contrastRandDataTemp = contrastRandDataTemp+double(contrastDataTemp<actualCdataTemp);
            case 'Both'
              contrastRandDataTemp = contrastRandDataTemp+double(abs(contrastDataTemp)>abs(actualCdataTemp));
          end
        end
      end
      if ~isempty(fTests)
        if iRand==1
          actualFdataTemp = fDataTemp;
          if testParams.randomizationTests
            fRandDataTemp = zeros(size(actualFdataTemp),precision); %will count how many randomization values are above the actual F value 
          end
        else
          fRandDataTemp = fRandDataTemp+double(fDataTemp>actualFdataTemp);
        end
      end

      %We compute the TFCE here only in the case we forced loading the whole volume at once
      if testParams.TFCE && testParams.randomizationTests
        
        if ~isempty(fTests)
          if isfield(d,'roiPositionInBox') 
            fDataTemp = reshapeToRoiBox(fDataTemp,d.roiPositionInBox); 
          end
          %NaNs in the data will be transformed to 0 by FSL, so for ROIs, 
          %TFCE is applied to the smallest box including all the ROIs, but replacing non-computed data by 0
          tfceData = applyFslTFCE(fDataTemp,'',0);   
          tfceData(isnan(fDataTemp)) = NaN; %put NaNs back in place
          if iRand == 1 %if it is the actual data
            fDataTfceCount = NaN(size(fDataTemp),precision);     %these will count how many randomization values are above the actual TFCE values voxelwise
            fDataTfceCount(~isnan(fDataTemp)) = 1; %we count the actual data in
            maxTfceF = zeros(nRand,size(fTests,1));       %and these will store the distribution of max
            fDataTFCE{scanNum} = tfceData;           %keep the TFCE transform
          else
            fDataTfceCount = fDataTfceCount+ double(tfceData>fDataTFCE{scanNum});
          end
          maxTfceF(iRand,:) = permute(max(max(max(fDataTFCE{scanNum}))),[1 4 2 3]);
        end
        
        %same for T-tests
        if computeTtests && ~isempty(contrasts)
          if isfield(d,'roiPositionInBox') 
            tDataTemp = reshapeToRoiBox(tDataTemp,d.roiPositionInBox);
          end
          tfceData = applyFslTFCE(tDataTemp,'',0);
          tfceData(isnan(tDataTemp)) = NaN; %put NaNs back in place
          if iRand == 1 %if it is the actual data
            tDataTfceCount = NaN(size(tDataTemp),precision);     %
            tDataTfceCount(~isnan(tDataTemp)) = 1;
            maxTfceT = zeros(nRand,size(contrasts,1));
            tDataTFCE{scanNum} = tfceData; 
          else
            tDataTfceCount = tDataTfceCount+ double(tfceData>tDataTFCE{scanNum});
          end
          maxTfceT(iRand,:) = permute(max(max(max(tDataTFCE{scanNum}))),[1 4 2 3]);
        end
        
      end
    end
    if testParams.randomizationTests 
       mrCloseDlg(hWaitBar);
    end

    if strcmp(params.analysisVolume,'Loaded ROI(s)') %reshape data into the subsetbox
      d.ehdr = reshapeToRoiBox(d.ehdr,d.roiPositionInBox);
      d.ehdrste = reshapeToRoiBox(d.ehdrste,d.roiPositionInBox);
      d.r2 = reshapeToRoiBox(d.r2,d.roiPositionInBox);
      d.rss = reshapeToRoiBox(d.rss,d.roiPositionInBox);
      if ~isempty(contrasts)
        actualCdataTemp = reshapeToRoiBox(actualCdataTemp,d.roiPositionInBox);
        if params.computeTtests
          actualTdataTemp = reshapeToRoiBox(actualTdataTemp,d.roiPositionInBox);
        end
        if testParams.randomizationTests
          contrastRandDataTemp = reshapeToRoiBox(contrastRandDataTemp,d.roiPositionInBox);
        end
      end
      if ~isempty(fTests) 
        actualFdataTemp = reshapeToRoiBox(actualFdataTemp,d.roiPositionInBox);
        if numel(d.rdf)>1
          d.rdf = reshapeToRoiBox(d.rdf,d.roiPositionInBox);
          d.mdf = reshapeToRoiBox(d.mdf,d.roiPositionInBox);
        end
        if testParams.randomizationTests
          fRandDataTemp = reshapeToRoiBox(fRandDataTemp,d.roiPositionInBox);
        end
      end
      d = rmfield(d,'roiPositionInBox');
    end

    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    ehdrste = cat(3,ehdrste,d.ehdrste);
    r2{scanNum} = cat(3,r2{scanNum},d.r2);
    rss = cat(3,rss,d.rss);
  %      tss = cat(3,tss,d.tss);
    if ~isempty(contrasts)
      if length(testParams.componentsToTest)==1 || strcmp(testParams.componentsCombination,'Add')
        contrastData{scanNum} = cat(3,contrastData{scanNum},actualCdataTemp);
      end
      if params.computeTtests
        tData{scanNum} = cat(3,tData{scanNum},actualTdataTemp);
      end
      if testParams.randomizationTests
        contrastRandData{scanNum} = cat(3,contrastRandData{scanNum},contrastRandDataTemp);
      end
    end
    if ~isempty(fTests)
      fData{scanNum} = cat(3,fData{scanNum},actualFdataTemp);
      if numel(d.rdf)>1
         rdf = cat(3,rdf,d.rdf);
         mdf = cat(3,mdf,d.mdf);
      end
      if testParams.randomizationTests
        fRandData{scanNum} = cat(3,fRandData{scanNum},fRandDataTemp);
      end
    end

  end
  clear('contrastDataTemp', 'tDataTemp', 'fDataTemp','contrastRandDataTemp','fRandDataTemp','actualCdataTemp','actualFdataTemp')
  
  %compute the TFCE if randomization test not run on it
  if testParams.TFCE && ~testParams.randomizationTests
    if ~isempty(fTests)
      fDataTFCE{scanNum} = NaN(size(fData{scanNum}),precision);     
      fDataTFCE{scanNum} = applyFslTFCE(fData{scanNum}); 
      fDataTFCE{scanNum}(isnan(fData{scanNum})) = NaN; %put NaNs back in place
    end
    if computeTtests && ~isempty(contrasts)
      tDataTFCE{scanNum} = NaN(size(tData{scanNum}),precision);    
      tDataTFCE{scanNum} = applyFslTFCE(tData{scanNum}); 
      tDataTFCE{scanNum}(isnan(tData{scanNum})) = NaN; %put NaNs back in place
    end
  end

  if params.computeTtests && ~isempty(contrasts)
    if testParams.parametricTests
      if ismember(testParams.parametricTestOutput,{'P value','Z value'})
         %convert to P value
         tData{scanNum} = 1 - cdf('t', double(tData{scanNum}), d.rdf); %here use doubles to deal with small Ps
         if strcmp(testParams.tTestSide,'Both')
            tData{scanNum} = 2*tData{scanNum};
         end
        if strcmp(testParams.parametricTestOutput,'Z value')
           %probabilities output by cdf will be 1 if (1-p)<1e-16, which will give a Z value of +/-inf 
           tData{scanNum} = pToZ(tData{scanNum},1e-16,precision); %this will be handled by pToZ
         end
      end
      if testParams.randomizationTests
        contrastRandData{scanNum} = contrastRandData{scanNum}/nRand;
        if strcmp(testParams.randomizationTestOutput,'Z value')
          contrastRandData{scanNum} = pToZ(contrastRandData{scanNum},1/nRand,precision); 
        end
        %compute TFCE thresholds
        if testParams.TFCE 
           for iContrast = 1:size(contrasts,1)
              sorted_max_tfce = sort(maxTfceT(:,iContrast));
              thresholdTfceT(iContrast,scanNum) = sorted_max_tfce(max(1,floor((1-randAlpha)*nRand)));
           end
          tRandTfce{scanNum} = tDataTfceCount/nRand;
          if strcmp(testParams.randomizationTestOutput,'Z value')
            tRandTfce{scanNum} = pToZ(tRandTfce{scanNum},1/nRand,precision); 
          end
        end
      end
    end
  end
  
  if ~isempty(fTests)
    %make parametric probability maps
    if testParams.parametricTests 
      if ismember(testParams.parametricTestOutput,{'P value','Z value'})
        if numel(d.rdf)>1
          d.rdf = rdf;
          d.mdf = mdf;
        else
          rdf = d.rdf*ones(size(rss));
          mdf = repmat(permute(d.mdf,[1 3 4 2]),[size(rss) 1]);
        end   
        fData{scanNum} = 1 - cdf('f', double(fData{scanNum}), mdf, repmat(rdf,[1 1 1 size(fTests,1)]));  
        clear('mdf','rdf')
        if strcmp(testParams.parametricTestOutput,'Z value')
          %probabilities output by cdf are 1 if (1-p)<1e-16, which will give a Z value of +/-inf
          fData{scanNum} = pToZ(fData{scanNum},1e-16,precision); %this will be handled by pToZ
        end
      end
      if testParams.randomizationTests 
        fRandData{scanNum} = fRandData{scanNum}/nRand;
        if strcmp(testParams.randomizationTestOutput,'Z value')
          fRandData{scanNum} = pToZ(fRandData{scanNum},1/nRand,precision); 
        end
        %compute TFCE thresholds
        if testParams.TFCE 
          for iFtest = 1:size(fTests,1)
            sorted_max_tfce = sort(maxTfceF(:,iFtest));
            thresholdTfceF(iFtest,scanNum) = sorted_max_tfce(max(1,floor((1-randAlpha)*nRand)));
          end
          fRandTfce{scanNum} = fDataTfceCount/nRand;
          if strcmp(testParams.randomizationTestOutput,'Z value')
            fRandTfce{scanNum} = pToZ(fRandTfce{scanNum},1/nRand,precision); 
          end
        end
      end
    end
  end


  % save eventRelated parameters
  glmAnal.d{scanNum}.hrf = downsample(d.hrf, d.supersampling);
  glmAnal.d{scanNum}.hdrlen = d.hdrlen;
  glmAnal.d{scanNum}.actualhrf = d.hrf;
  glmAnal.d{scanNum}.trsupersampling = d.supersampling;
  glmAnal.d{scanNum}.ver = d.ver;
  glmAnal.d{scanNum}.filename = d.filename;
  glmAnal.d{scanNum}.filepath = d.filepath;
  glmAnal.d{scanNum}.dim = [scanDims numVolumes];
  glmAnal.d{scanNum}.nhdr = d.nhdr;
  glmAnal.d{scanNum}.tr = d.tr;
  glmAnal.d{scanNum}.stimNames = d.stimNames;
  glmAnal.d{scanNum}.EVnames = testParams.EVnames;
  glmAnal.d{scanNum}.stimToEVmatrix = testParams.stimToEVmatrix;
  glmAnal.d{scanNum}.scm = d.scm;
  glmAnal.d{scanNum}.expname = d.expname;
  glmAnal.d{scanNum}.fullpath = d.fullpath;
  glmAnal.d{scanNum}.ehdr = NaN([scanDims size(ehdr,4) size(ehdr,5)],precision);
  glmAnal.d{scanNum}.ehdr(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdr;
  glmAnal.d{scanNum}.ehdrste = NaN([scanDims size(ehdrste,4) size(ehdrste,5)],precision);
  glmAnal.d{scanNum}.ehdrste(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = ehdrste;
  glmAnal.d{scanNum}.rss = NaN(scanDims,precision);
  glmAnal.d{scanNum}.rss(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:,:) = rss; %JB
  clear('ehdr','ehdrste','rss');
  if ~isempty(contrasts)
    glmAnal.d{scanNum}.rdf = d.rdf; %JB
    glmAnal.d{scanNum}.contrasts = contrasts;
  end
  if ~isempty(fTests)
    glmAnal.d{scanNum}.fTests = fTests;
    glmAnal.d{scanNum}.rdf = d.rdf; %JB
    glmAnal.d{scanNum}.mdf = d.mdf; %JB
  end

  stimvol = d.stimvol;
  for i=1:length(stimvol)
    stimvol{i} = unique(ceil(stimvol{i}/d.supersampling));
  end
  glmAnal.d{scanNum}.stimvol = stimvol;

end
toc


%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
glmAnal.name = params.saveName;
if strcmp(params.hrfModel,'hrfDeconvolution')
  glmAnal.type = 'deconvAnal';
elseif isempty(fTests) && ~params.computeTtests
  glmAnal.type = 'glmAnal';
else
  glmAnal.type = 'glmAnalStats';
end
glmAnal.groupName = params.groupName;
glmAnal.function = 'eventRelatedGlm';
glmAnal.reconcileFunction = 'defaultReconcileParams';
glmAnal.mergeFunction = 'defaultMergeParams';
glmAnal.guiFunction = 'eventRelatedGlmGUI';
glmAnal.params = params;
glmAnal.date = dateString;
thisView = viewSet(thisView,'newAnalysis',glmAnal);

%--------------------------------------------------------- Output overlay structures
% create generic parameters 
overlay.groupName = params.groupName;
overlay.function = 'eventRelatedGlm';
overlay.reconcileFunction = 'defaultReconcileParams';
overlay.date = dateString;
overlay.params = cell(1,viewGet(thisView,'nScans'));
% colormap is made with a little bit less on the dark end
overlay.colormap = hot(312);
overlay.colormap = overlay.colormap(end-255:end,:);
overlay.alpha = 1;
if strcmp(params.hrfModel,'hrfDeconvolution')
  overlay.interrogator = 'eventRelatedPlot';
else
  overlay.interrogator = 'glmContrastPlot';
end
overlay.mergeFunction = 'defaultMergeParams';
overlay.colormapType = 'normal';
overlay.range = [0 1];
overlay.clip = [0 1];
overlay.alphaOverlay='';
overlay.alphaOverlayExponent=1;
overlay.data = cell(1,viewGet(thisView,'nScans'));
for scanNum = params.scanNum
   overlay.data{scanNum} = NaN(scanDims,precision); %to make values outside the box transparent
end


%------------------------------------------------------ save the r2 overlay
thisOverlay = overlay;
thisOverlay.name = 'r2';
thisOverlay.colormapType = 'setRangeToMax';
for scanNum = params.scanNum
   thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) = r2{scanNum};
   thisOverlay.params{scanNum} = scanParams{scanNum};
end
clear('r2');
thisView = viewSet(thisView,'newoverlay',thisOverlay);

% for output 
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
if nargout > 1
  for i = 1:length(d)
    glmAnal.d{i}.r2 = thisOverlay.data{i}; %(this kind of doesn't make sense anymore since we're interested in more than just r2)
  end
  % make d strucutre
  if length(glmAnal.d) == 1
    d = glmAnal.d{1};
  else
    d = glmAnal.d;
  end
end


if (testParams.parametricTests && strcmp(testParams.parametricTestOutput,'P value')) || ...
      (testParams.parametricTests && strcmp(testParams.randomizationTestOutput,'P value'))
  % colormap is made with a little bit less on the dark end and only for the 10 smallest values
  stat_colormap = flipud(hot(11));
  stat_colormap = stat_colormap(1:10,:);
end

%--------------------------------------------- save the contrast overlay(s)

if ~isempty(contrasts)
  %this is to mask the beta values by the probability/Z maps     
  betaAlphaOverlay = cell(length(testParams.EVnames),1);
  betaAlphaOverlayExponent = 1;
  if params.computeTtests
    if testParams.parametricTests
      thisOverlay = overlay;
      switch(testParams.parametricTestOutput)
        case 'T/F value'                                                %T maps
          thisOverlay.colormap = jet(256);
          max_abs_T = max(max(max(max(max(abs(cell2mat(tData)))))));
          thisOverlay.range = [-max_abs_T max_abs_T];
          namePrefix = 'T ( ';
          betaAlphaOverlayExponent = 0;
        case 'P value'                                                  %statistical maps
          thisOverlay = overlay;
          thisOverlay.colormap = stat_colormap;
          namePrefix = 'p (T ';
          betaAlphaOverlayExponent = -1;      %set the alpha overlay to the statistical value with an inverse masking
        case 'Z value'                                                  %Z maps
          thisOverlay = overlay;
          %thisOverlay.range = [0 max(max(max(max(max(cell2mat(tData))))))];
          thisOverlay.range = [0 8.2096];
          namePrefix = 'Z (T ';
          betaAlphaOverlayExponent = .5;      %set the alpha overlay to the Z value with masking
      end

      thisOverlay.clip = thisOverlay.range;
      for iContrast = 1:size(contrasts,1)
        thisOverlay.name = [namePrefix testParams.EVnames{iContrast} ')'];
        if betaAlphaOverlayExponent
          betaAlphaOverlay{iContrast} = thisOverlay.name;
        end
        for scanNum = params.scanNum
          thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             tData{scanNum}(:,:,:,iContrast);
          thisOverlay.params{scanNum} = scanParams{scanNum};
        end
        thisView = viewSet(thisView,'newoverlay',thisOverlay);
      end

      clear('tData');

      if testParams.TFCE
         % T TFCE overlay
        thisOverlay = overlay;
        thisOverlay.range(1) = min(min(min(min(min(cell2mat(tDataTFCE))))));
        thisOverlay.range(2) = min(max(max(max(max(cell2mat(tDataTFCE))))));
        thisOverlay.clip = thisOverlay.range;
        for iContrast = 1:size(contrasts,1)
           thisOverlay.name = ['T TFCE (' testParams.EVnames{iContrast} ')'];
           if testParams.randomizationTests
              thisOverlay.name = [thisOverlay.name ' - threshold(p=' num2str(randAlpha) ') = ' num2str(thresholdTfceT(iContrast,scanNum))];
           end
           for scanNum = params.scanNum
              thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
                 tDataTFCE{scanNum}(:,:,:,iContrast);
              thisOverlay.params{scanNum} = scanParams{scanNum};
           end
           thisView = viewSet(thisView,'newoverlay',thisOverlay);
        end
        clear('tDataTFCE');
      end
    end
    
    if testParams.randomizationTests
      thisOverlay = overlay;
      switch(testParams.randomizationTestOutput)
        case 'P value'                                                  %statistical maps
          thisOverlay = overlay;
          thisOverlay.colormap = stat_colormap;
          namePrefix = 'p rand';
          betaAlphaOverlayExponent = -1;      %set the alpha overlay to the statistical value with an inverse masking
        case 'Z value'                                                  %Z maps
          thisOverlay = overlay;
          %thisOverlay.range = [0 max(max(max(max(max(cell2mat(contrastRandData))))))];
          thisOverlay.range = [0 8.2096];
          namePrefix = 'Z rand';
          betaAlphaOverlayExponent = .5;      %set the alpha overlay to the Z value
      end
      thisOverlay.clip = thisOverlay.range;
      for iContrast = 1:size(contrasts,1)
        thisOverlay.name = [namePrefix ' (T ' testParams.EVnames{iContrast} ')'];
        if betaAlphaOverlayExponent
          betaAlphaOverlay{iContrast} = thisOverlay.name;
        end
        for scanNum = params.scanNum
          thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             contrastRandData{scanNum}(:,:,:,iContrast);
          thisOverlay.params{scanNum} = scanParams{scanNum};
        end
        thisView = viewSet(thisView,'newoverlay',thisOverlay);
      end
      clear('contrastRandData');

      if testParams.TFCE
        thisOverlay = overlay;
        switch(testParams.randomizationTestOutput)
          case 'P value'                                                  %statistical maps
            thisOverlay = overlay;
            thisOverlay.colormap = stat_colormap;
            namePrefix = 'p rand';
          case 'Z value'                                                  %Z maps
            thisOverlay = overlay;
            %thisOverlay.range = [0 max(max(max(max(max(cell2mat(tRandTfce))))))];
            thisOverlay.range = [0 8.2096];
            namePrefix = 'Z rand';
        end

        thisOverlay.clip = thisOverlay.range;
        for iFtest = 1:size(contrasts,1)
          thisOverlay.name = [namePrefix ' (T TFCE ' testParams.EVnames{iFtest} ')'];
          for scanNum = params.scanNum
            thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
               tRandTfce{scanNum}(:,:,:,iFtest);
            thisOverlay.params{scanNum} = scanParams{scanNum};
          end
          thisView = viewSet(thisView,'newoverlay',thisOverlay);
        end
        clear('tRandTfce');
      end

    end
  end
%--------------------------------------------- save the contrast beta weights overlay(s) (ehdr if no contrast)

  if length(testParams.componentsToTest)==1 || strcmp(testParams.componentsCombination,'Add')
    %find the values for the scale of the beta overlays
    thisOverlay = overlay;
    ordered_abs_betas = cell2mat(contrastData);
    ordered_abs_betas = ordered_abs_betas(~isnan(ordered_abs_betas));
    ordered_abs_betas = sort(abs(ordered_abs_betas));
    beta_perc95 = 0; 
    beta_perc95 = max(beta_perc95,ordered_abs_betas(round(numel(ordered_abs_betas)*.95))); %take the 95th percentile for the min/max
    thisOverlay.range = [-beta_perc95 beta_perc95];
    min_beta = min(min(min(min(min(cell2mat(contrastData))))));
    max_beta = max(max(max(max(max(cell2mat(contrastData))))));
    thisOverlay.clip = [min_beta max_beta];
    thisOverlay.colormap = jet(256);
    for iContrast = 1:length(testParams.EVnames)
       thisOverlay.alphaOverlay=betaAlphaOverlay{iContrast};
       thisOverlay.name = testParams.EVnames{iContrast};
       for scanNum = params.scanNum
          thisOverlay.data{scanNum} = NaN(scanDims,precision); %to make values outside the box transparent
          thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             contrastData{scanNum}(:,:,:,iContrast);
          thisOverlay.params{scanNum} = scanParams{scanNum};
       end
       thisOverlay.alphaOverlayExponent=betaAlphaOverlayExponent;
       thisView = viewSet(thisView,'newoverlay',thisOverlay);
    end
    clear('contrastData');
  end
end

%----------------------------------------------- save the F-test overlay(s)
if ~isempty(fTests)
  if testParams.parametricTests
      thisOverlay = overlay;
    switch testParams.parametricTestOutput
      case 'T/F value'
        thisOverlay.range(1) = min(min(min(min(min(cell2mat(fData))))));
        thisOverlay.range(2) = min(max(max(max(max(cell2mat(fData))))));
        namePrefix = 'F ( ';
      case 'P value'
        thisOverlay = overlay;
        thisOverlay.colormap = stat_colormap;
        namePrefix = 'p (F ';
      case 'Z value'
        %thisOverlay.range = [0 max(max(max(max(cell2mat(fData)))))];
        thisOverlay.range = [0 8.2096];
        namePrefix = 'Z (F ';
    end
    thisOverlay.clip = thisOverlay.range;
    for iFtest = 1:size(fTests,1)
      %probability maps
      thisOverlay.name = [namePrefix num2str(fTests(iFtest,:)) ')'];
      for scanNum = params.scanNum
        thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
           fData{scanNum}(:,:,:,iFtest);
        thisOverlay.params{scanNum} = scanParams{scanNum};
      end
      thisView = viewSet(thisView,'newoverlay',thisOverlay);
    end
    clear('fData');

    if testParams.TFCE
       % F TFCE overlay
      thisOverlay = overlay;
      thisOverlay.range(1) = min(min(min(min(cell2mat(fDataTFCE)))));
      thisOverlay.range(2) = max(max(max(max(cell2mat(fDataTFCE)))));
      thisOverlay.clip = thisOverlay.range;
      for iFtest = 1:size(fTests,1)
         thisOverlay.name = ['F TFCE (' num2str(fTests(iFtest,:)) ')'];
         if testParams.randomizationTests
            thisOverlay.name = [thisOverlay.name ' - threshold(p=' num2str(randAlpha) ') = ' num2str(thresholdTfceF(iFtest,scanNum)) ')'];
         end
         for scanNum = params.scanNum
            thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) = ...
               fDataTFCE{scanNum}(:,:,:,iFtest);
            thisOverlay.params{scanNum} = scanParams{scanNum};
         end
         thisView = viewSet(thisView,'newoverlay',thisOverlay);
      end
      clear('fDataTFCE');
    end
  end
  
  if testParams.randomizationTests
    thisOverlay = overlay;
    switch(testParams.randomizationTestOutput)
      case 'P value'                                                  %statistical maps
        thisOverlay = overlay;
        thisOverlay.colormap = stat_colormap;
        namePrefix = 'p rand';
      case 'Z value'                                                  %Z maps
        thisOverlay = overlay;
        thisOverlay.range = [0 max(max(max(max(max(cell2mat(fRandData))))))];
        thisOverlay.range = [0 8.2096];
        namePrefix = 'Z rand';
    end

    thisOverlay.clip = thisOverlay.range;
    for iFtest = 1:size(fTests,1)
      thisOverlay.name = [namePrefix ' (F ' num2str(fTests(iFtest,:)) ')'];
      for scanNum = params.scanNum
        thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
           fRandData{scanNum}(:,:,:,iFtest);
        thisOverlay.params{scanNum} = scanParams{scanNum};
      end
      thisView = viewSet(thisView,'newoverlay',thisOverlay);
    end
    clear('fRandData');
    
    if testParams.TFCE
      thisOverlay = overlay;
      switch(testParams.randomizationTestOutput)
        case 'P value'                                                  %statistical maps
          thisOverlay = overlay;
          thisOverlay.colormap = stat_colormap;
          namePrefix = 'p rand';
        case 'Z value'                                                  %Z maps
          thisOverlay = overlay;
          %thisOverlay.range = [0 max(max(max(max(max(cell2mat(fRandTfce))))))];
          thisOverlay.range = [0 8.2096];
          namePrefix = 'Z rand';
      end

      thisOverlay.clip = thisOverlay.range;
      for iFtest = 1:size(fTests,1)
        thisOverlay.name = [namePrefix ' (F TFCE ' num2str(fTests(iFtest,:)) ')'];
        for scanNum = params.scanNum
          thisOverlay.data{scanNum}(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2)) =...
             fRandTfce{scanNum}(:,:,:,iFtest);
          thisOverlay.params{scanNum} = scanParams{scanNum};
        end
        thisView = viewSet(thisView,'newoverlay',thisOverlay);
      end
      clear('fRandTfce');
    end
    
  end
end

%-------------------------------------------------------- Save the analysis
saveAnalysis(thisView,glmAnal.name);

oneTimeWarning('nonZeroHrfStart',0);
oneTimeWarning('tfceOutputsZeros',0);
set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow
refreshMLRDisplay(viewGet(thisView,'viewNum'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sub routines

function outputData = reshapeToRoiBox(inputData,dataPosition)
%inputData must be size ([A 1 1 B], where A equals the number of non-zero values in dataPosition
%outputData will be size ([size(dataPosition) B])
outputData = NaN([numel(dataPosition) size(inputData,4) size(inputData,5)]);
inputData = permute(inputData,[1 4 5 2 3]); %remove singleton dimensions (y and z)
outputData(dataPosition>0,:,:) = inputData;
outputData = reshape(outputData,[size(dataPosition,1) size(dataPosition,2) size(dataPosition,3) size(outputData,2) size(outputData,3)]);


function Z = pToZ(p, minP, outputPrecision)
%convert P value to Z value, not allowing probabilities of 0 (replaced by minP)
Z = norminv(1-max(p,minP));  %we don't want Z values of +/-Inf
Z(Z<0) = 0; % we'renot interested in negative Z value
%if there was no round-off error from cdf, we could do the following:
%Z = max(-norminv(p),0);  %because the left side of norminv seems to be less sensitive to round-off errors,
%we get -norminv(x) instead of norminv(1-x). also we'renot interested in negative Z value
Z(isnan(p)) = NaN; %NaNs must remain NaNs (they became 0 when using max)
Z = eval([outputPrecision '(Z)']);


