
%   [thisView, params, uniqueLevels] = mlrGroupTtest(thisView,params,<'justGetParams'>)
%
%   goal: test contrast(s) against 0 across subjects according to conditions
%         specified in .mat file linked to a group scan, using a paired T-test.
%         A "group" scan is a scan in which each volume corresponds to some
%         single-subject estimate map for a given condition, concatenated across
%         multiple subjects and normalized to a common template space.
%         The mat file must contain at least one cell array of strings (factor) of length equal
%         to the number of volumes and specifying which volume corresponds to which
%         condition and/or subject.
%         Note: when subject-level OLS estimates are used, testing contrasts requires
%         an assumption of homoscedasticity of the within-subject variance. Also, the variance
%         maybe be overestimated, leading to conservative tests. However results by Mumford &
%         Nichols (Neuroimage 47, 2009) suggest that group-level T tests on OLS estimates are
%         fairly robust to violations of homoscdasticity and do not underperform too much
%
%   usage:
%     [~, params] = mlrGroupTtest(thisView, [],'justGetParams') %returns default parameters
%     ... % modify parameters
%     [~, params, uniqueLevels] = mlrGroupTtest(thisView,params,'justGetParams') % get levels
%     ... % modify params.contrasts
%     thisView = mlrSphericalNormGroup(thisView, params) % runs function with modified params and contrasts
%
%   Input parameters:
%   params.groupNum: group in which to run the analysis (default: current group)
%   params.analysisName: name of analysis where overlays will be saved (default: 'Group averages')
%   params.scanList: list of scans to process (default: all scans in group)
%   params.factors: factors whose levels will be used to code for the contrast
%   params.combinationMode: how factor levels will be combined. Options are 'marginal' (default) or 'interaction'
%                         'marginal': Each level of each factor will be computed
%                         'interaction': means corresponding to all level combinations across all factors will be computed
%   params.contrasts: each row specify which levels (or combination of levels) to include in the corresponding contrast, with
%                     each column corresponding to a given level or combination thereof, according to the combinationMode parameter
%                     using the order in which factors were specified, the order in which levels appear within each factor across scans
%   params.smoothingFWHM: FWHM of the Gaussian smoothing kernel in voxels of the specified base/scan space. Set to 0 for no smoothing (default).
%   params.smoothingSpace: space in which to smooth the data, 0 = current scan, any other number: base number (default: current base,
%                          unless smoothingFWHM>0, in which case, current base)
%                          This is only implemented for flat bases. Any other base (surface or volume) will be ignored and everything will be done in scan space
%   params.testSide: 'both', 'left' or 'right'
%   params.pThreshold: criterion for masking overlays, expressed as a probability value (default = 0.05)
%   params.testOutput: '-10log(P)' (default), 'P', 'Z'. P values are not corrected for multiple tests
%   params.fweAdjustment: default = false (uses default method in transformStatistic.m)
%   params.fdrAdjustment: default = true (uses default method in transformStatistic.m)
%   params.thresholdEstimates: whether to clip contrast estimate overlays using corresponding p-value (using most conservative correction)
%   params.outputContrastEstimates: output contrast estimates in addition to statistics (default = true)
%   params.outputContrastSte: output contrast standard errors in addition to statistics (default = false)
%   params.outputStatistic: output T statistic in addition to statistics (default = false)
%   params.outputSampleSize :  whether to output a map of the voxelwise number of subjects entering in the t-test for
%                               each average overlay (non-NaN values in subject overlays) (default = false)
%
%   Output: - contrast estimates
%           - P or T overlay are added to the specified analysis in the view 
%           - params structure (gives default parameters if 'justGetParams')
%           - uniqueLevels: list of unique levels/combinations of levels in the order to be used for contrast matrix,
%             according to given scan and parameter structure
%
%   author: julien besle (17/03/2021)

function [thisView, params, uniqueLevels] = mlrGroupTtest(thisView,params,varargin)

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end

if ieNotDefined('params')
  params = struct;
end

if fieldIsNotDefined(params,'groupNum')
  params.groupNum = viewGet(thisView,'curGroup');
end
nScans = viewGet(thisView,'nScans',params.groupNum);
if fieldIsNotDefined(params,'analysisName')
  params.analysisName = 'Group T-tests';
end
if fieldIsNotDefined(params,'scanList')
  params.scanList = 1:nScans;
end
if fieldIsNotDefined(params,'factors')
  params.factors = {};
end
if fieldIsNotDefined(params,'combinationMode')
  params.combinationMode = 'marginal';  % options are 'marginal' or 'interaction'
end
if fieldIsNotDefined(params,'contrasts')
  params.contrasts = [];
end
if fieldIsNotDefined(params,'smoothingFWHM')
  params.smoothingFWHM = 0;
end
if fieldIsNotDefined(params,'smoothingSpace')
  params.smoothingSpace = viewGet(thisView,'curbase');
end
if fieldIsNotDefined(params,'testSide')
  params.testSide = 'both';
end
if fieldIsNotDefined(params,'pThreshold')
  params.pThreshold = 0.05;
end
if fieldIsNotDefined(params,'testOutput')
  params.testOutput = '-log10(P)';
end
if fieldIsNotDefined(params,'fweAdjustment')
  params.fweAdjustment= false;
end
if fieldIsNotDefined(params,'fdrAdjustment')
  params.fdrAdjustment= true;
end
if fieldIsNotDefined(params,'thresholdEstimates')
  params.thresholdEstimates = true;
end
if fieldIsNotDefined(params,'outputStatistic')
  params.outputStatistic = false;
end
if fieldIsNotDefined(params,'outputContrastEstimates')
  params.outputContrastEstimates = true;
end
if fieldIsNotDefined(params,'outputContrastSte')
  params.outputContrastSte = false;
end
if fieldIsNotDefined(params,'outputSampleSize')
  params.outputSampleSize = false;
end

uniqueLevels = {};

%read log files associated with scans
cScan = 0;
noLinkedFile=false;
for iScan = params.scanList
  if iScan < 0 || iScan > nScans
    mrWarnDlg(sprintf('(mlrGroupAverage) Scan %d does not exist in group %d', iScan, params.groupNum));
    noLinkedFile = true;
  else
    cScan= cScan+1;
    logFileName = viewGet(thisView,'stimfilename',iScan,params.groupNum);
    if isempty(logFileName)
      mrWarnDlg(sprintf('(mlrGroupTtest) No mat file linked to scan %d, group %d', iScan, params.groupNum));
      noLinkedFile = true;
    else
      factors{cScan} = load(logFileName{1});
      if isempty(factors{cScan})
        mrWarnDlg(sprintf('(mlrGroupTtest) Cannot open file %s for scan %d, group %d', logFileName{1}, iScan, params.groupNum));
      end
      if iScan == params.scanList(1)
        commonFactors = fieldnames(factors{cScan});
        allFactors = commonFactors;
      else
        commonFactors = intersect(commonFactors,fieldnames(factors{cScan}));
        allFactors = union(allFactors,fieldnames(factors{cScan}));
      end
    end
  end
end
if noLinkedFile
  return;
end

if fieldIsNotDefined(params,'factors')
  params.factors = allFactors;
elseif ischar(params.factors)
  params.factors = {params.factors};
end
if ~ismember(params.combinationMode,{'marginal','interaction'})
  mrWarnDlg(sprintf('(mlrGroupTtest)Unknown combination mode ''%s''',params.combinationMode));
  return
end

if strcmp(params.combinationMode,'interaction')
  if ~all(ismember(params.factors,commonFactors))
    mrWarnDlg('(mlrGroupTtest) Cannot run t-tests because factors are missing in some scans');
    return;
  end
  whichFactors = {1: length(params.factors)};
else
  whichFactors = num2cell(1:length(params.factors));
end

if justGetParams && fieldIsNotDefined(params,'factors'), return; end

% CHECK THAT THERE IS EQUAL N FOR ALL LEVELS AND COMBINATION OF LEVELS HERE?

currentGroup = viewGet(thisView,'curGroup');
thisView = viewSet(thisView,'curGroup',params.groupNum);

compatibleLogfile = true;
for iScan = 1:length(params.scanList)
  tseriesPath{iScan} = viewGet(thisView,'tseriespathstr',params.scanList(iScan));
  hdr{iScan} = cbiReadNiftiHeader(tseriesPath{iScan});
  for iFactor = 1:length(params.factors)
    % check that the field exists for this scan
    if ~isfield(factors{iScan},params.factors{iFactor})
      mrWarnDlg(sprintf('(mlrGroupTtest) Variable ''%s'' does not exist in scan %d', params.factors{iFactor}, params.scanList(iScan)));
      compatibleLogfile = false;
    else
      % check that the number of volumes matches the number of elements in the factor variables
      if length(factors{iScan}.(params.factors{iFactor})) ~= hdr{iScan}.dim(5)
        mrWarnDlg(sprintf('(mlrGroupTtest) Scan %d: Mismatched number of volumes between .mat variable ''%s'' (%d) and time series file (%d)', ...
                            params.scanList(iScan),params.factors{iFactor},length(params.factors{iFactor}),hdr{iScan}.dim(5)));
        compatibleLogfile = false;
      end
      if size(factors{iScan}.(params.factors{iFactor}),1)==1
        factors{iScan}.(params.factors{iFactor}) = factors{iScan}.(params.factors{iFactor})'; % make sure the factor is a column cell array
      end
      levels{iScan}(:,iFactor) = factors{iScan}.(params.factors{iFactor});
    end
  end
  if iScan ==1
    allLevels = levels{iScan};
  else
    allLevels = [allLevels; levels{iScan}];
  end
end

if ~compatibleLogfile
  thisView = viewSet(thisView,'curGroup',currentGroup);
  return;
end

for iFactor = 1:length(params.factors)
  % get unique level numbers for each factor. This is necessary because unique.m with option 'rows'
  [~,~,allFactorLevelNums(:,iFactor)]= unique(allLevels(:,iFactor),'stable'); % does not support cell arrays
  % get corresponding unique level numbers for all volumes of each scan
  for iScan = 1:length(params.scanList)
    [~,levelNums{iScan}(:,iFactor)] = ismember(levels{iScan}(:,iFactor),unique(allLevels(:,iFactor),'stable'));
  end
end
nLevels = [];
nLevelsAcrossFactors = 0;
for iFactor = 1:length(whichFactors) % for each factor or combination of factors
  % count the unique levels or combination of levels
  [uniqueLevelNums,uniqueLevelIndices]=unique(allFactorLevelNums(:,whichFactors{iFactor}),'rows');
  % find the unique overlay number for each volume in each scan
  for iScan = 1:length(params.scanList)
    [~,whichOverlay{iScan}(:,iFactor)] = ismember(levelNums{iScan}(:,whichFactors{iFactor}),uniqueLevelNums,'rows');
    whichOverlay{iScan}(:,iFactor) = nLevelsAcrossFactors + whichOverlay{iScan}(:,iFactor);
  end
  % get corresponding unique level names
  for iLevel = 1:size(uniqueLevelNums,1)
    uniqueLevels{sum(nLevels)+iLevel} = [allLevels{uniqueLevelIndices(iLevel),whichFactors{iFactor}}];
  end
  nLevels(iFactor) = size(uniqueLevelNums,1);
  nLevelsAcrossFactors = nLevelsAcrossFactors + nLevels(iFactor);
end

if fieldIsNotDefined(params,'contrasts') || size(params.contrasts,2)~=sum(nLevels)
  params.contrasts = eye(sum(nLevels));
end

if justGetParams
  thisView = viewSet(thisView,'curGroup',currentGroup);
  return;
end

if all(params.smoothingFWHM==0) && params.smoothingSpace~=0
  mrWarnDlg('(mlrGroupTtest) Smoothing is set to 0, so all analyses will be done in current scan space');
  params.smoothingSpace = 0;
end
baseType = viewGet(thisView,'baseType',params.smoothingSpace);
if any(params.smoothingFWHM>0)
  if baseType==0
    mrWarnDlg('(mlrGroupTtest) Smoothing will be done in current scan space');
    params.smoothingSpace = 0;
  elseif baseType==2
    mrWarnDlg('(mlrGroupTtest) Smoothing not implemented for surfaces or vol, switching smoothing space to current scan');
    params.smoothingSpace = 0;
  end
end

%-------------------------------------- Compute contrasts and run t-tests
sampleSizesDiffer = false;
contrastNames = makeContrastNames(params.contrasts,uniqueLevels,params.testSide);
nContrasts = size(params.contrasts,1);
nOverlays = nContrasts* (1 + params.fweAdjustment + params.fdrAdjustment + ...
             params.outputStatistic + params.outputContrastEstimates + params.outputContrastSte + params.outputSampleSize);
minOverlay = inf(nOverlays,1);
maxOverlay = -1*minOverlay;
cScan = 0;
for iScan = 1:viewGet(thisView,'nScans')
  if ismember(iScan,params.scanList)
    cScan = cScan+1;
    
    % Check for equal Ns
    cLevel = 0;
    for iFactor = 1:length(whichFactors)
      for iOverlay = 1:nLevels(iFactor)
        cLevel = cLevel+1;
        sampleSize = nnz(ismember(whichOverlay{cScan}(:,iFactor),cLevel,'rows')); %for each volume (subject) in the scan matching this (combination of) levels(s)
        if cLevel == 1
          maxSampleSize = sampleSize;
        else
          if maxSampleSize~=sampleSize
            mrWarnDlg('(mlrGroupTtest) The number of subject overlays should be identical at all levels or combinations of levels');
            return;
          end
        end
      end
    end
    
    waitString = sprintf('(mlrGroupTtest) Running group-level t-tests for scan %d... ',iScan);
    if params.smoothingSpace > 0 && baseType == 1
      fprintf('%s\n',waitString);
    else
      hWaitBar = mrWaitBar(-inf,waitString);
    end
    
    % compute the contrast estimates, statistics and p values
    for iContrast = 1:nContrasts % for each contrast, need to read the data in
      if params.smoothingSpace == 0 || baseType ~= 1
        mrWaitBar( iContrast/nContrasts, hWaitBar);
      end
      nonZeroContrastLevels = find(params.contrasts(iContrast,:)~=0);
      subjectContrastEstimates = zeros(prod(hdr{cScan}.dim(2:4)),maxSampleSize);
      sampleSize = zeros(prod(hdr{cScan}.dim(2:4)),numel(nonZeroContrastLevels));
      cLevel = 0;
      dLevel = 0;
      for iFactor = 1:length(whichFactors)
        for iLevel = 1:nLevels(iFactor)
          cLevel = cLevel+1;
          if ismember(cLevel,nonZeroContrastLevels)
            dLevel = dLevel+1;
            volumes = find(ismember(whichOverlay{cScan}(:,iFactor),cLevel,'rows'))'; %for each volume in the scan matching this (combination of) condition(s)
            for iVolume = 1:length(volumes)
              data = cbiReadNifti(tseriesPath{cScan},{[],[],[],volumes(iVolume)},'double'); % read the data
              isNotNaN = ~isnan(data);
              % add non-NaN values to the appropriate overlay(s)
              subjectContrastEstimates(isNotNaN,iVolume) = subjectContrastEstimates(isNotNaN,iVolume) + params.contrasts(iContrast,cLevel)*data(isNotNaN);
              sampleSize(:,dLevel) = sampleSize(:,dLevel) + isNotNaN(:);
            end
          end
        end
      end
      if numel(nonZeroContrastLevels)>1
        for i = 2:dLevel
          if nnz(diff(sampleSize(:,[1 i]),1,2))
            sampleSizesDiffer = true;
          end
        end
        if sampleSizesDiffer
          keyboard % Ns are not equal across levels at all voxels
        end
        sampleSize = reshape(sampleSize(:,1),hdr{cScan}.dim(2:4)');
      end
      subjectContrastEstimates(subjectContrastEstimates==0)=NaN; % replace zeros by NaNs to avoid infinite values later on

      % need to reshape first
      sampleSize = reshape(sampleSize,hdr{cScan}.dim(2:4)');
      subjectContrastEstimates = reshape(subjectContrastEstimates,[hdr{cScan}.dim(2:4)' maxSampleSize]);

      % apply smoothing
      if any(params.smoothingFWHM>0)

        waitString = sprintf('(mlrGroupTtest) Smoothing data for scan %d',iScan);
        if params.smoothingSpace > 0 && baseType == 1
          fprintf('%s, contrast %d\n',waitString, iContrast);
          base2scan = viewGet(thisView,'base2scan',iScan,[],params.smoothingSpace);
          [subjectContrastEstimates, ~, baseCoordsMap] = getBaseSpaceOverlay(thisView, subjectContrastEstimates,iScan,params.smoothingSpace,'linear');
          %pre-compute coordinates map to put values back from flat base to scan space
          scanDims = viewGet(thisView,'dims',iScan);
          %make a coordinate map of which scan voxel each base map voxel corresponds to (convert base coordmap to scan coord map)
          flat2scan = inverseBaseCoordMap(baseCoordsMap,scanDims,base2scan);
        else
          mrWaitBar( iContrast/nContrasts, hWaitBar, waitString);
        end

        for iVolume = 1:maxSampleSize
          subjectContrastEstimates(:,:,:,iVolume) = spatialSmooth(subjectContrastEstimates(:,:,:,iVolume),params.smoothingFWHM);
        end

        if params.smoothingSpace > 0 && baseType == 1
          % transform data back into scan space
          subjectContrastEstimates = applyInverseBaseCoordMap(flat2scan,scanDims,subjectContrastEstimates);
        end
      end

      % compute contrast estimates and std error
      contrastEstimates = nansum(subjectContrastEstimates,4)./sampleSize;
      contrastSte = nansum((subjectContrastEstimates-repmat(contrastEstimates,[1 1 1 maxSampleSize])).^2,4) ./ ... % sum of squared errors
                               (sampleSize-1) ./ ... % divided by N-1
                               sqrt(sampleSize);    % divided by sqrt(N)
      %compute p
      t = contrastEstimates./contrastSte;
      p = nan(size(t));
      switch(params.testSide)
        case 'both' % two-tailed
          p(t>=0) = 2 * (1 - cdf('t', double(t(t>=0)), sampleSize(t>=0)-1)); %here use doubles to deal with small Ps
          p(t<0) = 2 * cdf('t', double(t(t<0)), sampleSize(t<0)-1);
        case 'left' % left tailed
          p = cdf('t', double(t), sampleSize-1);
        case 'right' % right-tailed
          p = 1 - cdf('t', double(t), sampleSize-1);
      end
      % we do not allow probabilities of 0 and replace them by minP
      % (this can occur because cdf cannot return values less than 1e-16)
      p(~isnan(p)) = max(p(~isnan(p)),1e-16);

      waitString = sprintf('(mlrGroupTtest) Correcting p-values for scan %d',iScan);
      if params.smoothingSpace > 0 && baseType == 1
        fprintf('%s, contrast %d\n',waitString, iContrast);
      else
        mrWaitBar( iContrast/nContrasts, hWaitBar, waitString);
      end
      outputPrecision = mrGetPref('defaultPrecision');
      [p, fdrAdjustedP, fweAdjustedP] = transformStatistic(p,outputPrecision, params);

      nOutputs=0;
      if params.outputContrastEstimates
        overlays(nOutputs*nContrasts+iContrast).data{iScan} = cast(contrastEstimates,outputPrecision);
        overlays(nOutputs*nContrasts+iContrast).name = contrastNames{iContrast};
        % get min and max
        minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(nOutputs*nContrasts+iContrast),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(nOutputs*nContrasts+iContrast),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        nOutputs = nOutputs+1;
        nOutputContrast = nOutputs;
      else
        nOutputContrast = 0;
      end

      if params.outputContrastSte
        overlays(nOutputs*nContrasts+iContrast).data{iScan} = cast(contrastSte,outputPrecision);
        overlays(nOutputs*nContrasts+iContrast).name = ['Std error: ' contrastNames{iContrast}];
        % get min and max
        minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(nOutputs*nContrasts+iContrast),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(nOutputs*nContrasts+iContrast),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        nOutputs = nOutputs+1;
      end

      overlays(nOutputs*nContrasts+iContrast).data{iScan} = p;
      overlays(nOutputs*nContrasts+iContrast).name = [params.testOutput ': ' contrastNames{iContrast}];
      % get min and max
      minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(nOutputs*nContrasts+iContrast),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
      maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(nOutputs*nContrasts+iContrast),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
      nOutputs = nOutputs+1;
      nOutputP = nOutputs;

      if params.fweAdjustment
        overlays(nOutputs*nContrasts+iContrast).data{iScan} = fweAdjustedP;
        overlays(nOutputs*nContrasts+iContrast).name = ['FWE-corrected ' params.testOutput ': ' contrastNames{iContrast}];
        % get min and max
        minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(nOutputs*nContrasts+iContrast),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(nOutputs*nContrasts+iContrast),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        nOutputs = nOutputs+1;
        nOutputFweP = nOutputs;
      else
        nOutputFweP = 0;
      end

      if params.fdrAdjustment
        overlays(nOutputs*nContrasts+iContrast).data{iScan} = fdrAdjustedP;
        overlays(nOutputs*nContrasts+iContrast).name = ['FDR-corrected ' params.testOutput ': ' contrastNames{iContrast}];
        % get min and max
        minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(iLevel),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(iLevel),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        nOutputs = nOutputs+1;
        nOutputFdrP = nOutputs;
      else
        nOutputFdrP = 0;
      end

      if params.outputStatistic
        overlays(nOutputs*nContrasts+iContrast).data{iScan} = cast(t,outputPrecision);
        overlays(nOutputs*nContrasts+iContrast).name = ['T: ' contrastNames{iContrast}];
        % get min and max
        minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(iLevel),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(iLevel),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        nOutputs = nOutputs+1;
      end

      if params.outputSampleSize
        if sampleSizesDiffer || iContrast == 1
          if ~sampleSizesDiffer
            overlays(nOutputs*nContrasts+iContrast).name = 'Sample size';
          else
            overlays(nOutputs*nContrasts+iContrast).name = ['N: ' contrastNames{iContrast}];
          end
          overlays(nOutputs*nContrasts+iContrast).data{iScan} = cast(sampleSize,outputPrecision);
          % get min and max
          minOverlay(nOutputs*nContrasts+iContrast) = min(minOverlay(iLevel),min(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
          maxOverlay(nOutputs*nContrasts+iContrast) = max(maxOverlay(iLevel),max(overlays(nOutputs*nContrasts+iContrast).data{iScan}(:)));
        end
        nOutputs = nOutputs+1;
        nOutputSampleSize = nOutputs;
      else
        nOutputSampleSize = 0;
      end
    end

    if params.smoothingSpace == 0 || baseType ~= 1
      mrCloseDlg(hWaitBar);
    else
     fprintf('(mlrGroupTtest) Scan %d done\n\n',iScan);
    end
  else
    for iOverlay = 1:nOverlays
      overlays(iOverlay).data{iScan} = [];
    end
  end
end

%add overlays' missing fields
switch(params.testOutput)
  case 'P'
    clipThreshold = [0 params.pThreshold];
    alphaOverlayExponent = -1;
  case 'Z'
    clipThreshold = [norminv(1-params.pThreshold) inf];
    alphaOverlayExponent = .5;
  case '-log10(P)'
    clipThreshold = [-log10(params.pThreshold) inf];
    alphaOverlayExponent = .5;
end
for iOutput = 1:nOutputs
  for iContrast = 1:nContrasts
    iOverlay = (iOutput-1)*nContrasts + iContrast;
    overlays(iOverlay).range = [minOverlay(iOverlay) maxOverlay(iOverlay)];
    overlays(iOverlay).groupName = viewGet(thisView,'groupName');
    overlays(iOverlay).params = params;
    overlays(iOverlay).type = 'Group t-test';
    overlays(iOverlay).function = 'mlrGroupTtest';
    overlays(iOverlay).interrogator = '';
    switch(iOutput)
      case {nOutputP,nOutputFweP,nOutputFdrP}
        overlays(iOverlay).clip(1) = max(overlays(iOverlay).range(1),clipThreshold(1));
        overlays(iOverlay).clip(2) = min(overlays(iOverlay).range(2),clipThreshold(2));
        switch params.testOutput
          case 'P'
            overlays.colormap = statsColorMap(256);
            overlays(iOverlay).colorRange = [0 1];
          case 'Z'
            overlays(iOverlay).colorRange = [0 norminv(1-1e-16)]; %1e-16 is smallest non-zero P value output by cdf in getGlmStatistics (local functions T2p and F2p)
          case '-log10(P)'
            overlays(iOverlay).colorRange = [0 -log10(1e-16)];
        end
      case nOutputContrast
        if params.thresholdEstimates
          if params.fweAdjustment
            overlays(iOverlay).alphaOverlay = overlays((nOutputFweP-1)*nContrasts + iContrast).name;
          elseif params.fdrAdjustment
            overlays(iOverlay).alphaOverlay = overlays((nOutputFdrP-1)*nContrasts + iContrast).name;
          else
            overlays(iOverlay).alphaOverlay = overlays((nOutputP-1)*nContrasts + iContrast).name;
          end
          overlays(iOverlay).alphaOverlayExponent = alphaOverlayExponent;
        end
    end
    if ~ismember(iOutput, [nOutputP nOutputFweP nOutputFdrP nOutputSampleSize]) % for overlays other and P-values and sample size
      allScanData = []; % determine the 1st-99th percentile range
      for iScan = params.scanList
        allScanData = [allScanData;overlays(iOverlay).data{iScan}(~isnan(overlays(iOverlay).data{iScan}))];
      end
      allScanData = sort(allScanData);
      overlays(iOverlay).colorRange = allScanData(round([0.01 0.99]*numel(allScanData)))';
    end
  end
end

% remove any overlay with no data (this happens when asking for sample sizes and they are identical at all levels)
overlays(isinf(minOverlay)) = [];

% set or create analysis
analysisNum = viewGet(thisView,'analysisNum',params.analysisName);
if isempty(analysisNum)
  thisView = newAnalysis(thisView,params.analysisName);
else
  thisView = viewSet(thisView,'curAnalysis',analysisNum);
end
% add overlays to view
thisView = viewSet(thisView,'newOverlay',overlays);
thisView = viewSet(thisView,'clipAcrossOverlays',false);
