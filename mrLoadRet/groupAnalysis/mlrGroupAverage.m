
%   [thisView, params, uniqueLevels] = mlrGroupAverage(thisView,params,<'justGetParams'>)
%
%   goal: averages volumes in "group" scan across subjects according to conditions
%         (or combination of conditions) specified in .mat file linked to the scan.
%         Optionally, computes averages only for a subset of (combinations of) conditions,
%         or linear combinations of these averages (using contrasts parameter)
%         A "group" scan is a scan in which each volume corresponds to some
%         single-subject estimate map for a given condition, concatenated across
%         multiple subjects and normalized to a common template space.
%         The mat file must contain at least one cell array of strings (factor) of length equal
%         to the number of volumes and specifying which volume corresponds to which
%         condition and/or subject. Multiple cell arrays of equal length can be used
%         to describe more complex factorial designs and calculate averages for combinations
%         of conditions.
%
%   usage:
%     [~,params] = mlrGroupAverage(thisView, [],'justGetParams') %returns default parameters
%     ... % modify parameters
%     [~, params, uniqueLevels] = mlrGroupAverage(thisView,params,'justGetParams') % get levels (optional)
%     ... % modify params.contrasts
%     thisView = mlrSphericalNormGroup(thisView, params) % runs function with modified params
%
%   parameters:
%   params.groupNum: group in which to run the analysis (default: current group)
%   params.analysisName: name of analysis where overlays will be saved (default: 'Group averages')
%   params.scanList: list of scans to process (default: all scans in group)
%   params.factors: factors whose levels will be used to average (default: all factors found in mat files associated with scans
%   params.averagingMode: how levels will be combined. Options are 'marginal' (default) or 'interaction'
%                         'marginal': all marginal means corresponding to each level of each facotr will be computed
%                         'interaction': means corresponding to all level combinations across all factors will be computed
%   params.contrasts: each row specifies which levels (or combination of levels) to include in the corresponding contrast, with
%                     each column corresponding to a given level or combination thereof, according to the combinationMode parameter
%                     using the order in which factors were specified, the order in which levels appear within each factor across scans
%   params.outputSampleSize :  whether to output a map of the the voxelwise number of subjects entering in the average for
%                               each average overlay (non-NaN values in subject overlays) (default = false)
%
%   author: julien besle (10/08/2020)

function [thisView, params, uniqueLevels] = mlrGroupAverage(thisView,params,varargin)

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
  params.analysisName = 'Group averages';
end
if fieldIsNotDefined(params,'scanList')
  params.scanList = 1:nScans;
end
if fieldIsNotDefined(params,'factors')
  params.factors = {};
end
if fieldIsNotDefined(params,'averagingMode')
  params.averagingMode = 'marginal';  % options are 'marginal' or 'interaction'
end
if fieldIsNotDefined(params,'contrasts')
  params.contrasts = [];
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
if ~ismember(params.averagingMode,{'marginal','interaction'})
  mrWarnDlg(sprintf('(mlrGroupAverage)Unknown combination mode ''%s''',params.averagingMode));
  return;
end

if strcmp(params.averagingMode,'interaction')
  if ~all(ismember(params.factors,commonFactors))
    mrWarnDlg('(mlrGroupAverage) Cannot compute averages because factors are missing in some scans');
    return;
  end
  whichFactors = {1: length(params.factors)};
else
  whichFactors = num2cell(1:length(params.factors));
end

if justGetParams && fieldIsNotDefined(params,'factors'), return; end

currentGroup = viewGet(thisView,'curGroup');
thisView = viewSet(thisView,'curGroup',params.groupNum);

compatibleLogfile = true;
for iScan = 1:length(params.scanList)
  tseriesPath{iScan} = viewGet(thisView,'tseriespathstr',params.scanList(iScan));
  hdr{iScan} = cbiReadNiftiHeader(tseriesPath{iScan});
  for iFactor = 1:length(params.factors)
    % check that the field exists for this scan
    if ~isfield(factors{iScan},params.factors{iFactor})
      mrWarnDlg(sprintf('(mlrGroupAverage) Variable ''%s'' does not exist in scan %d', params.factors{iFactor}, params.scanList(iScan)));
      compatibleLogfile = false;
    else
      % check that the number of volumes matches the number of elements in the factor variables
      if length(factors{iScan}.(params.factors{iFactor})) ~= hdr{iScan}.dim(5)
        mrWarnDlg(sprintf('(mlrGroupAverage) Scan %d: Mismatched number of volumes between .mat variable ''%s'' (%d) and time series file (%d)', ...
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


%-------------------------------------- Compute averages (or contrasts)
contrastNames = makeContrastNames(params.contrasts,uniqueLevels,'no test');
nContrasts = size(params.contrasts,1);
nonZeroContrastLevels = find(any(params.contrasts~=0));
nnzContrastLevels = numel(nonZeroContrastLevels);
minOverlay = inf(nContrasts+nnzContrastLevels,1);
maxOverlay = -1*minOverlay;
cScan = 0;
outputPrecision = mrGetPref('defaultPrecision');
for iScan = 1:viewGet(thisView,'nScans')
  if ismember(iScan,params.scanList)
    hWaitBar = mrWaitBar(-inf,sprintf('(mlrGroupAverage) Computing averages for scan %d... ',iScan));

    cScan = cScan+1;
  
    levelsData = zeros(prod(hdr{cScan}.dim(2:4)),nnzContrastLevels);
    sampleSize = zeros(prod(hdr{cScan}.dim(2:4)),nnzContrastLevels);
    cLevel = 0;
    dLevel = 0;
    for iFactor = 1:length(whichFactors)
      for iLevel = 1:nLevels(iFactor) %for each overlay
        cLevel = cLevel + 1;
        if ismember(cLevel,nonZeroContrastLevels)
          dLevel = dLevel+1;
          mrWaitBar( dLevel/nnzContrastLevels, hWaitBar);
          for iVolume = find(ismember(whichOverlay{cScan}(:,iFactor),cLevel,'rows'))' %for each volume in the scan matching this (combination of) condition(s)
            data = cbiReadNifti(tseriesPath{cScan},{[],[],[],iVolume},'double'); % read the data
            isNotNaN = ~isnan(data);
            % add non-NaN values to the appropriate overlay(s)
            levelsData(isNotNaN,dLevel) = levelsData(isNotNaN,dLevel) + data(isNotNaN);
            sampleSize(:,dLevel) = sampleSize(:,dLevel) + isNotNaN(:);
          end
          % divide by the number of added overlays
          levelsData(:,dLevel) = levelsData(:,dLevel)./sampleSize(:,dLevel);
        end
      end
    end
    
  end

  for iContrast = 1:nContrasts
    if ismember(iScan,params.scanList)
      nonZeroLevels = find(params.contrasts(iContrast,nonZeroContrastLevels)); % need to only use levels involved in this particular contrast to avoid NaNs in the other levels 
      overlays(iContrast).data{iScan} = cast(reshape(levelsData(:,nonZeroLevels)*params.contrasts(iContrast,nonZeroContrastLevels(nonZeroLevels))',...
                                                                 hdr{cScan}.dim(2:4)'),outputPrecision);
      overlays(iContrast).name = contrastNames{iContrast};
      % get min and max
      minOverlay(iContrast) = min(minOverlay(iContrast),min(overlays(iContrast).data{iScan}(:)));
      maxOverlay(iContrast) = max(maxOverlay(iContrast),max(overlays(iContrast).data{iScan}(:)));
    else
      overlays(iContrast).data{iScan} = [];
    end
  end
  
  if params.outputSampleSize
    for iLevel = 1:nnzContrastLevels
      if ismember(iScan,params.scanList)
        overlays(nContrasts+iLevel).data{iScan} = cast(reshape(sampleSize(:,iLevel),hdr{cScan}.dim(2:4)'),outputPrecision);
        overlays(nContrasts+iLevel).name = ['N: ' uniqueLevels{nonZeroContrastLevels(iLevel)}];
        % get min and max
        minOverlay(nContrasts+iLevel) = min(minOverlay(nContrasts+iLevel),min(overlays(nContrasts+iLevel).data{iScan}(:)));
        maxOverlay(nContrasts+iLevel) = max(maxOverlay(nContrasts+iLevel),max(overlays(nContrasts+iLevel).data{iScan}(:)));
      else
        overlays(nContrasts+iLevel).data{iScan} = [];
      end
    end
  end

  if ismember(iScan,params.scanList)
    mrCloseDlg(hWaitBar);
  end

end

%add overlays' missing fields
for iOutput = 1:nContrasts + params.outputSampleSize*nnzContrastLevels
  overlays(iOutput).range = [minOverlay(iOutput) maxOverlay(iOutput)];
  overlays(iOutput).groupName = viewGet(thisView,'groupName');
  overlays(iOutput).params = params;
  overlays(iOutput).type = 'Group average';
  overlays(iOutput).function = 'mlrGroupAverage';
  overlays(iOutput).interrogator = '';

  if iOutput<=nContrasts
    allScanData = []; % determine the 1st-99th percentile range
    for iScan = params.scanList
      allScanData = [allScanData;overlays(iOutput).data{iScan}(~isnan(overlays(iOutput).data{iScan}))];
    end
    allScanData = sort(allScanData);
    overlays(iOutput).colorRange = allScanData(round([0.01 0.99]*numel(allScanData)))';
  end
end

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
