
%   [thisView,params] = mlrGroupAverage(thisView,params,<'justGetParams'>)
%
%   goal: averages volumes in "group" scan across subjects according to conditions
%         specified in .mat file linked to the scan.
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
%   params.outputSampleSize :  whether to output a map of the the voxelwise number of subjects entering in the average for
%                               each average overlay (non-NaN values in subject overlays) (default = false)
%
%   author: julien besle (10/08/2020)

function [thisView,params] = mlrGroupAverage(thisView,params,varargin)

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
if fieldIsNotDefined(params,'averagingMode')
  params.averagingMode = 'marginal';  % options are 'marginal' or 'interaction'
elseif ~ismember(params.averagingMode,{'marginal','interaction'})
  mrErrorDlg(sprintf('(mlrGroupAverage)Unknown combination mode ''%s''',params.averagingMode));
end
if fieldIsNotDefined(params,'outputSampleSize')
  params.outputSampleSize = false;
end

if justGetParams, return; end


if strcmp(params.averagingMode,'interaction')
  if ~all(ismember(params.factors,commonFactors))
    mrErrorDlg('(mlrGroupAverage) Cannot compute averages because factors are missing in some scans');
  end
  whichFactors = {1: length(params.factors)};
else
  whichFactors = num2cell(1:length(params.factors));
end

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


minOverlay = inf(sum(nLevels),1);
maxOverlay = -1*inf(sum(nLevels),1);
maxCount = 0;
cScan = 0;
for iScan = 1:viewGet(thisView,'nScans')
  if ismember(iScan,params.scanList)
    fprintf('(mlrGroupAverage) Computing averages for scan %d... ',iScan);
    cScan = cScan+1;
  end
  cOverlay = 0;
  for iFactor = 1:length(whichFactors)
    for iOverlay = 1:nLevels(iFactor) %for each overlay
      cOverlay = cOverlay + 1;
      if ismember(iScan,params.scanList)
        overlays(cOverlay).data{iScan} = zeros(hdr{cScan}.dim(2:4)'); % initialize scans with zeros
        volumeCount = zeros(hdr{cScan}.dim(2:4)');
        for iVolume = find(ismember(whichOverlay{cScan}(:,iFactor),cOverlay,'rows'))' %for each volume in the scan matching this (combination of) condition(s)
          data = cbiReadNifti(tseriesPath{cScan},{[],[],[],iVolume},'double'); % read the data
          isNotNaN = ~isnan(data);
          % add non-NaN values to the appropriate overlay(s)
          overlays(cOverlay).data{iScan}(isNotNaN) = ...
            overlays(cOverlay).data{iScan}(isNotNaN) + data(isNotNaN);
          volumeCount = volumeCount + isNotNaN;
        end
        % divide by the number of added overlays
        overlays(cOverlay).data{iScan} = cast(overlays(cOverlay).data{iScan}./volumeCount,mrGetPref('defaultPrecision'));
        % get min and max
        minOverlay(cOverlay) = min(minOverlay(cOverlay),min(overlays(cOverlay).data{iScan}(:)));
        maxOverlay(cOverlay) = max(maxOverlay(cOverlay),max(overlays(cOverlay).data{iScan}(:)));
        if params.outputSampleSize
          overlays(sum(nLevels)+cOverlay).data{iScan} = volumeCount;
          maxCount = max(maxCount,max(volumeCount(:)));
        end
      else
        overlays(cOverlay).data{iScan} = [];
        if params.outputSampleSize
          overlays(sum(nLevels)+cOverlay).data{iScan} = [];
        end
      end
    end
  end
  if ismember(iScan,params.scanList)
    fprintf('Done\n');
  end
end

%add overlays' missing fields
for iOverlay = 1:sum(nLevels)*(1+params.outputSampleSize)
  
  if iOverlay<=sum(nLevels)
    overlays(iOverlay).name = uniqueLevels{iOverlay};
    overlays(iOverlay).range = [minOverlay(iOverlay) maxOverlay(iOverlay)];
  else
    overlays(iOverlay).name = sprintf('N (%s)',uniqueLevels{iOverlay-sum(nLevels)});
    overlays(iOverlay).range = [0 maxCount];
  end
  overlays(iOverlay).groupName = viewGet(thisView,'groupName');
  overlays(iOverlay).params = params;
  overlays(iOverlay).type = 'Group average';
  overlays(iOverlay).function = 'mlrGroupAverage';
  overlays(iOverlay).interrogator = '';
  
  allScanData = []; % determine the 1st-99th percentile range
  for iScan = params.scanList
    allScanData = [allScanData;overlays(iOverlay).data{iScan}(~isnan(overlays(iOverlay).data{iScan}))];
  end
  allScanData = sort(allScanData);
  overlays(iOverlay).colorRange = allScanData(round([0.01 0.99]*numel(allScanData)))';
  
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
