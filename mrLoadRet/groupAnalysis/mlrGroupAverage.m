
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
if fieldIsNotDefined(params,'analysisName')
  params.analysisName = 'Group averages';
end
if fieldIsNotDefined(params,'scanList')
  params.scanList = 1:viewGet(thisView,'nScans');
end

%read log files associated with scans
cScan = 0;
for iScan = params.scanList
  cScan= cScan+1;
  logFileName = viewGet(thisView,'stimfilename',iScan);
  if isempty(logFileName{1})
    mrWarnDlg(sprintf('(mlrGroupAverage) No mat file linked to scan %d', iScan));
  end
  factors{cScan} = load(logFileName{1});
  if isempty(factors{cScan})
    mrWarnDlg(sprintf('(mlrGroupAverage) Cannot open file %s for scan %d', logFileName, iScan));
  end
  if iScan == 1
    commonFactors = fieldnames(factors{cScan});
    allFactors = commonFactors;
  else
    commonFactors = intersect(commonFactors,fieldnames(factors{cScan}));
    allFactors = union(allFactors,fieldnames(factors{cScan}));
  end
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

if justGetParams, return; end


if strcmp(params.averagingMode,'interaction')
  if ~all(ismember(params.factors,commonFactors))
    mrErrorDlg('(mlrGroupAverage) Cannot compute averages because factors are missing in some scans');
  end
  whichFactors = {1: length(params.factors)};
else
  whichFactors = num2cell(1:length(params.factors));
end

thisView = viewSet(thisView,'curGroup',params.groupNum);

for iScan = 1:length(params.scanList)
  tseriesPath{iScan} = viewGet(thisView,'tseriespathstr',params.scanList(iScan));
  hdr{iScan} = cbiReadNiftiHeader(tseriesPath{iScan});
  for iFactor = 1:length(params.factors)
    % check that the field exists for this scan
    if ~isfield(factors{iScan},params.factors{iFactor})
      mrErrorDlg(sprintf('(mlrGroupAverage) Variable ''%s'' does not exist in scan %d', params.factors{iFactor}, params.scanList(iScan)));
    end
    % check that the number of volumes matches the number of elements in the factor variables
    if length(factors{iScan}.(params.factors{iFactor})) ~= hdr{iScan}.dim(5)
      mrErrorDlg(sprintf('(mlrGroupAverage) Scan %d: Mismatched number of volumes between .mat variable ''%s'' (%d) and time series file (%d)', ...
                          params.scanList(iScan),params.factors{iFactor},length(params.factors{iFactor}),hdr{iScan}.dim(5)));
    end
    if size(factors{iScan}.(params.factors{iFactor}),1)==1
      factors{iScan}.(params.factors{iFactor}) = factors{iScan}.(params.factors{iFactor})'; % make sure the factor is a column cell array
    end
    levels{iScan}(:,iFactor) = factors{iScan}.(params.factors{iFactor});
  end
  if iScan ==1
    allLevels = levels{iScan};
  else
    allLevels = [allLevels; levels{iScan}];
  end
end

for iFactor = 1:length(params.factors)
  % get unique level numbers for each factor. This is necessary because unique.m with option 'rows'
  [~,~,allFactorLevelNums(:,iFactor)]= unique(allLevels(:,iFactor),'stable'); % do not support cell arrays
  % get corresponding unique level numbers for all volumes of each scan
  for iScan = 1:length(params.scanList)
    [~,levelNums{iScan}(:,iFactor)] = ismember(levels{iScan}(:,iFactor),unique(allLevels(:,iFactor),'stable'));
  end
end
nLevels = 0;
for iFactor = 1:length(whichFactors) % for each factor or combination of factors
  % count the unique levels or combination of levels
  [uniqueLevelNums,uniqueLevelIndices]=unique(allFactorLevelNums(:,whichFactors{iFactor}),'rows');
  % find the unique overlay number for each volume in each scan
  for iScan = 1:length(params.scanList)
    [~,~,whichOverlay{iScan}(:,iFactor)] = unique(levelNums{iScan}(:,whichFactors{iFactor}),'rows');
    whichOverlay{iScan}(:,iFactor) = whichOverlay{iScan}(:,iFactor) + nLevels;
  end
  % get corresponding unique level names
  for iLevel = 1:size(uniqueLevelNums,1)
    uniqueLevels{nLevels+iLevel} = [allLevels{uniqueLevelIndices(iLevel),whichFactors{iFactor}}];
  end
  nLevels = nLevels + size(uniqueLevelNums,1);
end

% initialize scans with zeros
for iOverlay = 1:nLevels
  for iScan = 1:length(params.scanList)
    overlays(iOverlay).data{params.scanList(iScan)} = zeros(hdr{iScan}.dim(2:4)');
  end
end

minOverlay = inf(nLevels,1);
maxOverlay = -1*inf(nLevels,1);
for iScan = 1:length(params.scanList)
  volumeCount = zeros(nLevels,1);
  for iVolume = 1:hdr{iScan}.dim(5) %for each volume in the scan
    data = cbiReadNifti(tseriesPath{iScan},{[],[],[],iVolume},'double'); % read the data
    for iFactor = 1:length(whichFactors) % add it to the appropriate overlay(s)
      overlays(whichOverlay{iScan}(iVolume,iFactor)).data{params.scanList(iScan)} = overlays(whichOverlay{iScan}(iVolume,iFactor)).data{params.scanList(iScan)} + data;
      volumeCount(whichOverlay{iScan}(iVolume,iFactor)) = volumeCount(whichOverlay{iScan}(iVolume,iFactor)) + 1;
    end
  end
  for iOverlay = 1:nLevels %for each overlay
    % divide by the number of added overlays
    overlays(iOverlay).data{params.scanList(iScan)} = cast(overlays(iOverlay).data{params.scanList(iScan)}/volumeCount(iOverlay),mrGetPref('defaultPrecision'));
    % get min and max
    minOverlay(iOverlay) = min(minOverlay(iOverlay),min(overlays(iOverlay).data{params.scanList(iScan)}(:)));
    maxOverlay(iOverlay) = max(maxOverlay(iOverlay),max(overlays(iOverlay).data{params.scanList(iScan)}(:)));
  end
end

%add overlays' missing fields
for iOverlay = 1:nLevels
  overlays(iOverlay).name = uniqueLevels{iOverlay};
  overlays(iOverlay).groupName = viewGet(thisView,'groupName');
  overlays(iOverlay).params = params;
  overlays(iOverlay).range = [minOverlay(iOverlay) maxOverlay(iOverlay)];
  overlays(iOverlay).type = 'Group average';
  overlays(iOverlay).function = 'mlrGroupAverage';
  overlays(iOverlay).interrogator = '';
end

% set or create analysis
analysisNum = viewGet(thisView,'analysisName',params.analysisName);
if isempty(viewGet(thisView,'analysisNum',analysisNum))
  thisView = newAnalysis(thisView,params.analysisName);
else
  thisView = viewSet(thisView,'curAnalysis',analysisNum);
end
% add overlays to view
thisView = viewSet(thisView,'newOverlay',overlays);
