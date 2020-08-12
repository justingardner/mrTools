
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
end
if fieldIsNotDefined(params,'averagingMode')
  params.averagingMode = 'marginal';  % options are 'marginal' or 'interaction'
end

if justGetParams, return; end


if strcmp(params.averagingMode,'interaction')
  if ~all(ismember(params.factors,commonFactors))
    mrErrorDlg('(mlrGroupAverage) Cannot compute averages because factors are missing in some scans');
  end
  factorNames{1} = params.factors;
else
  for iFactor = 1:length(params.factors)
    factorNames{iFactor} = params.factors(iFactor);
  end
end

thisView = viewSet(thisView,'curGroup',params.groupNum);

for iScan = 1:length(params.scanList)
  tseriesPath{iScan} = viewGet(thisView,'tseriespathstr',params.scanList(iScan));
  hdr{iScan} = cbiReadNiftiHeader(tseriesPath{iScan});
  for iFactor = 1:length(factorNames)
    levels{iScan} = [];
    for jFactor = 1:length(factorNames{iFactor})
      % check that the number of volumes matches the number of elements in the factor variables
      if length(factors{iScan}.(factorNames{iFactor}{jFactor})) ~= hdr{iScan}.dim(5)
        mrErrorDlg(sprintf('(mlrGroupAverage) Scan %d: Mismatched number of volumes between .mat variable ''%s'' (%d) and time series file (%d)', ...
                            params.scanList(iScan),factorNames{iFactor}{jFactor},length(factors{iScan}.(factorNames{iFactor}{jFactor})),hdr{iScan}.dim(5)));
      end
      if size(factors{iScan}.(factorNames{iFactor}{jFactor}),1)==1
        factors{iScan}.(factorNames{iFactor}{jFactor}) = factors{iScan}.(factorNames{iFactor}{jFactor})'; % make sure the factor is a column cell array
      end
      levels{iScan} = [levels{iScan} factors{iScan}.(factorNames{iFactor}{jFactor})];
    end
    if iFactor == 1 &&  iScan ==1
      uniqueLevels = unique(levels{iScan},'rows','stable');
    else
      uniqueLevels = union(uniqueLevels,levels{iScan},'rows','stable');
    end
  end
end

nLevels = size(uniqueLevels,1);
for iOverlay = 1:nLevels
  for iScan = 1:length(params.scanList)
    overlays(iOverlay).data{params.scanList(iScan)} = zeros(hdr{iScan}.dim(2:4)');
  end
end

minOverlay = inf(nLevels,1);
maxOverlay = -1*inf(nLevels,1);
for iScan = 1:length(params.scanList)
  volumeCount = zeros(nLevels,1);
  for iVolume = 1:hdr{iScan}.dim(5)
    [~,whichOverlay] = ismember(levels{iScan}(iVolume),uniqueLevels);
    % add data
    overlays(whichOverlay).data{params.scanList(iScan)} = overlays(whichOverlay).data{params.scanList(iScan)} + cbiReadNifti(tseriesPath{iScan},{[],[],[],iVolume},'double');
    volumeCount(whichOverlay) = volumeCount(whichOverlay) + 1;
  end
  for iOverlay = 1:nLevels
    overlays(iOverlay).data{params.scanList(iScan)} = cast(overlays(iOverlay).data{params.scanList(iScan)}/volumeCount(whichOverlay),mrGetPref('defaultPrecision'));
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
