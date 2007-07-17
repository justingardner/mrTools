function newparams = motionCompReconcileParams(groupName,params)
% params = motionCompReconcileParams(groupName,[params])
%
% Checks for consistency between motionComp parameters and current
% status of group (nscans and tseries filenames). Reconciles params by
% matching tseries filenames.
%
% params: Optional initial parameters (see motionComp,
% motionCompBetweenScans, and motionCompWithinScan).
%
% djh 7/2006

groupNum = viewGet([],'groupNum',groupName);
nScans = viewGet([],'nscans',groupNum);

% Get tseries filenames for this group
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
  oldDescription = viewGet([],'description',scan,groupNum);
  descriptions{scan} = ['motion comp of ',groupName,' scan ',num2str(scan) ': ' oldDescription];
end

if ieNotDefined('params')
  % Use default params
  newparams.groupName = groupName;
  newparams.baseScan = 1;
  newparams.baseFrame = 'first';
  newparams.sliceTimeCorrection = 0;
  newparams.robust = 0;
  newparams.correctIntensityContrast = 0;
  newparams.crop = [];
  newparams.niters = 3;
  newparams.motionCompGroupName = 'MotionComp';
  newparams.interpMethod = 'nearest';
  newparams.targetScans = [1:nScans];
  newparams.tseriesfiles = tseriesfiles;
  newparams.descriptions = descriptions;
else
  % Set newparams according to params, reconciling with tseries files.
  newparams.groupName = params.groupName;
  newparams.baseScan = params.baseScan;
  newparams.baseFrame = params.baseFrame;
  newparams.sliceTimeCorrection = params.sliceTimeCorrection;
  newparams.robust = params.robust;
  newparams.correctIntensityContrast = params.correctIntensityContrast;
  newparams.crop = params.crop;
  newparams.niters = params.niters;
  newparams.motionCompGroupName = params.motionCompGroupName;
  newparams.interpMethod = params.interpMethod;
  % Find scans with tseries files that match those specified in
  % params.tseriesfiles. Use only those scans and the corresponding
  % params.
  if strcmp(params.tseriesfiles,'any')
    params.tseriesfiles = tseriesfiles;
    params.descriptions = descriptions;
  end
  match = [];
  targetScans = [];
  for s = params.targetScans
    m = find(strcmp(tseriesfiles{s},{params.tseriesfiles{:}}));
    if m
      targetScans = [targetScans,s];
      match = [match,m];
    end
  end
  newparams.targetScans = targetScans;
  newparams.tseriesfiles = {params.tseriesfiles{match}};
  newparams.descriptions = {params.descriptions{match}};
end
