function [newparams,newdata] = corAnalReconcileParams(groupName,params,data)
% params = corAnalReconcileParams(groupName,[params],[data])
%
% Checks for consistency between corAnal parameters and current status of
% group (nscans and tseries filenames). Reconciles params by matching
% tseries filenames. Also reorders data arrays along with params.
%
% params: Optional initial parameters (see corAnal).
% default if not passed: initialize new params with default values.
%
% data: cell array of data arrays
% default if not passed: then returns cell array of empty matrices
%
% djh 7/2006

groupNum = viewGet([],'groupNum',groupName);
nScans = viewGet([],'nscans',groupNum);

% Initialize newparams
newparams.groupName = groupName;
newparams.recompute = zeros(1,nScans);
newparams.ncycles = zeros(1,nScans);
newparams.detrend = cell(1,nScans);
newparams.spatialnorm = cell(1,nScans);
newparams.tseriesfile = cell(1,nScans);
for scan = 1:nScans
    newparams.detrend{scan} = 'None';
    newparams.spatialnorm{scan} = 'None';
    tseriesfile = viewGet([],'tseriesFile',scan,groupNum);
    newparams.tseriesfile{scan} = tseriesfile;
    newdata{scan} = [];
end

% Initialize newdata
if ieNotDefined('data')
    data = cell(1,nScans);
end
newdata = cell(1,nScans);

% Find scans with tseries files that match those specified in
% params.tseriesfile. Use only those scans and the corresponding params.
if ~ieNotDefined('params')
    for scan = 1:nScans
        tseriesfile = viewGet([],'tseriesFile',scan,groupNum);
        if strcmp(params.tseriesfile,'any')
            match = scan;
        else
            match = find(strcmp(tseriesfile,{params.tseriesfile{:}}));
        end
        if match
            newparams.recompute(scan) = params.recompute(match);
            newparams.ncycles(scan) = params.ncycles(match);
            newparams.detrend{scan} = params.detrend{match};
            newparams.spatialnorm{scan} = params.spatialnorm{match};
            newparams.tseriesfile{scan} = tseriesfile;
            newdata{scan} = data{match};
        end
    end
end
