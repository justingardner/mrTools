% mlrNumWorkers.m
%
%        $Id:$ 
%      usage: n = mlrNumWorkers(<askToStart>)
%         by: justin gardner
%       date: 01/03/12
%    purpose: Checks how many workers are running for parallel toolbox.
%             If no workers are running and parallel toolbox is available
%             asks users if they want to start matlabpool (unless askToStart=0).
%             If no, parallel toolbox available then returns 1
%
function n = mlrNumWorkers(askToStart)

% check arguments
if ~any(nargin == [0 1])
  help mlrNumWorkers
  return
end

if nargin < 1
  global mlrNumWorkersAskToStart;
  askToStart = false;
  % check if the global is set to not ask anymore
  if ~isequal(mlrNumWorkersAskToStart,false)
    askToStart = true;
  end
  % next time we come through, we won't ask to start again
  mlrNumWorkersAskToStart = false;
end

% check existence of Parallel Computing Toolbox
toolboxVersions = ver;
n = 1;
if any(strcmp({toolboxVersions.Name},'Parallel Computing Toolbox'))
  % get the pool size
  if verLessThan('matlab','8.4')
    n = matlabpool('size');
  else
    pool = gcp;
    n = pool.NumWorkers;
  end
  if n == 0
    n = 1;
    if askToStart
      if askuser('(mlrNumWorkers) You can speed up performance by starting a matlab pool of workers using the parallel computing toolbox. Do you wish to start it now?',0,1)
	%see what cluster profiles exist
	clusterProfiles = parallel.clusterProfiles;
	if length(clusterProfiles) == 1
	  % only one possibility, so select that cluster profile
	  clusterProfile = clusterProfiles{1};
	else
	  % multiple profiles, so let user decide which one to use
	  paramsInfo{1} = {'matlabpoolType',clusterProfiles,'Choose what kind of matlab pool to open'};
	  params = mrParamsDialog(paramsInfo);
	  if isempty(params),return,end
	  clusterProfile = params.matlabpoolType;
	end
	if verLessThan('matlab','8.4')
	  matlabpool('open',clusterProfile);
	  % get number of processors
	  n = matlabpool('size');
	else
	  pool = parpool(clusterProfile);
	  % get number of processors
	  n = pool.NumWorkers
	end
      end
    else
      disp(sprintf('(mlrNumWorkers) You can speed up performance by starting a matlab pool of workers using the parallel computing toolbox. Type: matlabpool open'));
    end
  end
end
