% mlrNumWorkers.m
%
%        $Id:$ 
%      usage: n = mlrNumWorkers(<numWorkers>)
%         by: justin gardner
%       date: 01/03/12
%    purpose: Checks how many workers are running for parallel toolbox.
%             If numWorkers is set to a number, it will try to run that
%             many workers. If it is set to 1 then it will ask user to
%             choose
%
function n = mlrNumWorkers(numWorkers)

% check arguments
if ~any(nargin == [0 1])
  help mlrNumWorkers
  return
end

% ask the user what to do
askToStart = true;

% if no arguments then just start with maximum
if nargin < 1
  global mlrNumWorkersAskToStart;
  askToStart = false;
  % check if the global is set to not ask anymore
  if ~isequal(mlrNumWorkersAskToStart,false)
    askToStart = true;
  end
  % next time we come through, we won't ask to start again
  mlrNumWorkersAskToStart = false;
  numWorkers = 1;
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
  % if not already running
  if n == 0
    n = 1;
    if askToStart
      % get number of workers to use
      if any(numWorkers == [0 1]);
	paramsInfo = {};
	if exist('parcluster') == 2
	  % check local cluster for maximum workers
	  myCluster = parcluster('local');
	  maxWorkers = myCluster.NumWorkers;
	  % put up dialog to select number of workers
	  paramsInfo{end+1} = {'numWorkers',min(max(maxWorkers-2,2),maxWorkers),'type=numeric','minmax',[1 maxWorkers],'incdec=[-1 1]','Set to how many workers you want to use for parallel computing'};
	else
	  paramsInfo{end+1} = {'numWorkers',4,'type=numeric','minmax=[1 inf]','incdec=[-1 1]','Set to how many workers you want to use for parallel computing'};
	end
	% put up dialog box
	if numWorkers
	  params = mrParamsDialog(paramsInfo,'Start a pool of workers?');
	  if isempty(params),return,end
	else
	  % use default parameters
	  params = mrParamsDefault(paramsInfo);
	end
	numWorkers = params.numWorkers;
      end

      % start the pool
      if verLessThan('matlab','8.4')
	matlabpool(numWorkers);
	% get number of processors
	n = matlabpool('size');
      else
	pool = parpool(numWorkers);
	% get number of processors
	n = pool.NumWorkers
      end
    else
      disp(sprintf('(mlrNumWorkers) You can speed up performance by starting a matlab pool of workers using the parallel computing toolbox. Type: matlabpool open'));
    end
  end
end
