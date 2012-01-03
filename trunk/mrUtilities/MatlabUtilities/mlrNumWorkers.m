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
  askToStart = true;
end

% check existence of Parallel Computing Toolbox
toolboxVersions = ver;
n = 1;
if any(strcmp({toolboxVersions.Name},'Parallel Computing Toolbox'))
  n = matlabpool('size');
  if n == 0
    n = 1;
    if askToStart
      if askuser('(mlrNumWorkers) You can speed up performance by starting a matlab pool of workers using the parallel computing toolbox. Do you wish to start it now',0,1)
	matlabpool('open');
	n = matlabpool('size');
      end
    else
      disp(sprintf('(mlrNumWorkers) You can speed up performance by starting a matlab pool of workers using the parallel computing toolbox. Type: matlabpool open'));
    end
  end
end
