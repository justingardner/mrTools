% pRFSplit.m
%
%       usage: pRF(v, scanNum, params, x,y,z,n, fit)
%          by: akshay jagadeesh
%        date: 06/20/2017
%     purpose: Split pRF data into chunks of voxels, so that they can be more easily computed.
%
%
function splits = pRFSplit(v, scanNum, params, x,y,z, n, fit, overlays)

%% Get current user and current session dir
curPath = pwd;
sherlockSessionPath = ['/share/PI/jlg/' curPath(findstr(curPath, 'data'):end)];
suid = params.pRFFit.suid;

% Set split directory and scripts directory
splitDir = 'Splits';
scriptsDir = 'Splits/Scripts';
numSplits = params.pRFFit.numSplits;
blockSize = ceil(n/numSplits);
whichSplit = 1;
scanDims = viewGet(v, 'scanDims');
pRFAnal = overlays.pRFAnal;

% Make Splits directory if it doesn't exist 
if exist(splitDir) ~= 7
  mkdir(splitDir);
end
if exist(scriptsDir) ~=7
  mkdir(scriptsDir);
end

system('echo "#\!/bin/bash" >! "Splits/Scripts/runAll.sh"');

for blockStart = 1:blockSize:n
  blockEnd = min(blockStart+blockSize-1,n);
  blockSize = blockEnd-blockStart+1;

  % Make a temporary ROI using the coordinates specified.
  loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
  loadROI.coords(1,1:blockSize) = x(blockStart:blockEnd);
  loadROI.coords(2,1:blockSize) = y(blockStart:blockEnd);
  loadROI.coords(3,1:blockSize) = z(blockStart:blockEnd);

  % Load time series
  loadROI = loadROITSeries(v, loadROI, scanNum, params.groupName);

  % reorder x,y,z coordinates since they can get scrambled in loadROITSeries
  x(blockStart:blockEnd) = loadROI.scanCoords(1,1:blockSize);
  y(blockStart:blockEnd) = loadROI.scanCoords(2,1:blockSize);
  z(blockStart:blockEnd) = loadROI.scanCoords(3,1:blockSize);
  % Keep the linear coords
  pRFAnal.d{scanNum}.linearCoords = [pRFAnal.d{scanNum}.linearCoords sub2ind(scanDims,x(blockStart:blockEnd),y(blockStart:blockEnd),z(blockStart:blockEnd))];

  % Add all the needed fields to a split struct
  %split.v = v;
  split.nVoxels = blockSize;
  split.scanCoords = loadROI.scanCoords;
  split.tSeries = loadROI.tSeries;
  split.stim = fit.stim;
  split.concatInfo = fit.concatInfo;
  %split.prefit = fit.prefit;
  split.pRFFitParams = params.pRFFit;
  split.paramsInfo = fit.paramsInfo;
  split.scanNum = scanNum;

  % Save split struct to the specified directory
  filename = sprintf('%s_split%d', params.saveName, whichSplit);
  saveFile = sprintf('%s/%s.mat', splitDir, filename);
  save(saveFile, 'split');
  disp(sprintf('Data split %d saved to %s', whichSplit, saveFile));

  % Call bash script to output a .sbatch file
  %disp('Generating bash scripts');
  %system(sprintf('sh ~/proj/mrTools/mrLoadRet/Plugin/pRF/Sherlock/generateBatchScripts.sh "%s" "%s" "%s" "%d"',params.saveName,sherlockSessionPath, suid, whichSplit));

  splits{whichSplit} = split;
  whichSplit = whichSplit+1;

end

% Save master split struct locally
prefit = fit.prefit;
disp('Saving master struct');
save(sprintf('Splits/%s_master.mat', params.saveName), 'fit', 'x', 'y', 'z', 'scanNum', 'overlays', 'pRFAnal', 'v', 'params', 'sherlockSessionPath', 'suid', 'prefit');

% Check if session directory exists on Sherlock - and make it otherwise.
[~,out] = system(sprintf('ssh %s@sherlock.stanford.edu "[ -d %s ] && echo exists || echo does not exist"', suid, sherlockSessionPath));
if ~strcmp(deblank(out), 'exists')
  disp('Session directory does not exist on Sherlock. Transferring session dir to Sherlock');
  system(sprintf('ssh %s@sherlock.stanford.edu "mkdir -p %s"', suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/Anatomy/* %s@sherlock.stanford.edu:%s/Anatomy/', curPath, suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/Etc/* %s@sherlock.stanford.edu:%s/Etc/', curPath, suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/%s/* %s@sherlock.stanford.edu:%s/%s/', curPath, params.groupName, suid, sherlockSessionPath, params.groupName)); 
  system(sprintf('rsync -q %s/mrSession.mat %s@sherlock.stanford.edu:%s/.', curPath, suid, sherlockSessionPath));
end

% Use rsync to transfer split structs to Sherlock
disp('Copying split structs (.mat) and batch scripts to sherlock server');
system(sprintf('rsync -q Splits/%s*.mat %s@sherlock.stanford.edu:%s/Splits/', params.saveName, suid, sherlockSessionPath));
system(sprintf('rsync -q %s/* %s@sherlock.stanford.edu:%s/%s', scriptsDir, suid, sherlockSessionPath, scriptsDir));

% Call batch submission scripts on Sherlock
%disp('Submitting batch scripts...');
%system(sprintf('ssh %s@sherlock.stanford.edu "cd %s/%s/; sh runAll.sh"', suid, sherlockSessionPath, scriptsDir));
