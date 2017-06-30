% pRFSplit.m
%
%       usage: pRF(v, scanNum, params, x,y,z,n, fit)
%          by: akshay jagadeesh
%        date: 06/20/2017
%     purpose: Split pRF data into chunks of voxels, so that they can be more easily computed.
%
%
function splits = pRFSplit(v, scanNum, params, x,y,z, n, fit)

splitDir = '~/Google Drive/data/mgldoublebars/s036020170331/Splits';
splitDir = './Splits';
numSplits = params.pRFFit.numSplits;
blockSize = ceil(n/numSplits);
whichSplit = 1;

if exist(splitDir) ~= 7
  mkdir(splitDir);
end
if exist('Splits/Scripts') ~=7
  mkdir('Splits/Scripts');
end

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

  % Add all the needed fields to a split struct
  split.v = v;
  split.nVoxels = blockSize;
  split.scanCoords = loadROI.scanCoords;
  split.tSeries = loadROI.tSeries;
  split.stim = fit.stim;
  split.concatInfo = fit.concatInfo;
  split.prefit = fit.prefit;
  split.pRFFitParams = params.pRFFit;
  split.paramsInfo = fit.paramsInfo;
  split.scanNum = scanNum;

  % Save split struct to the specified directory
  filename = sprintf('%ssplit%d', params.saveName, whichSplit);
  saveFile = sprintf('%s/%s.mat', splitDir, filename);
  save(saveFile, 'split');
  disp(sprintf('Data split %d saved to %s', whichSplit, saveFile));

  % Call bash script to output a .sbatch file
  disp('Generating bash scripts');
  system(sprintf('sh ~/proj/mrTools/mrLoadRet/Plugin/pRF/generateBatchScripts.sh "%s"',filename));

  splits{whichSplit} = split;
  whichSplit = whichSplit+1;

end
disp('Copying split structs (.mat) and batch scripts to sherlock server');
system('rsync -q Splits/* akshayj@sherlock.stanford.edu:/share/PI/jlg/data/mgldoublebars/s036020170331/Splits/');
system('rsync -q Splits/Scripts/* akshayj@sherlock.stanford.edu:/share/PI/jlg/data/mgldoublebars/s036020170331/Splits/Scripts/');
