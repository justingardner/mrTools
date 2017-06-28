% pRFSplit.m
%
%       usage: pRF(v, scanNum, params, x,y,z,n, fit)
%          by: akshay jagadeesh
%        date: 06/20/2017
%     purpose: Split pRF data into chunks of voxels, so that they can be more easily computed.
%
%
function splits = pRFSplit(v, scanNum, params, x,y,z, n, fit)

numSplits = params.pRFFit.numSplits;
blockSize = ceil(n/numSplits);
whichSplit = 1;

if exist('Splits') ~= 7
  mkdir('Splits')
end

for blockStart = 1:blockSize:n
  blockEnd = min(blockStart+blockSize-1,n);
  blockSize = blockEnd-blockStart+1;

  loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
  loadROI.coords(1,1:blockSize) = x(blockStart:blockEnd);
  loadROI.coords(2,1:blockSize) = y(blockStart:blockEnd);
  loadROI.coords(3,1:blockSize) = z(blockStart:blockEnd);

  loadROI = loadROITSeries(v, loadROI, scanNum, params.groupName);

  split.nVoxels = blockSize;
  split.scanCoords = loadROI.scanCoords;
  split.tSeries = loadROI.tSeries;
  split.stim = fit.stim;
  split.concatInfo = fit.concatInfo;
  split.prefit = fit.prefit;
  split.pRFFitParams = params.pRFFit;
  split.paramsInfo = fit.paramsInfo;
  split.scanNum = scanNum;

  save(sprintf('Splits/split%d.mat', whichSplit), 'split');
  disp(sprintf('Data split %d saved to Splits/split%d.mat', whichSplit, whichSplit));

  whichSplit = whichSplit+1;

  keyboard
end
