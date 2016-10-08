% trainPRF.m
%
%        $Id:$ 
%      usage: trainPRF(view, roiName)
%         by: akshay jagadeesh
%       date: 10/04/16
%    purpose: Given n MotionComp runs, this function computes a n-fold 
%             cross validation, running the pRF analysis n times. 
%             It generates and saves time series, model response, residuals,
%             and covariance matrices for each fold.
%

function trainPRF(view, roiName)

% Currently hardcoding number and names of scan.
numScans = 6
roiName = 'xVal'

for i = 1:numScans
  v = newView;

  %Get all but the i'th scan
  scanList = [1:i-1 i+1:numScans]

  %Average the other n-1 scans together
  v = viewSet(v, 'currentGroup', 'MotionComp');
  scanstr = sprintf('scanList=%s', mat2str(scanList));
  [v params] = averageTSeries(v, [], 'justGetParams=1', 'defaultParams=1', scanstr)
  %Filename of averaged scan is: tseries-avg-<SCANLIST>.nii
  scanListStr = mat2str(scanList); scanListStr = scanListStr(2:end-1); scanListStr=scanListStr(find(~isspace(scanListStr)));
  %params.fileName = sprintf('tseries-avg-%s.nii', fileName)
  %params = averageTSeriesGUI('groupName', 'MotionComp');
  averageTSeries(v, params)

  %Get group info about the average we just created
  gInfo = groupInfo;
  gAverageInfo = groupInfo('Averages');
  numAverageScans = gInfo(3).numScans;
  sInfo = scanInfo(numAverageScans, 'Averages');
      %% sInfo.description --> contains info about which scans were averaged
      %% sInfo.Filename --> contains filename in which this average is saved
      %% sInfo.OriginalN --> contains filename of original (replace N with 1 - 5)

  v = newView;
  v = viewSet(v, 'currentGroup', 'Averages');

  % Two ways of calculating ROI coordinates
  % Method 1: loadROI.coords
  v = loadROI(v, sprintf('%s.mat', roiName))
  roiCoords = v.ROIs(1).coords.'
  roiCoords = roiCoords(:, 1:3)

  % Method 2: loadROITSeries.scancoords
  roiCoords2 = loadROITSeries(v, 'xVal')
  roiCoords2 = roiCoords2.scanCoords.'

  %Run pRF Analysis on average of n-1 scans
  %pRF(v, params)
  [v, params] = pRF(v, [], 'justGetParams=1');
      %%set dispStimScan to the last most recently computed scan
      %%set saveName to pRF_ROINAME_SCANLIST
      %%select diffOfGamma
  pRF(v, params)

  % Load analysis by filename
  analysisFileName = sprintf('pRF_%s_%s_.mat', roiName, scanListStr)
  v = loadAnalysis(v, ['pRFAnal/' analysisFileName])
  scanNum = v.analyses{1}.params.scanNum;
 
  %Get the residual and covariance matrix
         %%% scanNum: scan number on which analysis was conducted
  [residual, covMat, tSeries, modelResponse] = pRFNoise(v, scanNum, coordinates, analysisFileName)

  %%
  %% TO DO: Get the time series data for the left-out scan
  %%

  %Save residual, covariance matrix, time series, and model response into a structure
  save(sprintf('SaveData/pRF_fold%i_%s', i, scanListStr), 'residual', 'covMat', 'tSeries', 'modelResponse')
  
end


function test
