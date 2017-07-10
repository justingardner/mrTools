% pRFRunSplits.m
%
%      Usage: pRFRunSplits(varargin)
%         by: akshay jagadeesh
%       date: 06/26/2017
%
function splits = pRFRunSplits(splitFileName)

% Load the split file
l = load(sprintf('Splits/%s', splitFileName));
s = l.split;
% split contains: nVoxels, scanCoords, tSeries, stim, concatInfo, prefit, pRFFitParams

% Call pRFFit on the voxels

tic
for i = 1:s.nVoxels
  x = s.scanCoords(1,:); y = s.scanCoords(2,:); z = s.scanCoords(3,:);

  fit = pRFFit(s.v, s.scanNum, x(i),y(i),z(i), 'stim', s.stim, 'concatInfo', s.concatInfo, 'prefit', s.prefit, 'fitTypeParams', s.pRFFitParams, 'tSeries', s.tSeries(i,:)', 'paramsInfo', s.paramsInfo);

  if ~isempty(fit)
    thisR2(i) = fit.r2;
    thisPolarAngle(i) = fit.polarAngle;
    thisEccentricity(i) = fit.eccentricity;
    thisRfHalfWidth(i) = fit.std;
    rawParams(:,i) = fit.params(:);
    r(i,:) = fit.r;
  end

end
toc
splits.scanCoords = s.scanCoords;
splits.r2 = thisR2;
splits.polarAngle = thisPolarAngle;
splits.eccentricity = thisEccentricity;
splits.rfHalfWidth = thisRfHalfWidth;
splits.params = rawParams;
splits.r = r;
save(sprintf('Splits/Analysis/%s_Anal.mat', splitFileName(1:end-4)), 'splits');
