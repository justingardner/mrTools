% pRFRunSplits.m
%
%      Usage: pRFRunSplits(varargin)
%         by: akshay jagadeesh
%       date: 06/26/2017
%
function splits = pRFRunSplits(splitFileName)

% Load the split file
l = load(splitFileName);
s = l.split;
% split contains: nVoxels, scanCoords, tSeries, stim, concatInfo, prefit, pRFFitParams



% Call pRFFit on the voxels

for i = 1:s.nVoxels
  x = scanCoords(1,:); y = scanCoords(2,:); z = scanCoords(3,:);

  fit = pRFFit(v, s.scanNum, x(i),y(i),z(i), 'stim', s.stim, 'concatInfo', s.concatInfo, 'prefit', s.prefit, 'fitTypeParams', s.pRFFitParams, 'tSeries', s.tSeries(i,:), 'paramsInfo', s.paramsInfo);

  if ~isempty(fit)
    thisR2(i) = fit.r2;
    thisPolarAngle(i) = fit.polarAngle;
    thisEccentricity(i) = fit.eccentricity;
    thisRfHalfWidth(i) = fit.std;
    rawParams(:,i) = fit.params(:);
    r(i,:) = fit.r;
  end

end
