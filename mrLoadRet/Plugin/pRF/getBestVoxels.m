function cords = getBestVoxels(scanNum, bestN, roiName, analysis)

if(ieNotDefined('scanNum'))
  disp(sprintf('(getBestVoxels) No inputs given, using preset defaults'));
  scanNum = 1;
  bestN = 90;
  roiName = 'both_pRF';
  analysis = 'pRF.mat'
end

v = newView;
v = viewSet(v, 'currentGroup', 'Concatenation');
groupNum = viewGet(v, 'currentGroup');
%v = viewSet(v, 'currentGroup', 'Averages');
lV1 = loadROITSeries(v, roiName, scanNum, groupNum, 'straightXform=1', 'loadType=none');

v = loadAnalysis(v, ['pRFAnal/' analysis]);

r2 = viewGet(v, 'overlaydata', scanNum, 1, 1);
lV1 = getSortIndex(v, lV1, r2);

sortIndex = lV1{1}.sortindex;
r2 = lV1{1}.r2;
scanLinCoords = lV1{1}.scanLinearCoords;

scanDims = viewGet(v, 'scanDims');

[i,j,k] = ind2sub(scanDims, scanLinCoords(sortIndex(1:bestN)));
cords = [i;j;k;ones(1,bestN)];

%fits = pRF_crossValidate(cords);
%fit.cords = cords;
end
