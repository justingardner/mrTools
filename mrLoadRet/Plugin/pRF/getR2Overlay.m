function r2Overlay = getR2Overlay(modelName)

v = newView;
v = viewSet(v, 'curGroup', 'Averages');
v = viewSet(v, 'curScan', 1);
prfModel = sprintf('pRFAnal/%s', modelName);
v = loadAnalysis(v, prfModel);
r2OverlayNum = find(strcmp('r2', viewGet(v, 'overlayNames')));
v = viewSet(v, 'curOverlay', r2OverlayNum);
r2Overlay = viewGet(v, 'overlayData', 1);

validR2 = r2Overlay(~isnan(r2Overlay));
fprintf('Mean r2 for well-defined voxels: %d', mean(validR2));
