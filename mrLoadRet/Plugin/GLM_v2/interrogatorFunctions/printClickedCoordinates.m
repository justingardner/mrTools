% dummyInterrogator
%
%      usage: [  ] = printClickedCoordinates(thisView,overlayNum,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2023-10-06
%
%    purpose: print the coordinates of the mouse-clicked point in the main mrLoadRet window in various coordinate systems
%             this is especially useful for surfaces, where coordinates are not displayed in the mrLoadRet window when hovering the mouse

function printClickedCoordinates(thisView,overlayNum,scanNum,x,y,z,roi)

scanCoords = viewGet(thisView,'mouseDownScanCoords');
baseCoords = viewGet(thisView,'mouseDownBaseCoords');
talCoords = viewGet(thisView,'mouseDownTalCoords');
mniCoords = viewGet(thisView,'mouseDownMniCoords');

fprintf('Clicked coordinates: ')
if any(~isnan(scanCoords))
  fprintf('\tScan: %s', num2str(scanCoords,'%d '))
end
if any(~isnan(baseCoords))
  fprintf('\tBase: %s', num2str(baseCoords,'%d '))
end
if any(~isnan(talCoords))
  fprintf('\tTalairach: %s', num2str(talCoords,'%.1f '))
end
if any(~isnan(mniCoords))
  fprintf('\tMNI: %s', num2str(mniCoords,'%.1f '))
end
fprintf('\n')

return;
