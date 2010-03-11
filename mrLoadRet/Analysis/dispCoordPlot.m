%
%      usage: dispCoordPlot()
%         by: justin gardner
%       date: 11/08/07
%    purpose: interrogator function that displays desired voxel in
%             a different view
%
function retval = dispCoordPlot(v,overlayNum,scan,x,y,s,roi,varargin)

% get the view number
viewNum = viewGet(v,'viewNum');

% convert scan coordinates to base coordinates
scanCoords = [x y s 1]';
scan2base = inv(viewGet(v,'base2scan'));
baseCoords = round(scan2base*scanCoords);

% get magnet coordinates
scan2mag = inv(viewGet(v,'scan2mag'));
magCoords = round(scan2mag*scanCoords);

% display the coordinates
disp(sprintf('(dispCoordPlot) scanCoords = [%i %i %i]',x,y,s));
disp(sprintf('(dispCoordPlot) baseCoords = [%i %i %i]',baseCoords(1),baseCoords(2),baseCoords(3)));
disp(sprintf('(dispCoordPlot) magCoords = [%i %i %i]',magCoords(1),magCoords(2),magCoords(3)));

