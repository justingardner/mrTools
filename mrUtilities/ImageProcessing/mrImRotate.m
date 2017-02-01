% mlrImRotate.m
%
%      usage: B = mrImRotate(A, rotate, interpMethod,	bbox, extrapVal)
%         by: julien besle
%       date: 23/75/2016
%    purpose: rotates matrix/image A by 'rotate' degrees
%             (replacement for imrotate from the Image Processing toolbox)
%             gives identical results to imrotate for interpMethod = nearest neighbour
%             and bbox='loose', and very similar results for other options
%             Also, added option extrapVal

function B = mrImRotate(A, rotate, interpMethod, bbox, extrapVal)

if ieNotDefined('bbox')
  bbox = 'loose';
end
if ieNotDefined('interpMethod')
  interpMethod = 'nearest';
end
if ieNotDefined('extrapVal')
  extrapVal = 0;
end
switch(bbox)
  case 'bilinear'
    bbox = 'linear';
  case 'bicubic'
    bbox = 'cubic';
end

originalSizeA = size(A);
if strcmp(bbox,'loose')
  %pad symmetrically around the image
  margin = ceil((max(originalSizeA)*sqrt(2)-originalSizeA)/2);
  paddedA = extrapVal*ones(margin*2 +originalSizeA);
  paddedA(margin(1)+1:margin(1)+originalSizeA(1),margin(2)+1:margin(2)+originalSizeA(2))=A;
  A=paddedA;
end

[x1,y1] = meshgrid(1:size(A,2),1:size(A,1));
x0 = (size(A,2)+1)/2;
y0 = (size(A,1)+1)/2;

theta = rotate/180*pi;
x2=cos(theta)*(x1-x0) - sin(theta)*(y1-y0) + x0;
y2=sin(theta)*(x1-x0) + cos(theta)*(y1-y0) + y0;

if strcmp(bbox,'loose')
  originalx2 = x2-margin(2);
  originaly2 = y2-margin(1);
  outsidePixels = originalx2<0.5 | originalx2>originalSizeA(2)+0.5 | ...
                  originaly2<0.5 | originaly2>originalSizeA(1)+0.5;
  columnsToKeep = any(~outsidePixels,1);
  rowsToKeep = any(~outsidePixels,2);
  x2 = x2(rowsToKeep,columnsToKeep);
  y2 = y2(rowsToKeep,columnsToKeep);
end

B = interp2(x1,y1,A,x2,y2,interpMethod,extrapVal);

