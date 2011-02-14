% orange.m
%
%        $Id: orange.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = orange(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR orange
%

function colorMap = orange(numberColors)

hsvRed = rgb2hsv(color2RGB('orange'));
colorMap = repmat(hsvRed,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors/2+.5;
colorMap = hsv2rgb(colorMap);