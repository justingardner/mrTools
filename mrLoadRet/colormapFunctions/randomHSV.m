% randomHSV.m
%
%      usage: colorMap = randomHSV(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap with colors randomly drawn from the HSV color map (without replacement)
%

function colorMap = randomHSV(numberColors)

colorMap = hsv(numberColors);
colorMap = colorMap(randperm(numberColors,numberColors),:);