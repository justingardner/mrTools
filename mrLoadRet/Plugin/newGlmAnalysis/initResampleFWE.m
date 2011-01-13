% initResampleFWE.m
%
%        $Id$
%      usage: [indexSort,indexReorder,isNotNaN] = initResampleFWE(actual,params,isNotNaN)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: 
%

function [indexSort,indexReorder,isNotNaN] = initResampleFWE(actual,params,isNotNaN)

if strcmp(params.resampleFWEadjustment,'Step-down')
  %find the sorting index for actual non-nan T values (independently for each contrast)
  
  if ieNotDefined('isNotNaN')
    isNotNaN = ~any(isnan(actual),2); %if some voxels are nan only for some tests, they are excluded from this analysis. But that should not happen often, if ever
  end
  actual(~isNotNaN,:) = NaN;
  [sortedActual,indexSort] = sort(actual,1);
  indexSort = indexSort(~any(isnan(sortedActual),2),:); 
  [temp,indexReorder] = sort(indexSort,1);
else
  indexSort = [];
  indexReorder = [];
  isNotNaN = [];
end
