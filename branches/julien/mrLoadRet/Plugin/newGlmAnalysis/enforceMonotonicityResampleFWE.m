% enforceMonotonictyResampleFWE.m
%
%        $Id$
%      usage: count = enforceMonotonictyResampleFWE(p,indexSortActual,indexReorderActual,actualIsNotNaN,params)
%         by: julien besle, 
%       date: 14/01/2011
%    purpose: enforces ordered p-value monotonicity in Step-Down bootstrap-based FWE control
%             see algorithm 4.1, step 6, Westfall and Young (1993) p117

function p = enforceMonotonicityResampleFWE(p,indexSortActual,indexReorderActual,actualIsNotNaN,params)

if strcmp(params.resampleFWEadjustment,'Step-down')
  %sort p using provided sort index (sorting of the non-nan actual statistics)
  %note that these values are sorted from the smallest to the largest
  sortedP = NaN(size(indexSortActual));
  for i=1:size(p,2)
    sortedP(:,i) = p(indexSortActual(:,i),i);
%       sortedActual(:,i) = actual(indexSortActual(:,i),i);                           % DEBUG
  end
  %compute consecutive maxima (from largest to smallest p value)
  for i = size(sortedP,1)-1:-1:1
    sortedP(i,:) = max(sortedP(i,:),sortedP(i+1,:));
  end
  for i=1:size(p,2)
    p(actualIsNotNaN,i) = sortedP(indexReorderActual(:,i),i);
%       reorderedActual(actualIsNotNaN,i) = sortedActual(indexReorderActual(:,i),i);  % DEBUG
  end
  %p(~actualIsNotNaN,:)=NaN;
end