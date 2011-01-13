% resampleFWE.m
%
%        $Id$
%      usage: count = resampleFWE(resampled,actual,indexSortActual,params,d)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: 
%

function count = resampleFWE(resampled,actual,indexSortActual,indexReorderActual,actualIsNotNaN,params,d)

switch(params.resampleFWEadjustment)
  case 'Single Step'
    count = double(repmat(max(resampled,[],1),[prod(d.dim(1:3)) 1])>actual);
    
  case 'Step-down'
    %sort test statistic using provided sort index (sorting of the non-nan actual data)
    sortedResampled = NaN(size(indexSortActual));
    for i=1:size(resampled,2)
      sortedResampled(:,i) = resampled(indexSortActual(:,i),i);
%       sortedActual(:,i) = actual(indexSortActual(:,i),i);                           % DEBUG
    end
    %compute consecutive maxima
    for i = 2:size(sortedResampled,1)
      sortedResampled(i,:) = max(sortedResampled(i,:),sortedResampled(i-1,:));
    end
    count = zeros(size(actual));
%     reorderedActual = NaN(size(actual));                                            % DEBUG
    for i=1:size(count,2)
      count(actualIsNotNaN,i) = sortedResampled(indexReorderActual(:,i),i);
%       reorderedActual(actualIsNotNaN,i) = sortedActual(indexReorderActual(:,i),i);  % DEBUG
    end
    count = double(count>actual);
    count(~actualIsNotNaN,:)=NaN;
end