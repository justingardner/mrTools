% estimateTrueH0Ratio.m
%
%        $Id$
%      usage: adjustedP = estimateTrueH0Ratio(p,params,)
%         by: julien besle, 
%       date: 17/01/2011
%    purpose: estimates the proportion of null hypotheses (see Sarkar, 2009, submitted)
%

function [trueH0Ratio,lambda] = estimateTrueH0Ratio(p, params)

isNotNan = ~isnan(p);
%sizePdata = size(p);
p = p(isNotNan);
[p,sortingIndex] = sort(p);
numberH0 = length(p);


switch(params.adaptiveFweMethod)
  case 'Threshold'
    %Storey et al. (2004) Sarkar, 2009, submitted
    trueH0Ratio = (numberH0 - nnz(p<params.adaptiveFweThreshold)+1)/numberH0/(1-params.adaptiveFweThreshold);
    trueH0Ratio = min(trueH0Ratio,1); %this ratio could be >1 if no test is rejected at the threshold
                                  %in this case, we set the ratio to 1 and use the non-adaptive procedure
    lambda = params.adaptiveFweThreshold;
    
  case 'Adaptive'
    %Schweder and Spjotvoll (1982), Hochberg and Benjamini (1990), Benjamini and Hochberg, (2000)
    %compute the estimated number of true H0 when considering the number of rejections at all possible levels
    numberTrueH0 = (numberH0+1-(1:numberH0)')./(1-p);
    %find the first estimated number that is greater than the previous one
    k = find(diff(numberTrueH0)>0,1,'first');
    numberTrueH0 = ceil(min(numberTrueH0(k+1),numberH0));
    lambda = p(k);
    %recompute the adjusted p-value using the estimated ratio of true H0
    trueH0Ratio = numberTrueH0/numberH0;
    
end
