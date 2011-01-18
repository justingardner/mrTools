% initResampleFWE.m
%
%        $Id$
%      usage: [indexSortActual,indexReorderActual,d.actualIsNotNaN] = initResampleFWE(actual,params,d.actualIsNotNaN)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: 
%

function d = initResampleFWE(actual,params,rdf,mdf)


d.actualIsNotNaN = ~isnan(actual); 
d.numberFalseH0 = zeros(1,size(actual,2));
d.numberTrueH0 = sum(d.actualIsNotNaN,1);

if strcmp(params.resampleFWEadjustment,'Step-down') || ~strcmp(params.adaptiveFweMethod,'None')
  %find the sorting index for actual non-nan T values (independently for each contrast)
  
  actual(~d.actualIsNotNaN) = NaN;
  [sortedActual,d.indexSortActual] = sort(actual,1);
  d.indexSortActual = mat2cell(d.indexSortActual,size(d.indexSortActual,1),ones(size(d.indexSortActual,2),1));
  d.indexReorderActual = cell(size(d.indexSortActual));
  for iTest = 1:length(d.indexSortActual)
    d.indexSortActual{iTest} = d.indexSortActual{iTest}(~isnan(sortedActual(:,iTest))); 
   [temp,d.indexReorderActual{iTest}] = sort(d.indexSortActual{iTest},1);
  end
else
  d.indexSortActual = {};
  d.indexReorderActual = {};
  d.actualIsNotNaN = [];
end


if ~strcmp(params.adaptiveFweMethod,'None')
  if nargin==3 %this is a T-test
    %convert to P value
    p = 1 - cdf('t', double(actual), rdf); %here use doubles to deal with small Ps
    if strcmp(params.tTestSide,'Both')
      p = 2*p;
    end
  elseif nargin==4 %this is an F-test
    p = 1 - cdf('f', double(actual), repmat(mdf,size(actual,1),1), repmat(rdf,size(actual)));  
  end
  for iTest = 1:size(p,2)
    trueH0ratio = estimateTrueH0Ratio(p(:,iTest),params);
    d.numberFalseH0(iTest) = round((1-trueH0ratio)*length(d.indexSortActual{iTest}));
  end
  d.numberTrueH0 = d.numberTrueH0 - d.numberFalseH0;
end


% %old version
% if strcmp(params.resampleFWEadjustment,'Step-down')
%   %find the sorting index for actual non-nan T values (independently for each contrast)
%   
%   if ieNotDefined('d.actualIsNotNaN')
%     d.actualIsNotNaN = ~any(isnan(actual),2); %if some voxels are nan only for some tests, they are excluded from this analysis. But that should not happen often, if ever
%   end
%   actual(~d.actualIsNotNaN,:) = NaN;
%   [sortedActual,indexSortActual] = sort(actual,1);
%   indexSortActual = indexSortActual(~any(isnan(sortedActual),2),:); 
%   [temp,indexReorderActual] = sort(indexSortActual,1);
% else
%   indexSortActual = [];
%   indexReorderActual = [];
%   d.actualIsNotNaN = [];
% end
