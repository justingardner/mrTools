% enforceMonotonictyResampleFWE.m
%
%        $Id$
%      usage: count = enforceMonotonictyResampleFWE(p,d.indexSortActual,d.indexReorderActual,d.actualIsNotNaN,params)
%         by: julien besle, 
%       date: 14/01/2011
%    purpose: enforces ordered p-value monotonicity in Step-Down bootstrap-based FWE control
%             see algorithm 4.1, step 6, Westfall and Young (1993) p117

function p = enforceMonotonicityResampleFWE(p,d,params)

if strcmp(params.resampleFWEadjustment,'Step-down')
  %sort p using provided sort index (sorting of the non-nan actual statistics)
  %note that these values are sorted from the smallest to the largest

  %     reorderedActual = NaN(size(p));                                            % DEBUG
  for i=1:size(p,2)
    sortedP = p(d.indexSortActual{i},i);
  %       sortedActual = actual(d.indexSortActual{i});                                % DEBUG

    %compute consecutive maxima (from largest to smallest p value)
    for j = size(sortedP,1)-1:-1:1
      sortedP(j) = max(sortedP(j),sortedP(j+1));
    end
    p(d.actualIsNotNaN(:,i),i) = sortedP(d.indexReorderActual{i});
%     reorderedActual(d.actualIsNotNaN(:,i),i) = sortedActual(d.indexReorderActual{i});  % DEBUG
  end
  %p(~d.actualIsNotNaN,:)=NaN;
end


% %old version
% if strcmp(params.resampleFWEadjustment,'Step-down')
%   %sort p using provided sort index (sorting of the non-nan actual statistics)
%   %note that these values are sorted from the smallest to the largest
%   sortedP = NaN(size(d.indexSortActual));
%   for i=1:size(p,2)
%     sortedP(:,i) = p(d.indexSortActual(:,i),i);
% %       sortedActual(:,i) = actual(d.indexSortActual(:,i),i);                           % DEBUG
%   end
%   %compute consecutive maxima (from largest to smallest p value)
%   for i = size(sortedP,1)-1:-1:1
%     sortedP(i,:) = max(sortedP(i,:),sortedP(i+1,:));
%   end
%   for i=1:size(p,2)
%     p(d.actualIsNotNaN,i) = sortedP(d.indexReorderActual(:,i),i);
% %       reorderedActual(d.actualIsNotNaN,i) = sortedActual(d.indexReorderActual(:,i),i);  % DEBUG
%   end
%   %p(~d.actualIsNotNaN,:)=NaN;
% end