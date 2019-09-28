% [fdrAdjustedP,fweAdjustedP] = multipleTestsAdjustment(p, <fdrAdjust, fweAdjust>)
%
%   adjusts p values using False Discovery Rate Step-up method and Hommel Bonferroni correction
%   optional inputs: set fdrAdjust or fweAdjust to 0 to skip either method
%     
% jb 15/03/2012
%
% $Id: maskAwithB.m 2172 2011-06-20 12:49:44Z julien $ 
function [fdrAdjustedP,fweAdjustedP] = multipleTestsAdjustment(p, fdrAdjust, fweAdjust)


if ~ismember(nargin,[1 2 3])
   help multipleTestsAdjustment;
   return
end
if ieNotDefined('fdrAdjust')
  fdrAdjust=1;
end
if ieNotDefined('fweAdjust')
  fweAdjust=1;
end
 params.fdrAdjustment= fdrAdjust;
 params.fweAdjustment= fweAdjust;
 
[~, fdrAdjustedP, fweAdjustedP] = transformStatistic(p,[],params);