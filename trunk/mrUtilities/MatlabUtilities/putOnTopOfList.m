% putOnTopOfList.m
%
%      usage: putOnTopOfList(topVal,list)
%         by: justin gardner
%       date: 04/03/07
%    purpose: puts the topVal on top of the list
%
% e.g.
%putOnTopOfList('topItem',{'oneItem','twoItem','topItem','threeItem'});
%
function outList = putOnTopOfList(topVal,inList)

% check arguments
if ~any(nargin == [2])
  help putOnTopOfList
  return
end

outList{1} = topVal;
inList = setdiff(inList,topVal);
for i = 1:length(inList)
  outList{end+1} = inList{i};
end


