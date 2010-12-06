% copyFields.m
%
%      usage: struct1 = copyFields(struct1,struct2)
%         by: julien besle 
%       date: 03/12/2010
%    purpose: copies (and replaces) all the fields from one structure to another
%              $Id$

function struct2 = copyFields(struct1,struct2)


if ieNotDefined('struct2')
  struct2=[];
end

fieldNames=fieldnames(struct1);
for iField = 1:length(fieldNames)
  struct2.(fieldNames{iField})=eval(['struct1.' fieldNames{iField}]);
end
