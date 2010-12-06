% mrParamsCopyFields.m
%
%      usage: params1 = mrParamsCopyFields(params1,params2)
%         by: julien besle 
%       date: 04/12/2010
%    purpose: copies (and replaces) all the fields from one params structure to another
%              $Id$

function params2 = mrParamsCopyFields(params1,params2)

fieldNames=fieldnames(params1);
for iField = 1:length(fieldNames)
  for iParams = 1:length(params1.paramInfo)
    if strcmp(params1.paramInfo{iParams}{1},fieldNames{iField})
      thisParamInfo = params1.paramInfo{iParams};
      break;
    end
  end
  if isfield(params2,fieldNames{iField}) && isfield(params2,'paramInfo')
    params2.(fieldNames{iField})=eval(['params1.' fieldNames{iField}]);
    for jParams = 1:length(params2.paramInfo)
      if strcmp(params2.paramInfo{jParams}{1},fieldNames{iField})
        params2.paramInfo{jParams}=thisParamInfo;
        break;
      end
    end
  else
    params2.(fieldNames{iField})=eval(['params1.' fieldNames{iField}]);
    if isfield(params2,'paramInfo')
      params2.paramInfo{end+1} = thisParamInfo;
    else
      params2.paramInfo{1} = thisParamInfo;
    end
  end
end
