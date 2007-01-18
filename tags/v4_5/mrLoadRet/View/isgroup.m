function val = isscan(group)
% function val = isview(group)
% 
% djh, 2007

val = 1;

if ieNotDefined('group')
    val = 0;
    return
end

if ~isstruct(group)
	val = 0;
	return
end

requiredFields = {'name','scanParams'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(group,fieldName)
		mrWarnDlg(['Invalid group, missing field: ',fieldName]);
		val = 0;
	end
end