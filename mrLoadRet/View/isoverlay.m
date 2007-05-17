function val = isscan(overlay)
% function val = isview(overlay)
% 
% djh, 2007

val = 1;

if ieNotDefined('overlay')
    val = 0;
    return
end

if ~isstruct(overlay)
	val = 0;
	return
end

requiredFields = {'alpha','clip','colormap','data','date','function',...
	'groupName','interrogator','name','params','range','reconcileFunction'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(overlay,fieldName)
		% mrWarnDlg(['Invalid overlay, missing field: ',fieldName]);
		val = 0;
	end
end