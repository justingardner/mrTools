function val = isscan(analysis)
% function val = isview(analysis)
% 
% djh, 2007

val = 1;

if ieNotDefined('analysis')
    val = 0;
    return
end

if ~isstruct(analysis)
	val = 0;
	return
end

requiredFields = {'curOverlay','date','function','groupName','guiFunction','name',...
	'overlays','params','reconcileFunction','type'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(analysis,fieldName)
		mrWarnDlg(['Invalid analysis, missing field: ',fieldName]);
		val = 0;
	end
end