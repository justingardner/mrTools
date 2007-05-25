function val = isview(view)
% function val = isview(view)
% 
% djh, 2004

mrGlobals

if ieNotDefined('view')
    val = 0;
    return
end

if ~isstruct(view)
	val = 0;
	return
end

% Check that it has the required fields
requiredFields = {'viewNum','viewType','baseVolumes','curBase','curGroup',...
    'analyses','curAnalysis','ROIs','curROI'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(view,fieldName)
		% mrWarnDlg(['Invalid view, missing field: ',fieldName]);
		val = 0;
		return
	end
end

% Check that viewNum is a number
if ~isnumeric(view.viewNum);
    val = 0;
    return
end

% Confirm that MLR.views{viewNum} and view have the same fields
names1 = fieldnames(MLR.views{view.viewNum});
names2 = fieldnames(view);
if length(names1) == length(names2)
  val = all(strcmp(names1,names2));
else
  val = 0;
end
