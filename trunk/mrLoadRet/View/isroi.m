function val = isbase(roi)
% function val = isview(roi)
% 
% djh, 2007

val = 1;

if ieNotDefined('roi')
    val = 0;
    return
end

if ~isstruct(roi)
	val = 0;
	return
end

requiredFields = {'color','coords','date','name','viewType','voxelSize','xform'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(roi,fieldName)
		% mrWarnDlg(['Invalid ROI, missing field: ',fieldName]);
		val = 0;
	end
end