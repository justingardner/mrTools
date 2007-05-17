function val = isbase(baseAnatomy)
% function val = isview(baseAnatomy)
% 
% djh, 2007

val = 1;

if ieNotDefined('baseAnatomy')
    val = 0;
    return
end

if ~isstruct(baseAnatomy)
	val = 0;
	return
end

requiredFields = {'clip','data','hdr','name','permutationMatrix','range'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(baseAnatomy,fieldName)
		% mrWarnDlg(['Invalid baseAnatomy, missing field: ',fieldName]);
		val = 0;
	end
end