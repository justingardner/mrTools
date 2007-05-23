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

requiredFields = {'function','groupName','name','params','range'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(overlay,fieldName)
		% mrWarnDlg(['Invalid overlay, missing field: ',fieldName]);
		val = 0;
	end
end

% Optional fields and defaults
optionalFields = {'date',dateString;
  'type',overlay.name;
  'alpha',1;
  'clip',overlay.range;
  'colormap',jet(256);
  'interrogator','mrDefaultInterrogator';
  'reconcileFunction','defaultReconcileFunction';
  'mergeFunction','defaultMergeFunction';
  'data',[]};
for f = 1:length(optionalFields)
	fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(analysis,fieldName)
    analysis.fieldName = default;
  end
end

if val ==1
  val = orderfields(overlay);
end