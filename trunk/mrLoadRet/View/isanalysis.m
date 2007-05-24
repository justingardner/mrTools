function val = isanalysis(analysis)
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

% Required fields
requiredFields = {'function','groupName','guiFunction','name','params'};
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(analysis,fieldName)
		% mrWarnDlg(['Invalid analysis, missing field: ',fieldName]);
		val = 0;
	end
end

% Optional fields and defaults
optionalFields = {'date',datestr(now);
  'type',analysis.name;
  'overlays',[];
  'curOverlay',[];
  'reconcileFunction','defaultReconcileParams';
  'mergeFunction','defaultMergeParams';
  'd',[]};
for f = 1:length(optionalFields)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(analysis,fieldName)
    % use eval args to set the fields properly
    varargin{1} = sprintf('analysis.%s',fieldName);
    varargin{2} = default;
    eval(evalargs(varargin));
  end
end

if val ==1
  val = orderfields(analysis);
end