function [tf overlay] =  isoverlay(overlay)
% function [tf overlay] =  isoverlay(overlay)
%
% Checks to see if it is a valid overlay structure. Can be called with
% either one or two output arguments:
%
% tf =  isoverlay(overlay)
% [tf overlay] =  isoverlay(overlay)
%
% tf is logical 1 (true) if overlay is a valid overlay structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid overlay structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  requiredFields = {'function','groupName','name','params','range'};
  optionalFields = {'date',datestr(now);
    'type',overlay.name;
    'alpha',1;
    'clip',overlay.range;
    'colormap',jet(256);
    'interrogator','mrDefaultInterrogator';
    'reconcileFunction','defaultReconcileParams';
    'mergeFunction','defaultMergeParams';
    'data',[]};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'function','groupName','name','params','range',...
    'date','type','alpha','clip','colormap',...
    'interrogator','reconcileFunction','mergeFunction','data'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('overlay')
    tf = false;
    return
end
if ~isstruct(overlay)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(overlay,fieldName)
		% mrWarnDlg(['Invalid overlay, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(overlay,fieldName)  
    % use eval args to set the fields properly
    varargin{1} = sprintf('overlay.%s',fieldName);
    varargin{2} = default;
    eval(evalargs(varargin),1);
  end
end
overlay = orderfields(overlay);
