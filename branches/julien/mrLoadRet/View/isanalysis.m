function [tf analysis] =  isanalysis(analysis)
% function [tf analysis] =  isanalysis(analysis)
%
% Checks to see if it is a valid analysis structure. Can be called with
% either one or two output arguments:
%
% tf =  isanalysis(analysis)
% [tf analysis] =  isanalysis(analysis)
%
% tf is logical 1 (true) if analysis is a valid analysis structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid analysis structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  requiredFields = {'function','groupName','guiFunction','name','params'};
  optionalFields = {'date',datestr(now);
    'type',analysis.name;
    'overlays',[];
    'curOverlay',[];
    'reconcileFunction','defaultReconcileParams';
    'mergeFunction','defaultMergeParams';
    'd',[];
    'clipAcrossOverlays',1};
else
  % Return 0 if the analysis structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'function','groupName','guiFunction','name','params',...
    'date','type','overlays','curOverlay','reconcileFunction','mergeFunction','d','clipAcrossOverlays'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('analysis')
    tf = false;
    return
end
if ~isstruct(analysis)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(analysis,fieldName)
		% mrWarnDlg(['Invalid analysis, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(analysis,fieldName)  
    % use eval args to set the fields properly
    varargin{1} = sprintf('analysis.%s',fieldName);
    varargin{2} = default;
    eval(evalargs(varargin),1);
  end
end
analysis = orderfields(analysis);
