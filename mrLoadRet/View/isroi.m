function [tf roi] =  isroi(roi)
% function [tf roi] =  isroi(roi)
%
% Checks to see if it is a valid roi structure. Can be called with
% either one or two output arguments:
%
% tf =  isroi(roi)
% [tf roi] =  isroi(roi)
%
% tf is logical 1 (true) if roi is a valid roi structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid roi structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the roi with optional
  % fields is valid.
  requiredFields = {'name','viewType','voxelSize','xform'};
  optionalFields = {'coords',[];
    'date',datestr(now);
    'color','blue';
    'notes',''};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'color','coords','date','name','viewType','voxelSize','xform'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('roi')
    tf = false;
    return
end
if ~isstruct(roi)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(roi,fieldName)
		% mrWarnDlg(['Invalid roi, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(roi,fieldName)  
    % use eval args to set the fields properly
    varargin{1} = sprintf('roi.%s',fieldName);
    varargin{2} = default;
    eval(evalargs(varargin),1);
  end
end
roi = orderfields(roi);
