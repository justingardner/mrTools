function [tf base] =  isbase(base)
% function [tf base] =  isbase(base)
%
% Checks to see if it is a valid base structure. Can be called with
% either one or two output arguments:
%
% tf =  isbase(base)
% [tf base] =  isbase(base)
%
% tf is logical 1 (true) if base is a valid base structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid base structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  requiredFields = {'data','hdr','name','permutationMatrix'};
  optionalFields = {'range',[min(base.data(:)) max(base.data(:))];
		    'clip',defaultClip(base.data);
		    'coordMap',[];
		    'rotate',[];
		    'curSlice',[];
		    'sliceOrientation',[]};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the base structure it is invalid).
  requiredFields = {'clip','coordMap','curSlice','data','hdr','name','permutationMatrix','range','rotate','sliceOrientation'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ieNotDefined('base')
    tf = false;
    return
end
if ~isstruct(base)
	tf = false;
	return
end

% Check required fields
for f = 1:length(requiredFields)
	fieldName = requiredFields{f};
	if ~isfield(base,fieldName)
		% mrWarnDlg(['Invalid base, missing field: ',fieldName]);
		tf = false;
	end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(base,fieldName)  
    base.(fieldName) = default;
  end
end
base = orderfields(base);



function clip = defaultClip(image)
% Choose default clipping based on histogram
histThresh = length(image(:))/1000;
[cnt, val] = hist(image(:),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
clip = [clipMin,clipMax];
