% mlrImageIsHeader.m
%
%        $Id:$ 
%      usage: [tf h] mlrImageIsHeader(h)
%         by: justin gardner
%       date: 08/19/11
%    purpose: Test whether the passed in strucutre is an mlrImageHeader or not
%
function [tf h] = mlrImageIsHeader(h)

tf = false;
% check arguments
if ~any(nargin == [1])
  help mlrImageIsHeader
  return
end

if (nargout == 2)
  % Add optional fields and return true if the image header with optional fields is valid
  requiredFields = {'dim'};
  optionalFields = {'qform_code',0;
		    'qform44',[]'
		    'sform_code',0;
		    'sform44',[];
		    'vol2mag',[];
		    'vol2tal',[];
		    'base',[];
		    'hdr',[];
		    'ext','';
		    'type','unknown';
		    'filename','';
		   };
  % The fields nDim and pixdim will be gereated below if necessary
else
  % Return 0 if the image header structure is missing any fields required or
  % optional (since w/out changing the header structure it is invalid).
  requiredFields = {'dim','nDim','pixdim','qform_code','qform44','sform_code','sform44','vol2mag','vol2tal','base','hdr','ext','type','filename'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ~isstruct(h)
  tf = false;
  return
end

% Check required fields
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(h,fieldName)
    tf = false;
    return;
  end
end

% set nDim
if ~isfield(h,'nDim')
  h.nDim = length(h.dim);
end

% set pixdim to same length as dim with all ones
if ~isfield(h,'pixdim')
  h.pixdim = ones(size(h.dim));
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(h,fieldName)  
    h.(fieldName) = default;
  end
end

% make sure we have row arrays
h.dim = h.dim(:)';
h.pixdim = h.pixdim(:)';

% always have the same sort order of fields
h = orderfields(h);
