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
% $Id$	
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  requiredFields = {'groupName','name','params','range'};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  requiredFields = {'function','groupName','name','params','range',...
                    'date','type','alpha','clip','colormap',...
                    'interrogator','reconcileFunction','mergeFunction',...
		    'data','alphaOverlay','alphaOverlayExponent'};
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

%if not valid, no point in trying to make it valid
if ~tf
  return
end

%do not allow empty range
if isempty(overlay.range) || any(size(overlay.range)~=[1 2])
  overlay.range(1) = -1;  
  overlay.range(2) = 1; 
end

%only set optional fields once required fields have been checked
if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  optionalFields = {'function','unknown';
                    'date',datestr(now);
                    'type',overlay.name;
                    'alpha',1;
                    'clip',overlay.range;
                    'colorRange',overlay.range;
                    'colormap',jet(256);
                    'interrogator','timecoursePlot';
                    'reconcileFunction','defaultReconcileParams';
                    'mergeFunction','defaultMergeParams';
                    'colormapType', 'normal';
                    'data',[];
		    'alphaOverlay',[];
		    'alphaOverlayExponent',1;
		   };
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the analysis structure it is invalid).
  optionalFields = {'colorRange',overlay.range};
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

%Do not allow empty clip and color ranges and force range to include clip range (for 'New' overlay range behaviour only)
if tf
  if isempty(overlay.clip) || any(size(overlay.clip)~=[1 2])
    overlay.clip = overlay.range;
  end
  if isempty(overlay.colorRange) || any(size(overlay.colorRange)~=[1 2])
    overlay.colorRange = overlay.range;
  end
  if strcmp(mrGetPref('overlayRangeBehaviour'),'New')
    overlay.range(1) = min(overlay.clip(1),overlay.range(1));
    overlay.range(2) = max(overlay.clip(2),overlay.range(2));
  end
end
  
