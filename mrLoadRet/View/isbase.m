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
		    'rotate',0;
		    'tilt',0;
		    'curSlice',[];
		    'sliceOrientation',[],;
		    'type',[],;
		    'gamma',1,;
		    'vol2tal',[];
		    'vol2mag',[];
                    'talInfo',[];
		    'originalOrient',[];
		    'xformFromOriginal',[];
        'curCorticalDepth',[]};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the base structure it is invalid).
  requiredFields = {'clip','coordMap','curSlice','data','hdr','name','permutationMatrix','range','rotate','sliceOrientation','type','gamma','tilt','vol2tal','vol2mag','talInfo'};
  optionalFields = {'curCorticalDepth'};
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

% hack to change talPoints field to talinfo
if isfield(base,'talPoints')
  disp(sprintf('(isbase) Changing talPoints field to talInfo'));
  if ~isfield(base,'talInfo') || isempty(base.talInfo)
    base.talInfo = base.talPoints;
  end
  % remove the legacy field name
  base = rmfield(base,'talPoints');
end

% Check required fields
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(base,fieldName)
    % mrWarnDlg(['Invalid base, missing field: ',fieldName]);
    tf = false;
    return
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

% remove any fields that are not required or optional
if nargout == 2
  baseFieldNames = fieldnames(base);
  for f = 1:length(baseFieldNames)
    % check required fields
    if ~any(strcmp(baseFieldNames{f},requiredFields))
      % check optional fields, (only check first column of
      % field names, not the default values...
      match = strcmp(baseFieldNames{f},optionalFields);
      if ~any(match(:,1))
	disp(sprintf('(isbase) Removing unnecessary field %s from base',baseFieldNames{f}));
	base = rmfield(base,baseFieldNames{f});
      end
    end
  end
end

% order the fields
base = orderfields(base);

% auto set the base type if it is set to empty
if isempty(base.type)
  base.type = ~isempty(base.coordMap);
end

% validate coordMap field
if ~isempty(base.coordMap) && any(base.type==[1 2]);
  [tf base.coordMap] = validateCoordMap(base.coordMap,base.type);
  if ~tf,disp(sprintf('(isbase) Invalid baseCoordMap'));end
end

% make sure that name does not include full path
base.name = getLastDir(base.name);

%%%%%%%%%%%%%%%%%%%%%
%%   defaultClip   %%
%%%%%%%%%%%%%%%%%%%%%
function clip = defaultClip(image)
% Choose default clipping based on histogram
histThresh = length(image(:))/1000;
[cnt, val] = hist(image(:),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
clip = [clipMin,clipMax];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   validateCoordMap   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf coordMap] = validateCoordMap(coordMap,baseType)

% First, we fix old coordMaps to have correctly named fields
if isfield(coordMap,'flatDir')
  disp(sprintf('(isbase) Updating format of coordMap'));
  newCoordMap.path = coordMap.flatDir;
  newCoordMap.dims = coordMap.dims;
  newCoordMap.flatFileName = coordMap.flatFileName;
  newCoordMap.innerCoordsFileName = coordMap.innerFileName;
  newCoordMap.outerCoordsFileName = coordMap.outerFileName;
  newCoordMap.curvFileName = coordMap.curvFileName;
  newCoordMap.anatFileName = coordMap.anatFileName;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  coordMap = newCoordMap;
elseif isfield(coordMap,'innerFileName')
  disp(sprintf('(isbase) Updating format of coordMap'));
  if isfield(coordMap,'path')
    newCoordMap.path = coordMap.path;
  else
    newCoordMap.path = '';
  end
  newCoordMap.dims = coordMap.dims;
  newCoordMap.innerSurfaceFileName = coordMap.innerFileName;
  newCoordMap.innerCoordsFileName = coordMap.innerCoordsFileName;
  newCoordMap.outerSurfaceFileName = coordMap.outerFileName;
  newCoordMap.outerCoordsFileName = coordMap.outerCoordsFileName;
  newCoordMap.curvFileName = coordMap.curvFileName;
  newCoordMap.anatFileName = coordMap.anatFileName;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  newCoordMap.innerVtcs = coordMap.innerVtcs;
  newCoordMap.outerVtcs = coordMap.outerVtcs;
  newCoordMap.tris = coordMap.tris;
  if isfield(coordMap,'rotate')
    newCoordMap.rotate = coordMap.rotate;
  end
  coordMap = newCoordMap;
elseif isfield(coordMap,'inner')
  disp(sprintf('(isbase) Updating format of coordMap'));
  if isfield(coordMap,'path')
    newCoordMap.path = coordMap.path;
  else
    newCoordMap.path = '';
  end
  newCoordMap.innerSurfaceFileName = coordMap.inner;
  newCoordMap.innerCoordsFileName = coordMap.inner;
  newCoordMap.outerSurfaceFileName = coordMap.outer;
  newCoordMap.outerCoordsFileName = coordMap.outer;
  newCoordMap.curvFileName = coordMap.curv;
  newCoordMap.anatFileName = coordMap.anatomy;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  newCoordMap.innerVtcs = coordMap.innerVtcs;
  newCoordMap.outerVtcs = coordMap.outerVtcs;
  newCoordMap.tris = coordMap.tris;
  if isfield(coordMap,'rotate')
    newCoordMap.rotate = coordMap.rotate;
  end
  coordMap = newCoordMap;
end

if baseType == 1
  requiredFields = {'path','dims','flatFileName','innerCoordsFileName','outerCoordsFileName','curvFileName','anatFileName','innerCoords','outerCoords'};
  optionalFields = {'coords',[]};
elseif baseType == 2
  requiredFields = {'path','dims','innerSurfaceFileName','innerCoordsFileName','outerSurfaceFileName','outerCoordsFileName','curvFileName','anatFileName','innerCoords','outerCoords','innerVtcs','outerVtcs','tris'};
  optionalFields = {'coords',[];'rotate',0};
end  

% Check required fields
tf = true;
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(coordMap,fieldName)
    disp(sprintf('(isbase) Missing coordMap field: %s',fieldName));
    tf = false;
  end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(coordMap,fieldName)  
    coordMap.(fieldName) = default;
  end
end

% make sure the filenames are set correctly
if baseType == 2
  if strcmp(coordMap.innerCoordsFileName,'Same as surface')
    coordMap.innerCoordsFileName = coordMap.innerSurfaceFileName;
  end
  if strcmp(coordMap.outerCoordsFileName,'Same as surface')
    coordMap.outerCoordsFileName = coordMap.outerSurfaceFileName;
  end
end


  
