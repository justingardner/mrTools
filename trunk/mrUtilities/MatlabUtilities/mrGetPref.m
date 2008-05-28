function [value prefDefaults] = mrGetPref(pref)
%
% value = mrGetPref(pref)
%
% Replaces Matlab's getpref function. Gets a field from the global variable
% mrDEFAULTS.preferences, which is a structure with fields for each
% preference. Returns the value for that preference.
%
% Examples:
%   mrGetPref('verbose');
%   mrGetPref('site');
%   mrGetPref('niftiFileExtension');
%   mrGetPref('interpMethod');
%   mrGetPref('overwritePolicy');
%   mrGetPref('maxBlockSize');
%   mrGetPref('volumeDirectory');
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
%   mrGetPref with no arguments returns a list of known preference and known default values
%   [prefNames prefefaults] = mrGetPref;
%
% djh, 5/2007
% %	$Id$	

% with no arguments, return a list of possible preferences
prefNames = {'interpMethod','overwritePolicy','verbose','niftiFileExtension','roiPolygonMethod','systemInterrogatros','selectedROIColor','site','maxBlocksize','volumeDirectory','roiCacheSize','baseCacheSize','overlayCacheSize','defaultInterrogators','importROIPath'};

% set the defaults for preference we have defaults for
prefDefaults{length(prefNames)} = [];
prefDefaults{find(strcmp('interpMethod',prefNames))} = {'nearest','linear','spline','cubic'};
prefDefaults{find(strcmp('overwritePolicy',prefNames))} = {'Ask','Merge','Rename','Overwrite'};
prefDefaults{find(strcmp('verbose',prefNames))} = {'Yes','No'};
prefDefaults{find(strcmp('niftiFileExtension',prefNames))} = {'.img','.nii'};
prefDefaults{find(strcmp('roiPolygonMethod',prefNames))} = {'roipoly','getpts','getptsNoDoubleClick'};
prefDefaults{find(strcmp('selectedROIColor',prefNames))} = color2RGB;
prefDefaults{find(strcmp('selectedROIColor',prefNames))}{end+1} = 'none';
prefDefaults{find(strcmp('site',prefNames))} = 'NYU';
prefDefaults{find(strcmp('maxBlocksize',prefNames))} = 250000000;
prefDefaults{find(strcmp('volumeDirectory',prefNames))} = '';
prefDefaults{find(strcmp('roiCacheSize',prefNames))} = 100;
prefDefaults{find(strcmp('baseCacheSize',prefNames))} = 50;
prefDefaults{find(strcmp('overlayCacheSize',prefNames))} = 50;

if nargin == 0
  % return arguments
  value = prefNames;
  return
end

global mrDEFAULTS

% fix any caps differences
prefNum = find(strcmp(lower(pref),lower(prefNames)));
if ~isempty(prefNum)
  pref = prefNames{prefNum};
end

% read the preferences and figlocs
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

if isfield(mrDEFAULTS.prefs,pref)
    value = getfield(mrDEFAULTS.prefs,pref);
else
  % not set yet, take the top most possibility in the default
  % list, otherwise return empty
  if ~isempty(prefNum) && ~isempty(prefDefaults{prefNum})
    if iscell(prefDefaults{prefNum})
      value = prefDefaults{prefNum}{1};
    else
      value = prefDefaults{prefNum};
    end
  else
    value = [];
  end
end

% default value for selectedROIColor
if strcmp(pref,'selectedROIColor') && isempty(value)
  value = 'white';
end

% deal with interrogators
if strcmp(pref,'systemInterrogators')
  if isempty(value)
    value = {'timecoursePlot','makeFlat','searchForVoxel'};
  end
end
if strcmp(pref,'defaultInterrogators')
  if isempty(value)
    value = mrGetPref('systemInterrogators');
  else
    value = union(value,mrGetPref('systemInterrogators'));
  end
end
