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
if nargin == 0
  % set the defaults for preference we have defaults for
  prefDefaults{length(prefNames)} = [];
  prefDefaults{find(strcmp('interpMethod',prefNames))} = {'nearest','linear','spline','cubic'};
  prefDefaults{find(strcmp('overwritePolicy',prefNames))} = {'Ask','Merge','Rename','Overwrite'};
  prefDefaults{find(strcmp('verbose',prefNames))} = {'Yes','No'};
  prefDefaults{find(strcmp('niftiFileExtension',prefNames))} = {'.img','.nii'};
  prefDefaults{find(strcmp('roiPolygonMethod',prefNames))} = {'roipoly','getpts','getptsNoDoubleClick'};
  prefDefaults{find(strcmp('selectedROIColor',prefNames))} = color2RGB;
  prefDefaults{find(strcmp('selectedROIColor',prefNames))}{end+1} = 'none';
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
    value = [];
end

% default value for selectedROIColor
if strcmp(pref,'selectedROIColor') && isempty(value)
  value = 'white';
end

% deal with interrogators
if strcmp(pref,'systemInterrogators')
  value = {'timecoursePlot','makeFlat','searchForVoxel'};
end
if strcmp(pref,'defaultInterrogators')
  if isempty(value)
    value = mrGetPref('systemInterrogators');
  else
    value = union(value,mrGetPref('systemInterrogators'));
  end
end
