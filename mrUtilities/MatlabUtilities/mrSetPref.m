function preferences = mrSetPref(pref,value)
%
% preferences = mrSetPref(pref,value)
%
% Replaces Matlab's setpref function. Sets a field in the global variable
% mrDEFAULTS.preferences, which is a structure with fields for each
% preference. Returns the preferences structure.
%
% Examples:
%   mrSetPref('verbose','Yes');
%   mrSetPref('verbose','No');
%   mrSetPref('site','NYU');
%   mrSetPref('niftiFileExtension','.img');
%   mrSetPref('niftiFileExtension','.nii');
%   mrSetPref('interpMethod','nearest');
%      Options: 'nearest','linear','spline','cubic'
%   mrSetPref('overwritePolicy','Ask');
%      Options: 'Ask','Merge','Rename','Overwrite'
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
% djh, 5/2007

global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

% get pref names and defaults
[prefNames prefDefaults] = mrGetPref;

% check for a known preference
prefNum = find(strcmp(lower(pref),lower(prefNames)));
if isempty(prefNum)
  % print message for unknown preference
  disp(sprintf('(mrSetPref) Unknown preference %s',pref));
else
  % this will fix the caps on the prefs name
  pref = prefNames{prefNum};
  % check for a known default
  if ~isempty(prefDefaults{prefNum})
    prefDefaultNum = find(strcmp(lower(value),lower(prefDefaults{prefNum})));
    % print message if it is not known
    if isempty(prefDefaultNum)
      disp(sprintf('(mrSetPref) Value %s for preference %s is not known',value,pref));
    else
      % fix caps
      value = prefDefaults{prefNum}{prefDefaultNum};
    end
  end
end

mrDEFAULTS.prefs = setfield(mrDEFAULTS.prefs,pref,value);
preferences = mrDEFAULTS.prefs;
saveMrDefaults;