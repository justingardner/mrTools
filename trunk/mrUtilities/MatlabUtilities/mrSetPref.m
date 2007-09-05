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
% djh, 5/2007

global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

mrDEFAULTS.prefs = setfield(mrDEFAULTS.prefs,pref,value);
preferences = mrDEFAULTS.prefs;
saveMrDefaults;