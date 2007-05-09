function value = mrGetPref(pref,value)
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
% djh, 5/2007

global mrDEFAULTS

if isfield(mrDEFAULTS.prefs,pref)
    value = getfield(mrDEFAULTS.prefs,pref);
else
    value = [];
end
