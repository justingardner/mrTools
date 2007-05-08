function preferences = mrSetPref(pref,value)
%
% preferences = mrSetPref(pref,value)
%
% Replaces Matlab's setpref function. Sets a field in the global variable
% mrDEFAULTS.preferences, which is a structure with fields for each
% preference. Returns the preferences structure.
%
% Examples:
%   mrSetPref('verbose',1);
%   mrSetPref('verbose',0);
%   mrSetPref('site','NYU');
%   mrSetPref('niftiFileExtension','.img');
%   mrSetPref('niftiFileExtension','.nii');
%   mrSetPref('interpMethod','nearest');
%   mrSetPref('interpMethod','linear');
%   mrSetPref('interpMethod','cubic');
%   mrSetPref('overwritePolicy','ask');
%
% djh, 5/2007

