function mrEditPrefs()
%
% Dialog box for changing mrDEFAULTS.prefs
%
% jg & djh, 5/2007

% get interpTypes
interpTypes = {'nearest','linear','spline','cubic'};
if ~isempty(strcmp(mrGetPref('interpMethod'),interpTypes))
    interpTypes = putOnTopOfList(mrGetPref('interpMethod'),interpTypes);
end

% get overwritePolicy
overwritePolicy = {'Ask','Merge','Rename','Overwrite'};
if ~isempty(strcmp(mrGetPref('overwritePolicy'),overwritePolicy))
    overwritePolicy = putOnTopOfList(mrGetPref('overwritePolicy'),overwritePolicy);
end

% get current values for other prefs
site = mrGetPref('site');
verbose = mrGetPref('verbose');
maxBlocksize = mrGetPref('maxBlocksize');
volumeDirectory = mrGetPref('volumeDirectory');
niftiFileExtension = mrGetPref('niftiFileExtension');

% set up the dialog and ask the user to set parameters
prefParams = {{'site',site,'Where you are using this code'},...
    {'verbose',verbose,'minmax=[0 1]','incdec=[-1 1]','Set to 1 if you want to have dialog waitbars, set to 0 to have information printed to the terminal'},...
    {'interpMethod',interpTypes,'Type of interpolation to use. Normally this is set to nearest for nearest neighbor interpolation'},...
    {'maxBlocksize',maxBlocksize,'Size of chunks of data to analyze at a time. If you are running out of memory, set lower. A good starting value is 250000000','minmax=[0 inf]','incdec=[-10000000 10000000]'},...
    {'volumeDirectory',volumeDirectory,'The directory to default to when you load base anatomy from the Volume directory'},...
    {'overwritePolicy',overwritePolicy,'Method to use when analysis is going to overwrite an existing file'},...
    {'niftiFileExtension',niftiFileExtension,'Nifti file extension, usually .img'}};
prefParams = mrParamsDialog(prefParams);

% if they did not cancel then actually set the parameters
if ~isempty(prefParams)
    mrSetPref('site',prefParams.site);
    mrSetPref('verbose',prefParams.verbose);
    mrSetPref('interpMethod',prefParams.interpMethod);
    mrSetPref('maxBlocksize',prefParams.maxBlocksize);
    mrSetPref('volumeDirectory',prefParams.volumeDirectory);
    mrSetPref('overwritePolicy',prefParams.overwritePolicy);
    mrSetPref('niftiFileExtension',prefParams.niftiFileExtension);
end
