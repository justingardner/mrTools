function prefParams = mrEditPrefs()
%
% Dialog box for changing mrDEFAULTS.prefs
%
% jg & djh, 5/2007

% get interpTypes
interpTypes = {'nearest','linear','spline','cubic'};
if find(strcmp(mrGetPref('interpMethod'),interpTypes))
    interpTypes = putOnTopOfList(mrGetPref('interpMethod'),interpTypes);
end

% get overwritePolicy
overwritePolicy = {'Ask','Merge','Rename','Overwrite'};
if find(strcmp(mrGetPref('overwritePolicy'),overwritePolicy))
    overwritePolicy = putOnTopOfList(mrGetPref('overwritePolicy'),overwritePolicy);
end

% get verbose
verbose = {'Yes','No'};
if find(strcmp(mrGetPref('verbose'),verbose))
    verbose = putOnTopOfList(mrGetPref('verbose'),verbose);
end

% get niftiFileExtension
niftiFileExtension = {'.img','.nii'};
if find(strcmp(mrGetPref('niftiFileExtension'),niftiFileExtension))
    niftiFileExtension = putOnTopOfList(mrGetPref('niftiFileExtension'),niftiFileExtension);
end

% get current values for other prefs
site = mrGetPref('site');
maxBlocksize = mrGetPref('maxBlocksize');
volumeDirectory = mrGetPref('volumeDirectory');

% get values for cache sizes
roiCacheSize = mrGetPref('roiCacheSize');
baseCacheSize = mrGetPref('baseCacheSize');
overlayCacheSize = mrGetPref('overlayCacheSize');

% set up the dialog and ask the user to set parameters
prefParams = {{'site',site,'Where you are using this code'},...
    {'verbose',verbose,'Yes if you want to have dialog waitbars, No to have information printed to the terminal'},...
    {'interpMethod',interpTypes,'Type of interpolation to use. Normally this is set to nearest for nearest neighbor interpolation'},...
    {'maxBlocksize',maxBlocksize,'Size of chunks of data to analyze at a time. If you are running out of memory, set lower. A good starting value is 250000000','minmax=[0 inf]','incdec=[-10000000 10000000]'},...
    {'volumeDirectory',volumeDirectory,'The directory to default to when you load base anatomy from the Volume directory'},...
    {'overwritePolicy',overwritePolicy,'Method to use when analysis is going to overwrite an existing file'},...
    {'niftiFileExtension',niftiFileExtension,'Nifti file extension, usually .img'},...
    {'roiCacheSize',roiCacheSize,'Size of ROI cache, usually 3.','minmax=[0 inf]','incdec=[-1 1]'},...
    {'baseCacheSize',baseCacheSize,'Size of base image cache. Set to the number of base slices you want to be able to quickly view','minmax=[0 inf]','incdec=[-1 1]'},...
    {'overlayCacheSize',overlayCacheSize,'Size of overlay image cache. Set to the number of base slices you want to be able to quickly view','minmax=[0 inf]','incdec=[-1 1]'}};
prefParams = mrParamsDialog(prefParams,'Set mrTools preferences');

% if they did not cancel then actually set the parameters
if ~isempty(prefParams)
    mrSetPref('site',prefParams.site);
    mrSetPref('verbose',prefParams.verbose);
    mrSetPref('interpMethod',prefParams.interpMethod);
    mrSetPref('maxBlocksize',prefParams.maxBlocksize);
    mrSetPref('volumeDirectory',prefParams.volumeDirectory);
    mrSetPref('overwritePolicy',prefParams.overwritePolicy);
    mrSetPref('niftiFileExtension',prefParams.niftiFileExtension);
    mrSetPref('roiCacheSize',prefParams.roiCacheSize);
    mrSetPref('baseCacheSize',prefParams.baseCacheSize);
    mrSetPref('overlayCacheSize',prefParams.overlayCacheSize);
end
