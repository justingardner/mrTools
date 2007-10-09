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

% get roiPolygonMethod
roiPolygonMethod = {'roipoly','getpts'};
if find(strcmp(mrGetPref('roiPolygonMethod'),roiPolygonMethod))
  roiPolygonMethod = putOnTopOfList(mrGetPref('roiPolygonMethod'),roiPolygonMethod);
end

% get defaultInterrogators
systemInterrogators = {'timecoursePlot','makeFlat'};
defaultInterrogators = mrGetPref('defaultInterrogators');
if isempty(defaultInterrogators)
%  defaultInterrogators = 'mrDefaultInterrogator';
  defaultInterrogators = '';
end
% only show ones that are not systemInterrgators and
% make into a comma delimited list so that the user
% can edit it easily
if isstr(defaultInterrogators)
  defaultInterrogators = commaDelimitedToCell(defaultInterrogators);
end
defaultInterrogators = cellToCommaDelimited(setdiff(defaultInterrogators,systemInterrogators));

% get selectedROIColor
selectedROIColor = mrGetPref('selectedROIColor');
if isempty(selectedROIColor)
  selectedROIColor = 'white';
end
selectedROIColor = putOnTopOfList(selectedROIColor,{'yellow','magenta','cyan','red','green','blue','white','black'});

% get roiMatchMethod
roiMatchMethod = mrGetPref('roiMatchMethod');
if isempty(roiMatchMethod)
  roiMatchMethod = 'Normal';
end
roiMatchMethod = putOnTopOfList(roiMatchMethod,{'Normal','No transform'});

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
    {'selectedROIColor',selectedROIColor,'What color to use to draw the selected ROI.'},...
    {'roiPolygonMethod',roiPolygonMethod,'Method used to create ROI polygons. The default roipoly function calls the line drawing function which can be very slow if you have already drawn a buch of lines (i.e. have some ROIs displaying). If you choose getpts instead, you will not have the lines drawn between points as you draw the ROI, but it will be much faster.'},...
    {'roiMatchMethod',roiMatchMethod,'If you are transforming an ROI to an image which has the same voxel dimensions and the same transform, you can use ''No transform'' to have the code just pass the voxels back instead of going through the normal transformation process. While this is fast it may leave holes in the ROI. The ''Normal'' transformation process will add some partial volumed voxels, which is the default. It is a bit slower, but does not leave wholes in the ROI'},...
    {'defaultInterrogators',defaultInterrogators,'This is a comma separated list that contains interrogators that you can use.'},...
    {'roiCacheSize',roiCacheSize,'Size of ROI cache, usually 50.','minmax=[0 inf]','incdec=[-1 1]'},...
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
    mrSetPref('roiPolygonMethod',prefParams.roiPolygonMethod);
    mrSetPref('selectedROIColor',prefParams.selectedROIColor);
    mrSetPref('defaultInterrogators',union(commaDelimitedToCell(prefParams.defaultInterrogators),systemInterrogators));
    mrSetPref('roiCacheSize',prefParams.roiCacheSize);
    mrSetPref('roiMatchMethod',prefParams.roiMatchMethod);
    mrSetPref('baseCacheSize',prefParams.baseCacheSize);
    mrSetPref('overlayCacheSize',prefParams.overlayCacheSize);
end


