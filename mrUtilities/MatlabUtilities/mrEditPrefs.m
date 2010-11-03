function prefParams = mrEditPrefs()
%
% Dialog box for changing mrDEFAULTS.prefs
%
% jg & djh, 5/2007

% get pref names /defaults
[prefNames prefDefaults] = mrGetPref();

% cycle through preference name list/defaults
% and make variables that have their value
for i = 1:length(prefNames)
  % set a variable with pref to its value
  eval(sprintf('%s = mrGetPref(prefNames{i});',prefNames{i}));
  % if there is a defaults list, then
  % make the variable have a cell array with
  % the current setting + all the choices
  if ~isempty(prefDefaults{i}) && iscell(prefDefaults{i}) && ~iscell(prefDefaults{i}{1})
    % check to see if the preference is in the list
    eval(sprintf('%s = putOnTopOfList(mrGetPref(prefNames{i}),prefDefaults{i});',prefNames{i}));
  end
end

% get defaultInterrogators
systemInterrogators = mrGetPref('systemInterrogators');
defaultInterrogators = mrGetPref('defaultInterrogators');
if isempty(defaultInterrogators)
  defaultInterrogators = '';
end
% only show ones that are not systemInterrgators and
% make into a comma delimited list so that the user
% can edit it easily
if isstr(defaultInterrogators)
  defaultInterrogators = commaDelimitedToCell(defaultInterrogators);
end
defaultInterrogators = cellToCommaDelimited(setdiff(defaultInterrogators,systemInterrogators));

% set up the dialog and ask the user to set parameters
prefParams = {{'site',site,'Where you are using this code'},...
    {'verbose',verbose,'Yes if you want to have dialog waitbars, No to have information printed to the terminal'},...
    {'interpMethod',interpMethod,'Type of interpolation to use. Normally this is set to nearest for nearest neighbor interpolation'},...
    {'maxBlocksize',maxBlocksize,'Size of chunks of data to analyze at a time. If you are running out of memory, set lower. A good starting value is 250000000','minmax=[0 inf]','incdec=[-10000000 10000000]'},...
    {'volumeDirectory',volumeDirectory,'The directory to default to when you load base anatomy from the Volume directory'},...
    {'overwritePolicy',overwritePolicy,'Method to use when analysis is going to overwrite an existing file'},...
    {'niftiFileExtension',niftiFileExtension,'Nifti file extension, usually .img'},...
    {'selectedROIColor',selectedROIColor,'What color to use to draw the selected ROI.'},...
    {'roiContourWidth',roiContourWidth,'What line width to use to draw the selected ROI.'},...
    {'roiPolygonMethod',roiPolygonMethod,'Method used to create ROI polygons. The default roipoly function calls the line drawing function which can be very slow if you have already drawn a buch of lines (i.e. have some ROIs displaying). If you choose getpts instead, you will not have the lines drawn between points as you draw the ROI, but it will be much faster. getptsNoDoubleClick is similar to getpts but instead of double-click to end the selection you hit the return key (on some machines Matlab''s idea of what constitutes a double-click can be very slow)'},...
    {'defaultInterrogators',defaultInterrogators,'This is a comma separated list that contains interrogators that you can use.'},...
    {'roiCacheSize',roiCacheSize,'Size of ROI cache, usually 100.','minmax=[0 inf]','incdec=[-1 1]'},...
    {'baseCacheSize',baseCacheSize,'Size of base image cache. Set to the number of base slices you want to be able to quickly view','minmax=[0 inf]','incdec=[-1 1]'},...
    {'overlayCacheSize',overlayCacheSize,'Size of overlay image cache. Set to the number of base slices you want to be able to quickly view','minmax=[0 inf]','incdec=[-1 1]'}};

% open up dialog
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
  mrSetPref('roiContourWidth',prefParams.roiContourWidth);
  mrSetPref('defaultInterrogators',commaDelimitedToCell(prefParams.defaultInterrogators));
  mrSetPref('roiCacheSize',prefParams.roiCacheSize);
  mrSetPref('baseCacheSize',prefParams.baseCacheSize);
  mrSetPref('overlayCacheSize',prefParams.overlayCacheSize);
end


