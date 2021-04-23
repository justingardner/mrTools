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
%   [prefNames prefDefaults] = mrGetPref;
%
%   if you type mrGetPref alone, it will print out all known preferences
%   mrGetPref;
%
% djh, 5/2007
% %	$Id: mrGetPref.m 2890 2013-10-27 12:52:21Z justin $	

% with no arguments, return a list of possible preferences
prefNames = {'overwritePolicy','verbose','graphWindow','checkParamsConsistency'...
   'maxBlocksize','roiCacheSize','baseCacheSize','overlayCacheSize','defaultPrecision',...
   'defaultInterrogators','systemInterrogators','interrogatorPaths',...,
   'importROIPath','volumeDirectory','niftiFileExtension','fslPath',...
   'selectedROIColor','roiContourWidth','roiPolygonMethod',...,
   'interpMethod','corticalDepthBins','roiCorticalDepthDisplayRatio',...
   'multiSliceProjectionMethod','colorBlending','overlayRangeBehaviour','baseNaNsColor',...
   'pluginPaths','selectedPlugins',...
   'statisticalTestOutput',...
   'site','magnet','coil','pulseSequence',...
   'maxArrayWidthForParamsDialog','maxArrayHeightForParamsDialog',...
   'mlrVolDisplayControls','mlrVolOverlayAlpha','motionCompDefaultParams','colorNames',...
   'mlrPath','vistaPath','lastPath'...
   'overlayCombineTransformPaths','roiTransformPaths','colormapPaths'...
	    };

% set the defaults for preference we have defaults for. Note that the "find" in
% here is to make sure that the prefDefaults list matches the prefNames order
prefDefaults{length(prefNames)} = [];
prefDefaults{strcmp('overwritePolicy',prefNames)} = {'Ask','Merge','Rename','Overwrite'};
prefDefaults{strcmp('verbose',prefNames)} = {'No','Yes'};
prefDefaults{strcmp('graphWindow',prefNames)} = {'Replace','Make new'};
prefDefaults{strcmp('checkParamsConsistency',prefNames)} = {'Yes','No'};
prefDefaults{strcmp('maxBlocksize',prefNames)} = 250000000;
prefDefaults{strcmp('roiCacheSize',prefNames)} = 100;
prefDefaults{strcmp('baseCacheSize',prefNames)} = 50;
prefDefaults{strcmp('overlayCacheSize',prefNames)} = 50;
prefDefaults{strcmp('defaultPrecision',prefNames)} = 'double';
prefDefaults{strcmp('interrogatorPaths',prefNames)} = '';
prefDefaults{strcmp('volumeDirectory',prefNames)} = '';
prefDefaults{strcmp('niftiFileExtension',prefNames)} = {'.img','.nii'};
prefDefaults{strcmp('fslPath',prefNames)} = 'FSL not installed';
prefDefaults{strcmp('selectedROIColor',prefNames)} = color2RGB;
prefDefaults{strcmp('selectedROIColor',prefNames)}{end+1} = 'none';
prefDefaults{strcmp('roiContourWidth',prefNames)} = 1;
prefDefaults{strcmp('roiCorticalDepthDisplayRatio',prefNames)} = .5;
prefDefaults{strcmp('roiPolygonMethod',prefNames)} = {'getpts','roipoly','getptsNoDoubleClick'};
prefDefaults{strcmp('interpMethod',prefNames)} = {'nearest','linear','spline','cubic'};
prefDefaults{strcmp('corticalDepthBins',prefNames)} = 11;
prefDefaults{strcmp('multiSliceProjectionMethod',prefNames)} = {'Average','Maximum Intensity Projection'};
prefDefaults{strcmp('colorBlending',prefNames)} = {'Additive','Alpha blend','Contours'};
prefDefaults{strcmp('overlayRangeBehaviour',prefNames)} = {'Classic','New'};
prefDefaults{strcmp('baseNaNsColor',prefNames)} = {'Black','White','Transparent'};
prefDefaults{strcmp('pluginPaths',prefNames)} = '';
prefDefaults{strcmp('selectedPlugins',prefNames)} = '';
prefDefaults{strcmp('statisticalTestOutput',prefNames)} = {'P value','Z value','-log10(P) value'};
prefDefaults{strcmp('site',prefNames)} = 'NYU';
prefDefaults{strcmp('magnet',prefNames)} = {{'Allegra 3T','other'}};
prefDefaults{strcmp('coil',prefNames)} = {{'LifeService','Siemens birdcage','Nova birdcage','Nova surface','Nova quadrapus','Nova visual array','other'}};
prefDefaults{strcmp('pulseSequence',prefNames)} = {{'cbi_ep2d_bold','other'}};
prefDefaults{strcmp('maxArrayWidthForParamsDialog',prefNames)} = 25;
prefDefaults{strcmp('maxArrayHeightForParamsDialog',prefNames)} = 50;
prefDefaults{strcmp('mlrVolDisplayControls',prefNames)} = false;
prefDefaults{strcmp('mlrVolOverlayAlpha',prefNames)} = 0.8;
prefDefaults{strcmp('motionCompDefaultParams',prefNames)} = [];
prefDefaults{strcmp('colorNames',prefNames)} = {};
prefDefaults{strcmp('mlrPath',prefNames)} = '';
prefDefaults{strcmp('vistaPath',prefNames)} = '';
prefDefaults{strcmp('lastPath',prefNames)} = '';
prefDefaults{strcmp('overlayCombineTransformPaths',prefNames)} = '';
prefDefaults{strcmp('roiTransformPaths',prefNames)} = '';
prefDefaults{strcmp('colormapPaths',prefNames)} = '';

if nargin == 0
  if nargout > 0
    % return arguments
    value = prefNames;
  else
    % print out list of preferences
    for i = 1:length(prefNames)
      %get the preference
      prefValue =  mrGetPref(prefNames{i});
      % print it out
      if isnumeric(prefValue)
	disp(sprintf('%s: %s',prefNames{i},num2str(prefValue)));
      elseif isstr(prefValue)
	disp(sprintf('%s: %s',prefNames{i},prefValue));
      elseif iscell(prefValue)
	mrDisp(sprintf('%s:',prefNames{i}));
	for j = 1:length(prefValue)
	  if isnumeric(prefValue{j})
	    mrDisp(sprintf(' %s',num2str(prefValue{j})));
	  elseif isstr(prefValue{j})
	    mrDisp(sprintf(' ''%s''',prefValue{j}));
	  end
	end
	mrDisp(sprintf('\n'));
      end
    end
  end
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
