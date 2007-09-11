% loadMrDefaults
%
%      usage: loadMrDefaults()
%         by: justin gardner
%       date: 03/17/07
%    purpose: load default positions
%
function mrDefaults = loadMrDefaults()

mrDefaults = [];

% check arguments
if ~any(nargin == [0])
  help loadMrDefaults
  return
end

% load the defaults or set them to default values
if isfile('~/.mrDefaults.mat')
  mrDefaults = load('~/.mrDefaults.mat');
else
    mrDefaults.prefs = [];
    mrDefaults.figloc = [];
end

% If any of the defaults are unset, them set them
if ~isfield(mrDefaults.prefs,'interpMethod') | isempty(mrDefaults.prefs.interpMethod)
    mrDefaults.prefs.interpMethod = 'nearest';
end
if ~isfield(mrDefaults.prefs,'overwritePolicy') | isempty(mrDefaults.prefs.overwritePolicy)
    mrDefaults.prefs.overwritePolicy = 'ask';
end
if ~isfield(mrDefaults.prefs,'site') | isempty(mrDefaults.prefs.site)
    mrDefaults.prefs.site = 'NYU';
end
if ~isfield(mrDefaults.prefs,'verbose') | isempty(mrDefaults.prefs.verbose)
    mrDefaults.prefs.verbose = 'Yes';
end
if ~isfield(mrDefaults.prefs,'maxBlocksize') | isempty(mrDefaults.prefs.maxBlocksize)
    mrDefaults.prefs.maxBlocksize = 250000000;
end
if ~isfield(mrDefaults.prefs,'volumeDirectory') | isempty(mrDefaults.prefs.volumeDirectory)
    mrDefaults.prefs.volumeDirectory = pwd;
end
if ~isfield(mrDefaults.prefs,'niftiFileExtension') | isempty(mrDefaults.prefs.niftiFileExtension)
    mrDefaults.prefs.niftiFileExtension = '.img';
end
if ~isfield(mrDefaults.prefs,'roiCacheSize') | isempty(mrDefaults.prefs.roiCacheSize)
  mrDefaults.prefs.roiCacheSize = 50;
% this is added since when we went to single roi caches
% we need to have a larger roi cache size--this will convert
% everyone who had the default 3 to 50
elseif (mrDefaults.prefs.roiCacheSize == 3)
  mrDefaults.prefs.roiCacheSize = 50;
end  
if ~isfield(mrDefaults.prefs,'overlayCacheSize') | isempty(mrDefaults.prefs.overlayCacheSize)
  mrDefaults.prefs.overlayCacheSize = 20;
end
if ~isfield(mrDefaults.prefs,'baseCacheSize') | isempty(mrDefaults.prefs.baseCacheSize)
  mrDefaults.prefs.baseCacheSize = 20;
end

% check for any figloc that are strange
if ~isempty(mrDefaults.figloc)
  figlocNames = fieldnames(mrDefaults.figloc);
  for i = 1:length(figlocNames)
    % if it isn't of length four then just remove it 
    if length(mrDefaults.figloc.(figlocNames{i})) ~= 4
      mrDefaults.figloc = rmfield(mrDefaults.figloc,figlocNames{i});
    end
  end
end
