% mlrMake.m
%
%        $Id:$ 
%      usage: mlrMake()
%         by: justin gardner
%       date: 05/24/18
%    purpose: recompile files in mlr
%       e.g.: mlrMake - will try make all files
%             mlrMake filename.c - will make the named file
%
function retval = mlrMake(varargin)

% check arguments
if ~any(nargin == [0 1])
  help mlrMake
  return
end

% get the mrTools directory
mlrTop = fileparts(fileparts(which('mrLoadRet')));

% find filenames that need compiling
if nargin == 0
  cFiles = findFiles(mlrTop,'.c');
  ccFiles = findFiles(mlrTop,'.cpp');
  compiledFunctionList = {cFiles{:} ccFiles{:}};
else
  compiledFunctionList = findFiles(mlrTop,varargin{1})
end

% set mexopts file
optf = sprintf('-Dchar16_t=uint16_T -f %s',fullfile(mlrTop,'mrUtilities','make','mexopts.sh'));

% list of files that were compiled ok
successfullyCompiled = {};
failedToCompile = {};
for iFile = 1:length(compiledFunctionList)
  % check for file
  if isfile(compiledFunctionList{iFile})
    % display what we are doing
    disp(sprintf('(mlrMake) mex: %s',compiledFunctionList{iFile}));
    % mex the file
    try
      eval(sprintf('mex %s %s',optf,compiledFunctionList{iFile}));
      successfullyCompiled{end+1} = compiledFunctionList{iFile};
    catch
      disp(sprintf('(mlrMake) Error making %s',compiledFunctionList{iFile}));
      failedToCompile{end+1} = compiledFunctionList{iFile};
    end
  else 
    % display that we can't find file
    disp(sprintf('(mlrMake) Could not find file: %s', compiledFunctionList{iFile}));
  end
end

% list files that were not ok
disp(sprintf('(mglMake) Failed to compile'));
for iFile = 1:length(failedToCompile)
  disp(sprintf('%s',failedToCompile{iFile}));
end

% list files that were ok
disp(sprintf('(mglMake) Successfully compiled'));
for iFile = 1:length(successfullyCompiled)
  disp(sprintf('%s',successfullyCompiled{iFile}));
end

%%%%%%%%%%%%%%%%%%%
%    findFiles    %
%%%%%%%%%%%%%%%%%%%
function filenames = findFiles(dirname,matchName)

filenames = {};
d = dir(dirname);
for iFile = 1:length(d)
  if d(iFile).isdir && ~any(strcmp(d(iFile).name,{'.','..'}))
    % recursively search
    subdirMatch = findFiles(fullfile(dirname,d(iFile).name),matchName);
    % concat if not empty
    if ~isempty(subdirMatch)
      filenames = {filenames{:} subdirMatch{:}};
    end
  else
    % if it is an extension we are matching
    if (matchName(1) == '.')
      % get filename parts
      [pathsr,name,ext] = fileparts(d(iFile).name);
      % check extension
      if strcmp(lower(matchName),ext)
	filenames{end+1} = fullfile(dirname,d(iFile).name);
      end
    else
      % check filename
      if strcmp(lower(matchName),lower(d(iFile).name))
	filenames{end+1} = fullfile(dirname,d(iFile).name);
      end
    end
  end
end
