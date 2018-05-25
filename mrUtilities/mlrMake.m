% mlrMake.m
%
%        $Id:$ 
%      usage: mlrMake()
%         by: justin gardner
%       date: 05/24/18
%    purpose: 
%
function retval = mlrMake()

% check arguments
if ~any(nargin == [0])
  help mlrMake
  return
end

% get the mrTools directory
mlrTop = fileparts(fileparts(which('mrLoadRet')));

% find filenames that need compiling
cFiles = findFiles(mlrTop,'.c');
ccFiles = findFiles(mlrTop,'.cc');
compiledFunctionList = {cFiles{:} ccFiles{:}};

for iFile = 1:length(compiledFunctionList)
  % check for file
  if isfile(compiledFunctionList{iFile})
    % display what we are doing
    disp(sprintf('(mlrMake) mex: %s',compiledFunctionList{iFile}));
    % mex the file
    keyboard
    mex(compiledFunctionList{iFile});
  else 
    % display that we can't find file
    disp(sprintf('(mlrMake) Could not find file: %s', compiledFunctionList{iFile}));
  end
end


%%%%%%%%%%%%%%%%%%%
%    findFiles    %
%%%%%%%%%%%%%%%%%%%
function filenames = findFiles(dirname,extMatch)

filenames = {};
d = dir(dirname);
for iFile = 1:length(d)
  if d(iFile).isdir && ~any(strcmp(d(iFile).name,{'.','..'}))
    % recursively search
    subdirMatch = findFiles(fullfile(dirname,d(iFile).name),extMatch);
    % concat if not empty
    if ~isempty(subdirMatch)
      filenames = {filenames{:} subdirMatch{:}};
    end
  else
    % get filename parts
    [pathsr,name,ext] = fileparts(d(iFile).name);
    % check extension
    if strcmp(lower(extMatch),ext)
      filenames{end+1} = fullfile(dirname,d(iFile).name);
    end
  end
end
