% mlrPath.m
%
%        $Id:$ 
%      usage: mlrPath
%         by: justin gardner
%       date: 09/09/14 (from nott.m)
%    purpose: allows happy disambiguation of paths for people
%             with mrTools and vista installed. Examines what
%             paths you have installed and will shut down
%             conflicting paths from vista so that you can run
%             mrTools - also, will bring back those paths when
%             you quit and go back
%
function mlrPath(switchTo)

% try to get paths for vista and for mlr
vistaRoot = getRootPath('vista');
mlrRoot = getRootPath('mlr');

% try to guess vista from mlr if we didn't get it
if isempty(vistaRoot) && ~isempty(mlrRoot)
  vistaRoot = fullfile(fileparts(mlrRoot),'vistasoft');
  if ~isdir(vistaRoot), vistaRoot = []; end
end

% try to guess mlr from vista if we didn't get it
if isempty(mlrRoot) && ~isempty(vistaRoot)
  mlrRoot = fullfile(fileparts(vistaRoot),'mrTools');
  if ~isdir(mlrRoot), mlrRoot = []; end
end

% save as prefs
if ~isempty(mlrRoot), mrSetPref('mlrPath',mlrRoot);end
if ~isempty(vistaRoot), mrSetPref('vistaPath',vistaRoot);end

% get which one we are using
whichMLR = fileparts(fileparts(which('mrAlign')));

% decide which one to switch to.
if (nargin <= 0) || isempty(switchTo)
  if strcmp(whichMLR,vistaRoot) 
    switchTo = mlrRoot;
  elseif strcmp(whichMLR,mlrRoot) 
    switchTo = vistaRoot;
  end
end

% get the top level names of paths
[dump mlrRootName] = fileparts(mlrRoot);
if isempty(mlrRootName),mlrRootName = mlrRoot;end
mlrRootName = lower(mlrRootName);
[dump vistaName] = fileparts(mlrRoot);
if isempty(vistaName),vistaName = vistaRoot;end
vistaName = lower(vistaName);
[dump switchToName] = fileparts(switchTo);
if isempty(switchToName),switchToName = switchTo;end
switchToName = lower(switchToName);

% unknown switchTo
if ~any(strcmp(switchToName,{mlrRootName,vistaName}))
  disp(sprintf('(mlrLife) Currently using: %s',whichMLR));
  return
end  

% display what we are going to do
if verbose
  disp(sprintf('(vista) Switching from %s to %s',whichMLR,switchTo));
end

% first remove everyone
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath(mlrPath));
rmpath(genpath(vistaPath));
warning('on','MATLAB:rmpath:DirNotFound')

% switch path
if strcmp(lower(switchToName),'mrtools')
  % add mrTools
  addpath(genpath(mlrRoot));
  % selectively add some paths from vista
  pathNames = {'mrDiffusion','mrMesh','utilities','mrBOLD/Utilities','fileFilters','mrAnatomy','mrBOLD/UI','external/pyrTools'};
  for i = 1:length(pathNames)
    if isdir(pathNames{i})
      addpath(genpath(fullfile(vistaRoot,pathNames{i})));
    end
  end
elseif strcmp(lower(switchTo),'vista')
  % add vista
  addpath(genpath(vistaRoot));
  % add just the top level directory for mlr (so that we can have this
  % function
  addpath(mlrRoot);
end

%%%%%%%%%%%%%%%%%%%%%
%    getRootPath    %
%%%%%%%%%%%%%%%%%%%%%
function rootPath = getRootPath(packageName)

% get the function that returns the path of the package and pref name
functionName = sprintf('%sRootPath',packageName);
prefName = sprintf('%sPath',packageName);

% see if that function is on the path
if exist(functionName)
  rootPath = eval(functionName);
else
  rootPath = [];
  % try to get the path
  rootPath = mrGetPref(prefName);
end

