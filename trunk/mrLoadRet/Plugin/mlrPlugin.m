% mlrPlugin.m
%
%        $Id:$ 
%      usage: mlrPlugin(<v>)
%         by: justin gardner
%       date: 11/24/10
%    purpose: With no arguments, brings up plugin selection dialog. This function will
%             look in the directory mrLoadRet/Plugin for plugins (i.e. directories
%             that have a function names dirnamePlugin.m), you can also specify
%             a comma-delimited list of alternative plugin directories to also 
%             search for plugins with the preference 'pluginPaths' - i.e. if you had
%             a group plugin path you could do: mrSetPref('pluginPaths','/path/to/group/plugins');
% 
%             It saves in the preference 'selectedPlugins', a cell array with the names of 
%             allc urrently selected plugins.
%
%             With a view argument initializes all plugins that the user has selected for that view
%
function retval = mlrPlugin(v)

% check arguments
if ~any(nargin == [0 1])
  help mlrPlugin
  return
end

% default to view is empty
if nargin == 0,v = [];end

% check Plugin directory for possible plugins
% first, find where this function lives. Under that
% should be directories for plugins that are distributed
% with MLR
mlrPluginPath = fileparts(which('mlrPlugin'));
if ~isdir(mlrPluginPath)
  disp(sprintf('(mlrPlugin) Could not find default plugin directory'));
  return
end
% get the plugins from the default location within the MLR distribution
plugins = getPlugins(mlrPluginPath);

% check any directories found in the alternate plugin path
pluginPaths = commaDelimitedToCell(mrGetPref('pluginPaths'));
for i = 1:length(pluginPaths)
  plugins = getPlugins(pluginPaths{i},plugins);
end

% no plugins - display message and return
if isempty(plugins)
  disp(sprintf('(mlrPlugin) No plugin directories found. Plugins should be in directory mrLoadRet/Plugins'));
  return
end

% get which plugins user has selected on, this should be a list
% of plugin names
selectedPlugins = mrGetPref('selectedPlugins');
for i = 1:length(selectedPlugins)
  % get which number plugins are selected
  selected = find(strcmp(selectedPlugins{i},{plugins.name}));
  % and select them.
  if length(selected) == 1
    plugins(selected).selected = true;
  elseif length(selected)>1
    disp(sprintf('(mlrPlugin) Multiple (%i) plugins names %s were found. Only selecting the first one (%s).',length(selected),selectedPlugins{i},plugins(selected(1)).path));
    plugins(selected(1)).selected = true;
  end
end

% when run with no arguments...
if nargin == 0
  % keep the list of which plugins were selected
  previousSelectedPlugins = selectedPlugins;
  % put up a selection dialog box
  for i = 1:length(plugins)
    paramsInfo{i} = {plugins(i).name,plugins(i).selected,'type=checkbox',plugins(i).help};
  end
  params = mrParamsDialog(paramsInfo,'Choose plugins you wish to install');
  % if user didn't hit cancel, then save the new choices
  if ~isempty(params)
    % get user choices
    selectedPlugins = {};
    for i = 1:length(plugins)
      if params.(plugins(i).name)
	selectedPlugins{end+1} = plugins(i).name;
      end
    end
    if ~isempty(getMLRView) && ~isequal(previousSelectedPlugins,selectedPlugins)
      mrWarnDlg(sprintf('(mlrPlugin) Restart of MLR is required for change in plugin list to take effect'));
    end
    mrSetPref('selectedPlugins',selectedPlugins);
  end
% with a view argument, then run plugins
else
  % get which plugins are selected
  selectedPlugins = find([plugins.selected]);

  % and run their install command
  for i = 1:length(selectedPlugins)
    disp(sprintf('(mlrPlugin) Installing plugin %s',plugins(selectedPlugins(i)).name));
    eval(plugins(selectedPlugins(i)).installCommand);
  end
end
  

%%%%%%%%%%%%%%%%%%%%
%    getPlugins    %
%%%%%%%%%%%%%%%%%%%%
function plugins = getPlugins(pluginPath,plugins)

% empty plugins directory to start with
if nargin == 1,plugins = [];end

% check for plugin function in each directory
pluginDir = dir(pluginPath);
for i = 1:length(pluginDir)
  if pluginDir(i).isdir && (length(pluginDir(i).name) > 1) && (pluginDir(i).name(1)~='.')
    % check for proper file, that is there should be a file
    % called directoryPlugin.m which can be run to get help
    % info and install the plugin
    pluginName = sprintf('%sPlugin',pluginDir(i).name);
    pluginFullName = fullfile(pluginPath,pluginDir(i).name,sprintf('%s.m',pluginName));
    pluginFullName = mlrReplaceTilde(pluginFullName);
    if isfile(pluginFullName)
      % Make sure plugin function exists
      if isempty(which(pluginName))
	disp(sprintf('(mlrPlugin) Plugin function %s does not exist on path',pluginFullName));
	continue
      end
      % make sure that calling the function returns the correct one on the path
      if ~isequal(pluginFullName,which(pluginName))
	keyboard
	% just give a warning
	disp(sprintf('(mlrPlugin) Plugin function %s found in %s but should be in %s',pluginName,which(pluginName),pluginFullName));
      end
      % now keep info about the plugin
      plugins(end+1).name = pluginDir(i).name;
      plugins(end).path = fileparts(pluginFullName);
      plugins(end).help = eval(sprintf('%s(''help'')',pluginName));
      % make sure we got a string back for the help
      if ~isstr(plugins(end).help)
	disp(sprintf('(mlrPlugin) Plugin %s did not return help string correctly',pluginName));
	plugins(end).help = 'No help information for plugin';
      end
      plugins(end).installCommand = sprintf('%s(''install'',v);',pluginName);
      plugins(end).selected = 0;
    else
      disp(sprintf('(mlrPlugin) Plugin function %s not found in %s',pluginName,fileparts(pluginFullName)));
    end
  end
end
