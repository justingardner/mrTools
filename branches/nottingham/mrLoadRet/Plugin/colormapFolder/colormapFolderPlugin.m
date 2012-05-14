% colormapFolderPlugin.m
%
%        $Id$ 
%      usage: colormapFolderPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%    purpose: gathers a list of all colormaps in the colormap folder
%             and adds it to the default colormaps list
%

function retval = colormapFolderPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help colormapFolderPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(colormapFolderPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    %get colormaps in the colormapFunctions directory
    colorMapsFolder = [fileparts(which('colormapFolderPlugin')) '/ColormapFunctions/'];
    colorMapFiles =  dir([colorMapsFolder '*.m']);
    if ~isempty(colorMapFiles)
      colorMapList = cell(1,length(colorMapFiles));
      for iFile=1:length(colorMapFiles)
         colorMapList{iFile} = stripext(colorMapFiles(iFile).name);
      end
      % install default colormaps
      % that will show up when you do /Edit/Overlay
      mlrAdjustGUI(thisView,'add','colormap',colorMapList);
    else
      disp('(colormapFolderPlugin) No colormap function found in folder');
    end
    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Scans the mrLoadRet/Plugin/colormapFolder/colormapFunctions folder at MLR startup and adds all the functions to the default colormap list';
 otherwise
   disp(sprintf('(colormapFolderPlugin) Unknown command %s',action));
end

