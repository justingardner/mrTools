% interrogatorFolderPlugin.m
%
%        $Id$ 
%      usage: interrogatorFolderPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%    purpose: gathers a list of all interrogators in the interrogator folder
%             and adds it to the default interrogator list
%

function retval = interrogatorFolderPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help interrogatorFolderPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(interrogatorFolderPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    %get interrogators in the interrogatorFunctions directory
    interrogatorsFolder = [fileparts(which('interrogatorFolderPlugin')) '/InterrogatorFunctions/'];
    interrogatorFiles =  dir([interrogatorsFolder '*.m']);
    if ~isempty(interrogatorFiles)
      interrogatorList = cell(1,length(interrogatorFiles));
      for iFile=1:length(interrogatorFiles)
         interrogatorList{iFile} = stripext(interrogatorFiles(iFile).name);
      end
      % install default interrogators
      % that will show up when you do /Edit/Overlay
      mlrAdjustGUI(thisView,'add','interrogator',interrogatorList);
    else
      disp('(interrogatorFolderPlugin) No interrogator function found in folder');
    end
    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Scans the mrLoadRet/Plugin/interrogatorFolder/interrogatorFunctions folder at MLR startup and adds all the functions to the default interrogator list';
 otherwise
   disp(sprintf('(interrogatorFolderPlugin) Unknown command %s',action));
end

