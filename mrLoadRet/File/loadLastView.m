% loadLastView.m
%
%        $Id:$ 
%      usage: loadLastView()
%         by: justin gardner
%       date: 04/02/15
%    purpose: This function became necessary since in version 8.2, mathworks 
%             inexplicably busted figure handles. Now they are handles and
%             if you load an old variable with a handle, it causes massive
%             convulsions in matlab. It also tries to pop up a figure - which
%             is not what we want - eeks. Anyway, the solution is to load
%             part of the mrLastView variable (viewSettings) and check
%             its version - if it is an old style version then warn the user
%             that we can't load
%
function [v viewSettings]  = loadLastView(filename)

v = [];
viewSettings = [];

% check arguments
if ~any(nargin == [0 1])
  help loadLastView
  return
end

if nargin == 0;
  filename = 'mrLastView.mat';
end
filename = setext(filename,'mat');

if ~isfile(filename);
  mrWarnDlg('(loadLastView) Could not find %s',filename);
  return
end

% new version of matlab needs to check if we are trying to load an old file
if mlrGetMatlabVersion >= 8
  % load the viewSettings part and look to see if it is a new file
  check = load(filename,'viewSettings');
  if isfield(check,'viewSettings') && isfield(check.viewSettings,'version') && (check.viewSettings.version>=2.0)
    % then we are ok, load the view part
    l = load(filename,'view');
    % return them both
    if nargout == 1
      % return as single argument
      v = l;
    else
      % or as two
      v = l.view;
      viewSettings = check.viewSettings;
    end
    return
  else
    mrWarnDlg(sprintf('(loadLastView) The mrLastView found: %s is from an older version of matlab which allowed saving figure handles. The geniuses at Mathworks have busted that, so loading this file will no longer work. Moving mrLastView.mat to mrLastView.mat.old You will lose any rois that were loaded but not saved and mrLoadRet will start up without bases and analyses loaded. If you really need what was in the viewer we suggest running on an earlier version of matlab - you just then need to copy mrLastView.mat.old back to mrLastView.mat, open mrLoadRet and then quit - this will save the file back w/out the offending figure handles. Send complaints to Mathworks!',filename),'Yes');
    movefile(filename,sprintf('%s.old',filename));
    return
  end
else
  % otherwise just load
  l = load('mrLastView');
  if nargout == 1
    v = l;
  else
    v = l.view;
    viewSettings = l.viewSettings;
  end
end


