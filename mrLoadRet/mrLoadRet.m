% mrLoadRet, version 4.0 
%
% AUTHOR:   DJH
% PURPOSE:  Everything.
% DATE:     June, 2004
%
% This is the default mrLoadRet script that simply opens an
% inplaneView window.  You are encouraged to put a copy of this
% file in each of your data directories and modify the copied
% version of this script to customize it for that data set.
% Examples of a number of possible customizations are commented
% below.

function [v]= mrLoadRet()

if ~isfile('mrSession.mat')
  disp('(mrLoadRet) No mrSession.mat found in current directory');
  return
end

% Define and initialize global variable MLR.
mrGlobals

% check to make sure we are not being run from another data directory
% since mrLoadRet is designed only to be run on one data directory at a time
if isempty(strfind(pwd,viewGet([],'homeDir')))
  % check to see if there is a session mismatch
  [thisSession thisGroups] = loadSession(pwd);
  if ~isempty(thisSession) && (~isequal(thisSession,MLR.session) || ~isequal(thisGroups,MLR.groups))
    disp(sprintf('(mrLoadRet) Current path: %s does not match',pwd));
    disp(sprintf('(mrLoadRet) homeDir: %s in MLR global',viewGet([],'homeDir')));
    disp(sprintf('(mrLoadRet) If you are trying to run two mrLoadRet sessions on different'));
    disp(sprintf('(mrLoadRet) datasets, you should run two separate matlab processes instead'));
    return
  end
end

% Open inplane window
v = mrOpenWindow('Volume');
