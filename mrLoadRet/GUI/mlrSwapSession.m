% mlrSwapSession.m
%
%        $Id:$ 
%      usage: [lastSession currentSession] = mlrSwapSession()
%         by: justin gardner
%       date: 01/03/17
%    purpose: Temporarily swap to a different session. This is used
%             when you have mrLoadRet open to one session, but you temporarily
%             need to access a different session. It will store all the
%             global variables temporarily then you can open a new view 
%             in a different session (make sure that you don't try
%             to access the old session while you are swapped out). THen
%             you can close and return to the old session.
%
%             lastSession is the name of the session you swapped out of
%             currentSession is the name of the session you swapped into
%
%             Simple usage:
%
%             firstSession = mlrSwapSession('path/to/new/session');;
%
%             % do some stuff with the other session
%             v = newView;
%             ...
%
%             % return to this session
%             secondSession = mlrSwapSession(firstSession);
%
%             You can continue to go back and forth between the two sessions
%             by calling mlrSwapSession with whichever session you want
%             to use
%
%             Without any arguments, will swap back to the last session
%             you came from:
%
%             mlrSwapSession;
%             cd('/path/to/other/session');
%             v = newView;
%             mlrSwapSession;
%
%             You can also swap between multiple session:
%
%             mlrSwapSession('/path/to/first/session');
%             mlrSwapSession('/path/to/second/session');
%             mlrSwapSession('/path/to/third/session');
%
%
function [currentSessionName sessionName] = mlrSwapSession(sessionName)

% check arguments
if ~any(nargin == [0 1])
  help mlrSwapSession
  return
end

% default arguments
if nargin < 1
  sessionName = [];
end

% init global variable for storing session data
global gSwapSession
if isempty(gSwapSession)
  gSwapSession.lastSession = [];
  gSwapSession.sessionNames = {};
  gSwapSession.globals = {};
  gSwapSession.currentDir = {};
end

% get current session name
if mlrIsRunning(0)
  v = newView;
  currentSessionName = viewGet(v,'homeDir');
  deleteView(v);
else
  currentSessionName = [];
  % check to see if we have a stored session with current path
  storedSessionNum = first(find(strcmp(pwd,gSwapSession.sessionNames)));
  % clear it
  if ~isempty(storedSessionNum)
    gSwapSession.globals{storedSessionNum} = [];
    gSwapSession.currentDir{storedSessionNum} = [];
  end
end

% if sessionName is not set then get the last session that was swapped
if isempty(sessionName)
  sessionName = gSwapSession.lastSession;
end

if ~isempty(sessionName) && ~isdir(sessionName)
  disp(sprintf('(mlrSwapSession) Could not find sesion directory; %s',sessionName));
  return
end

% save what session we are on now
gSwapSession.lastSession = currentSessionName;
if ~isempty(currentSessionName)
  % get the session number to save under (i.e. if there is one
  % that was already saved then use that one)
  saveSessionNum = first(find(strcmp(currentSessionName,gSwapSession.sessionNames)));
  % otherwise add to the list of saved variables
  if isempty(saveSessionNum)
    saveSessionNum = length(gSwapSession.sessionNames)+1;
  end
  % save the session name
  gSwapSession.sessionNames{saveSessionNum} = currentSessionName;
  % and save the globals
  gSwapSession.globals{saveSessionNum} = storeAndClearMLRGlobals;
  % and save the current directory
  gSwapSession.currentDir{saveSessionNum} = pwd;
end

% if session name is empty, then we are done
if isempty(sessionName),return,end

% if not, then check to see if it is a stored session
loadSessionNum = first(find(strcmp(sessionName,gSwapSession.sessionNames)));
if ~isempty(loadSessionNum)
  % get the names of the stored variables
  if ~isempty(gSwapSession.globals{loadSessionNum})
    varNames = fieldnames(gSwapSession.globals{loadSessionNum});
    % and restore them
    for iVar = 1:length(varNames)
      eval(sprintf('global %s;',varNames{iVar}));
      eval(sprintf('%s = gSwapSession.globals{%i}.%s;',varNames{iVar},loadSessionNum,varNames{iVar}));
    end
    % and restore directory
    cd(gSwapSession.currentDir{loadSessionNum});
  end
else
  % switch to session
  cd(sessionName);
end
  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    storeAndClearMLRGlobals    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function globalVars = storeAndClearMLRGlobals

% list of global variables that need to be swapped out/in
mlrGlobalVariables = {'mrDEFAULTS','MLR'};

% get the global variables, save them in gSwapSession and clear them
for iVar = 1:length(mlrGlobalVariables)
  % get this global name
  thisGlobal = mlrGlobalVariables{iVar};
  % get the global variable
  eval(sprintf('global %s;',thisGlobal));
  % and keep its value
  eval(sprintf('globalVars.%s = %s;',thisGlobal,thisGlobal));
  % and clear it
  eval(sprintf('clear global %s;',thisGlobal));
end
