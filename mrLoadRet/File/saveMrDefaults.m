% saveMrDefaults
%
%      usage: saveMrDefaults()
%         by: justin gardner
%       date: 03/17/07
%    purpose: save default positions
%
function retval = saveMrDefaults()

% check arguments
if ~any(nargin == [0])
  help saveMrDefaults
  return
end

% get globals
mrGlobals;

% save figloc
if isfield(MLR,'figloc')
  figloc = MLR.figloc;
  save '~/.mrDefaults.mat' figloc -V6;
end

