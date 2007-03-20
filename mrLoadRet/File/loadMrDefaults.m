% loadMrDefaults
%
%      usage: loadMrDefaults()
%         by: justin gardner
%       date: 03/17/07
%    purpose: load default positions
%
function mrDefaults = loadMrDefaults()

mrDefaults = [];

% check arguments
if ~any(nargin == [0])
  help loadMrDefaults
  return
end

% load the defaults
if isfile('~/.mrDefaults.mat')
  mrDefaults = load('~/.mrDefaults.mat');
end

% check for any figloc that are strange
if isfield(mrDefaults,'figloc')
  figlocNames = fieldnames(mrDefaults.figloc);
  for i = 1:length(figlocNames)
    % if it isn't of lenght four then just remove it 
    if length(mrDefaults.figloc.(figlocNames{i})) ~= 4
      mrDefaults.figloc = rmfield(mrDefaults.figloc,figlocNames{i});
    end
  end
end
