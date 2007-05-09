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
global mrDEFAULTS;

% save figloc
if isfield(mrDEFAULTS,'figloc')
    figloc = mrDEFAULTS.figloc;
else
    figloc = [];
end

if isfield(mrDEFAULTS,'prefs')
    prefs = mrDEFAULTS.prefs;
else
    prefs = [];
end


save '~/.mrDefaults.mat' figloc prefs -V6;
