% makeEmptyMLRDir.m
%
%        $Id$ 
%      usage: makeEmptyMLRDir(dirname)
%         by: justin gardner
%       date: 01/08/09
%    purpose: Makes an empty MLR directory that can be used
%             to import groups into
%
function retval = makeEmptyMLRDir(dirname)

% check arguments
if ~any(nargin == [1])
  help makeEmptyMLRDir
  return
end

directories = {'Raw' fullfile('Raw','TSeries') 'Anatomy' 'Etc' 'Anal'};

% dirname is a file, abort
if isfile(dirname)
  mrWarnDlg(sprintf('(makeEmptyMLRDir) %s exists. Aborting',dirname));
  return
end

% existing directory. Ask user what to do
if isdir(dirname)
  if ~askuser(sprintf('(makeEmptyMLRDir) Directory %s exists, continue?',dirname))
    return
  end
else
  % make the directory
  mkdir(dirname)
end

% check for existing mrSession
if isfile(fullfile(dirname,'mrSession.mat'))
  disp(sprintf('(makeEmptyMLRDir) %s already has an mrSession.mat. Aborting',dirname));
  return
end

% make the directories
for dirnum = 1:length(directories)
  mkdir(fullfile(dirname,directories{dirnum}));
end



% set up session variable
magnet = mrGetPref('magnet');
coil = mrGetPref('coil');
pulseSequence = mrGetPref('pulseSequence');
description = '';
subject = '';
operator = '';

% setup params dialog
paramsInfo = {};
paramsInfo{end+1} = {'description',description,'type=string','Description of the session. Can be anything to help you remember what the session was'};
paramsInfo{end+1} = {'subject',subject,'type=string','Subject ID. Use an identifier that does not break the subject confidentiality'};
paramsInfo{end+1} = {'operator',operator,'type=string','Person who operated the scanner'};
paramsInfo{end+1} = {'magnet',magnet,'Choose which magnet you scanned on'};
paramsInfo{end+1} = {'coil',coil,'Choose which coil you used'};
paramsInfo{end+1} = {'pulseSequence',pulseSequence,'Choose which pulse sequence you scanned with'};
paramsInfo{end+1} = {'pulseSequenceText','','Optional: enter some text to describe or qualify the pulseSequence text. This will get appended to the pulseSequence name chosen above'};

% get the session params from the user
sessionParams = mrParamsDialog(paramsInfo,'Initialize session for mrTools');

% create session variable
session.mrLoadRetVersion = mrLoadRetVersion;
session.description = sessionParams.description;
session.subject = sessionParams.subject;
session.operator = sessionParams.operator;
session.magnet = sessionParams.magnet;
session.coil = sessionParams.coil;
session.protocol = sprintf('%s: %s',sessionParams.pulseSequence,sessionParams.pulseSequenceText);

% create groups variables
groups.name = 'Raw';
groups.scanParams = [];
[tf groups] = isgroup(groups);

% save the mrSession
eval(sprintf('save %s session groups',fullfile(dirname,'mrSession.mat')));

