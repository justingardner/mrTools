function saveSession(verbose)
% saveSession([verbose])
%
% verbose: ask before overwriting existing mrSession file
%    default: viewGet(MLR.views{1},'pref','verbose') or if that's not set 0
% djh 5/2005

mrGlobals;

if ieNotDefined('verbose')
    verbose = viewGet([],'pref','verbose');
    if isempty(verbose)
        verbose = 0;
    end
end

pathStr = fullfile(MLR.homeDir,'mrSession.mat');

if verbose
    if exist(pathStr,'file')
        questionString = 'mrSession.mat already exists. Do you want to overwrite it?';
        buttonName = questdlg(questionString, 'Warning', 'Yes', 'No', 'No');
        pause(.1);  % Prevent hanging
        if strcmp(buttonName, 'No')
            return
        end
    end
else
    disp('Warning: overwriting mrSession.mat');
end

session = MLR.session;
groups = MLR.groups;
save(pathStr,'session','groups');
