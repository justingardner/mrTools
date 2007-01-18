function saveSession(verbose)
% saveSession([verbose])
%
% verbose: ask before overwriting existing mrSession file
%    default: getpref('mrLoadRet','verbose') or that's not set 0
% djh 5/2005

mrGlobals;

if ~exist('verbose','var')
    if ispref('mrLoadRet','verbose')
        verbose = getpref('mrLoadRet','verbose');
    else
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
    