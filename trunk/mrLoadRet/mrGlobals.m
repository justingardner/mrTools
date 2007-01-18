% mrGlobals script
%
% Defines MLR as a global variable.
% If MLR is not yet initialized then do so.
% Runs as a script in the scope of the calling function.
%
% djh 6/2004

global MLR

% If MLR is not yet initialized then do so
if isempty(MLR)
    % Check Matlab version number
    expectedMatlabVersion = [7.2];  % Change this after testing Matlab upgrades
    version = ver('Matlab');
    matlabVersion = str2num(version.Version(1:3));        
    if ~ismember(matlabVersion, expectedMatlabVersion);
        mrWarnDlg(['mrLoadRet is intended for Matlab ',num2str(expectedMatlabVersion),...
                   '. You are running Matlab ',version.Version]);
    end
    
    % Load session and groups structures from mrSESSION.mat
    MLR.homeDir = pwd;
    [session, groups, mrLoadRetVersion] = loadSession(MLR.homeDir);
    MLR.version = mrLoadRetVersion;
    MLR.session = session;
    MLR.groups = groups;
    
    % Initialize MLR.views
    MLR.views = {};
    
    % Initialize graph window
    MLR.graphFigure = [];
    
    % Inform user that mrLoadRet has started up
    disp(['mrLoadRet ',num2str(MLR.version),', Matlab ',matlabVersion]);
    
    % Clean up
    clear expectedMatlabVersion version matlabVersion session groups mrLoadRetVersion
end
