% mrGlobals script
%
% Defines MLR as a global variable.
% If MLR is not yet initialized then do so.
% Runs as a script in the scope of the calling function.
%
% djh 6/2004

global MLR
global mrDEFAULTS

% If MLR is not yet initialized then do so
if isempty(MLR) || (isfield(MLR,'session') && isempty(MLR.session))

    % read the preferences and figlocs
    mrDEFAULTS = loadMrDefaults;

    % Check Matlab version number
    [mlrVersion, expectedMatlabVersion] = mrLoadRetVersion;
    version = ver('Matlab');
    matlabVersion = str2num(version.Version(1:3));
    if ~ismember(matlabVersion, expectedMatlabVersion);
        mrWarnDlg(['mrLoadRet is intended for Matlab ',num2str(expectedMatlabVersion),...
            '. You are running Matlab ',version.Version]);
    end

    % Initialize MLR
    MLR.version = mlrVersion;
    MLR.homeDir = pwd;

    % Load session and groups structures from mrSESSION.mat
    [session, groups] = loadSession(MLR.homeDir);
    % check session
    if isempty(session)
      disp(sprintf('(mrGlobals) Could not find mrSession in %s',MLR.homeDir));
    end
    MLR.session = session;
    MLR.groups = groups;

    % Initialize MLR.views
    MLR.views = {};

    % Initialize graph window
    MLR.graphFigure = [];

    % setup caches
    MLR.caches = {};
    
    % Inform user that mrLoadRet has started up
    oneTimeWarning('mrLoadRetVersion',['(mrGlobals) mrLoadRet ',num2str(MLR.version),', Matlab ',num2str(matlabVersion)],1);

    % Clean up
    clear expectedMatlabVersion version matlabVersion session groups mlrVersion 
end

