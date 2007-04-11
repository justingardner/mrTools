% mrGlobals script
%
% Defines MLR as a global variable.
% If MLR is not yet initialized then do so.
% Runs as a script in the scope of the calling function.
%
% djh 6/2004

global MLR

% If MLR is not yet initialized then do so
if isempty(MLR) || (isfield(MLR,'session') && isempty(MLR.session))
    
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
    MLR.session = session;
    MLR.groups = groups;
    
    % Initialize MLR.views
    MLR.views = {};
    
    % Initialize graph window
    MLR.graphFigure = [];

    % Init the figloc and prefs
    if ~isfield(MLR,'figloc')
      MLR.figloc = [];
    end
    if ~isfield(MLR,'prefs')
      MLR.prefs = [];
    end

    % read the default figlocs
    mrDefaults = loadMrDefaults;
    if isfield(mrDefaults,'figloc') && isempty(MLR.figloc)
      MLR.figloc = mrDefaults.figloc;
    end
    if isfield(mrDefaults,'prefs') && isempty(MLR.prefs)
      MLR.prefs = mrDefaults.prefs;
    end
    
    % load preferences, note that this will override
    % any preferences set my mrDefaults. This is intended
    % behavior. If you use setpref it means you really want it
    % that way. Ones saved in mrDefaults are preferences 
    % the last time you savedMrDefaults
    prefs = getpref('mrLoadRet');
    % if there are some preferences, put them into the MLR variable
    % for easy and quicker access
    if ~isempty(prefs)
      prefnames = fieldnames(prefs);
      for i = 1:length(prefnames)
	MLR.prefs.(prefnames{i}) = getpref('mrLoadRet',prefnames{i});
      end
    end
    
    % Inform user that mrLoadRet has started up
    disp(['mrLoadRet ',num2str(MLR.version),', Matlab ',num2str(matlabVersion)]);
    
    % Clean up
    clear expectedMatlabVersion version matlabVersion session groups mlrVersion mrDefaults
end

