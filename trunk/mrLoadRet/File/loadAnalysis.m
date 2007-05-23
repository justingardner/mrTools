function view = loadAnalysis(view,pathStr,startPathStr)
%
% view = loadAnalysis(view,[pathStr],[startPathStr])
%
% Loads an analysis and adds it to view.analyses.
%
% If pathStr is not specified, prompts user to select a file. If pathStr is
% specified, it is assumed to be relative to startPathStr which is assumed
% (if not specified) to be the current group's data directory. pathStr can
% be a string specifying a single analysis file or it can be a cell array
% of filenames to load multiple analyses at once.
%
% The file must contain a structure or structures, each with the following
% fields:
% - name: string
% - groupName: string group name
% - function: string function name by which it was computed
% - reconcileFunction: string function name that reconciles params
%   and data with tseries file names.
%      [newparams,newdata] = reconcileFunction(groupName,params,data)
% - guiFunction: string function name that allows user to specify
%   or change params.
%      params = guiFunction('groupName',groupName,'params',params)
% - params: structure specifying arguments to function
%      To recompute: view = function(view,params)
% - overlays: struct array of overlays
% - curOverlay: integer corresponding to currently selected overlay
% - date: specifies when it was computed
%
% In addition the 'name' field is set to the variable name to ensure
% consistency.
%
% djh, 1/9/98
% 6/2005, djh, update to mrLoadRet-4.0

mrGlobals

% Path to analyses
if ieNotDefined('startPathStr')
  startPathStr = viewGet(view,'dataDir');
end

% Complete pathStr
if ieNotDefined('pathStr')
    pathStr = getPathStrDialog(startPathStr,'Choose one or more analyses','*.mat','on');
else
  if iscell(pathStr)
    for p=1:length(pathStr)
      pathStr{p} = fullfile(startPathStr,pathStr{p});
    end
  else
    pathStr = {fullfile(startPathStr,pathStr)};
  end
end
if isempty(pathStr)
  return
end

if ~iscell(pathStr)
    pathStr = {pathStr};
end

% Load the file. Loop through the variables that were loaded and add
% each of them as a new analysis, setting analysis.fieldnames as we go.
for p = 1:length(pathStr)
	if exist(pathStr{p},'file')
		h = mrMsgBox(['Loading analysis: ',pathStr{p},'. Please wait']);
		s = load(pathStr{p});
		varNames = fieldnames(s);
		analysis = eval(['s.',varNames{1}]);
		analysis.name = varNames{1};
		% Add it to the view
		view = viewSet(view,'newAnalysis',analysis);
		mrCloseDlg(h);
	else
		mrWarnDlg(['Analysis ',pathStr{p},' not found.']);
	end
end

return;

% Test/debug
view = loadAnalysis(MLR.views{1},'co');
view = loadAnalysis(MLR.views{1},{'co','amp','ph'});
view = loadAnalysis(MLR.views{1});

