% mlrAnatDBPreferences.m
%
%        $Id:$ 
%      usage: mlrAnatDBPreferences()
%         by: justin gardner
%       date: 06/30/15
%    purpose: Bring up dialog box to set preference for mlrAnatDB. Note that this
%             is setup to take (but not use) arguments since it can be called from a callback
%
function mlrAnatDBPreferences(varargin)

% get repo locations
centralRepo = mrGetPref('mlrAnatDBCentralRepo');
localRepoTop = mrGetPref('mlrAnatDBLocalRepo');

% set defaults
if isempty(centralRepo),centralRepo = '';end
if isempty(localRepoTop),localRepoTop = '~/data/mlrAnatDB';end

% setup params info for mrParamsDialog
paramsInfo = {...
    {'mlrAnatDBCentralRepo',centralRepo,'Location of central repo, Typically on a shared server with an https address (or via ssh), but could be on a shared drive in the file structure.'}...
    {'mlrAnatDBLocalRepo',localRepoTop,'Location of local repo which is typically under a data directory - this will have local copies of ROIs and other data but can be removed the file system as copies will be stored in the central repo'}...
};

% and display the dialog
params = mrParamsDialog(paramsInfo);

% save params, if user did not hit cancel
if ~isempty(params)
  mrSetPref('mlrAnatDBCentralRepo',params.mlrAnatDBCentralRepo,false);
  mrSetPref('mlrAnatDBLocalRepo',params.mlrAnatDBLocalRepo,false);
end

