% saveSform.m
%
%        $Id$
%      usage: saveSform(sform)
%         by: justin gardner
%       date: 06/20/07
%    purpose: saves sform to mrLoadRet-4.5 mrSession variable and
%             changes header in epis
%
function retval = saveSform(sform)

% check arguments
if ~any(nargin == [1])
  help saveSform
  return
end

% first check this directory for mrSession.m
path = '';
if isfile('mrSession.mat')
  % if there is a session then ask if the user wants to export to this directory
  answer = questdlg(sprintf('Save alignment to %s?',getLastDir(pwd)),'Export');
  if strcmp(answer,'Cancel'),return,end
  if strcmp(answer,'Yes')
    path = pwd;
  end
end

% if path is not set then ask the user to choose a path
if isempty(path)
  [filename path] = uigetfile('*.mat','Choose a mrSession.mat file to save alignment to');
  if filename == 0,return,end
  if ~strcmp(filename,'mrSession.mat')
    mrWarnDlg(sprintf('(saveSform) %s is not an mrSession file',filename));
    return
  end
end

% start a view in the corresponding location
cd(path);
v = newView('Volume');
mrGlobals
% just make sure that the home dir matches
if ~strcmp(MLR.homeDir,path)
  answer = questdlg(sprintf('mrLoadRet is open on sesion %s? Ok to close and open on %s',getLastDir(MLR.homeDir),getLastDir(path)));
  if ~strcmp(answer,'Yes')
    mrWarnDlg(sprintf('(saveSform) Could not open a view to %s',getLastDir(path)));
    deleteView(v);
    return
  end
  % clear MLR and start over
  deleteView(v);
  clear global MLR;
  v = newView('Volume');
end

% ask the user which groups to export to
groupNames = viewGet(v,'groupNames');
for i = 1:length(groupNames)
  paramsInfo{i}{1} = groupNames{i};
  paramsInfo{i}{2} = 1;
  paramsInfo{i}{3} = 'type=checkbox';
  paramsInfo{i}{4} = sprintf('Export alignment to group %s',groupNames{i});
end
params = mrParamsDialog(paramsInfo,'Choose groups to export alignment to');
if isempty(params),return,end

% now go through and update the headers
for iGroup = 1:viewGet(v, 'numberofGroups')
  % see if we are supposed to update the group
  if params.(viewGet(v,'groupName',iGroup))
    for iScan = 1:viewGet(v, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(v, 'niftiHdr', iScan, iGroup);

      % check to see if it is a change (this is just to print user info)
      if ~isfield(curhdr,'sform44') || ~isequal(curhdr.sform44,sform)
	disp(sprintf('Nifti hdr for scan %i in %s group has been modified', iScan, viewGet(v, 'groupName',iGroup)));
      end

      % get the nifti filename
      filename = viewGet(v, 'tseriesPath', iScan, iGroup);
      % check if it is there
      if isfile(filename)
	% load the header
	hdr = cbiReadNiftiHeader(filename);
	% set the sform
	hdr = cbiSetNiftiSform(hdr,sform);
	% and write it back
	hdr = cbiWriteNiftiHeader(hdr,filename);
	% now save it in the session
	v = viewSet(v,'niftiHdr',hdr,iScan,iGroup);
      else
	disp(sprintf('(saveSform) Could not open file %s',filename));
      end
    end
  end
end

% save the session
saveSession
% also remove any base anatomies from mrLastView if it is
% there since those might have a different sform
if isfile('mrLastView.mat')
  disp(sprintf('(saveSform) Removing base anatomies from mrLastView'));
  load mrLastView
  view.baseVolumes = [];
  save mrLastView view viewSettings
end
deleteView(v);

