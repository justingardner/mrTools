% mrSaveView.m
%
%        $Id: mrQuit.m 1942 2010-12-16 18:14:41Z julien $ 
%      usage: mrSaveView(v)
%         by: justin gardner, taken out from mrQuit by julien besle
%       date: 07/11/08, 2011/08/05
%    purpose: saves view and view settings in session directory

function mrSaveView(v)

% remember figure location
fig = viewGet(v,'fignum');
if ~isempty(fig)
  mrSetFigLoc('mrLoadRetGUI',get(fig,'Position'));
end
% remember GUI settings
viewSettings.curBase = viewGet(v,'curBase');
viewSettings.rotate = viewGet(v,'rotate');
viewSettings.curScan = viewGet(v,'curScan');
viewSettings.curSlice = viewGet(v,'curSlice');
viewSettings.curGroup = viewGet(v,'curGroup');
viewSettings.sliceOrientation = viewGet(v,'sliceOrientation');
viewSettings.overlayMin = viewGet(v,'overlayMin');
viewSettings.overlayMax = viewGet(v,'overlayMax');
viewSettings.alpha = viewGet(v,'alpha');
viewSettings.showROIs = viewGet(v,'showROIs');
viewSettings.labelROIs = viewGet(v,'labelROIs');
viewSettings.roiGroup = viewGet(v,'roiGroupNames');
homeDir = viewGet(v,'homeDir');
try
  disppercent(-inf,sprintf('(mrSaveView) Saving %s/mrLastView',homeDir));
        % save the view in the current directory
  view = v;
  if getfield(whos('view'),'bytes')<2e9
    eval(sprintf('save %s view viewSettings -V6;',fullfile(homeDir,'mrLastView')));
  else
    mrWarnDlg('(mrSaveView) Variable view is more than 2Gb, using option -v7.3 to save');
    eval(sprintf('save %s view viewSettings -v7.3;',fullfile(homeDir,'mrLastView')));
  end
  % save .mrDefaults in the home directory
  disppercent(inf);
catch
  disppercent(inf);
  mrErrorDlg('(mrQuit) Could not save mrLastView.mat');
end
