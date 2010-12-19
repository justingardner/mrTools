% mrQuit.m
%
%        $Id$ 
%      usage: mrQuit(<saveMrLastView>,<v>)
%         by: justin gardner
%       date: 07/11/08
%    purpose: Function used to quit the mrLoadRet viewer. This is called from
%             the quit menu item, but if you call from the command line directly without 
%             any arguments it will close all open views and the viewer. Note that if you
%             have multiple viewers open on the data set, it will use the last view
%             opened to save settings in mrLastView:
%
%             mrQuit
%
%             Set saveMrLastView to 0 if you do not want to save a mrLastView. i.e.
%
%             mrQuit(0);
%
%
function retval = mrQuit(saveMrLastView,v)

% check arguments
if ~any(nargin == [0 1 2])
  help mrQuit
  return
end

mrGlobals;

% look for the open figure view
if ieNotDefined('v')
  v = [];
  % go through and look for the view with the figure
  if isfield(MLR,'views')
    for i = 1:length(MLR.views)
      if ~isempty(MLR.views{i}) & isfield(MLR.views{i},'figure') & ~isempty(MLR.views{i}.figure)
	if isempty(v)
	  v = MLR.views{i};
	else
	  disp(sprintf('(mrQuit) Multiple open views found. Using last view opened for mrLastView'))
	  v = MLR.views{i};
	end
      end
    end
  end
end

% if we are not saving, then set the view to save from to empty
if ~ieNotDefined('saveMrLastView') & ~saveMrLastView
  v = [];
end

if isfield(MLR,'views') && ~isempty(MLR.views)
  if ~isempty(v)
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
  end
  % close graph figure, remembering figure location
  if ~isempty(MLR.graphFigure)
    mrSetFigLoc('graphFigure',get(MLR.graphFigure,'Position'));
    close(MLR.graphFigure);
    MLR.graphFigure = [];
  end
  % close view figures
  % keep a local copy of everything since
  % the last view that is deleted will
  % clear the MLR global
  views = MLR.views;
  viewCount = 0;
  for viewNum = 1:length(views)
    view = views{viewNum};
    if isview(view)
      viewCount = viewCount+1;
      delete(view.figure);
    end
  end
  if viewCount > 1
    disp(sprintf('(mrQuit) Closing %i open views',viewCount));
  end
  drawnow
  % save mrLastView
  if ~isempty(v)
    try
      disppercent(-inf,sprintf('(mrQuit) Saving %s/mrLastView',homeDir));
						% save the view in the current directory
      view = v;
      if getfield(whos('view'),'bytes')<2e9
        eval(sprintf('save %s view viewSettings -V6;',fullfile(homeDir,'mrLastView')));
      else
        mrWarnDlg('(mrQuit) Variable view is more than 2Gb, using option -v7.3 to save');
        eval(sprintf('save %s view viewSettings -v7.3;',fullfile(homeDir,'mrLastView')));
      end
      % save .mrDefaults in the home directory
      disppercent(inf);
    catch
      disppercent(inf);
      mrErrorDlg('(mrQuit) Could not save mrLastView.mat');
    end
  end
  disppercent(-inf,sprintf('(mrQuit) Saving %s',mrDefaultsFilename));
  saveMrDefaults;
  disppercent(inf);
else
  if ~isempty(v)
    closereq;
  end
end
clear global MLR

