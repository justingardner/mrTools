function view = mrOpenWindow(viewType)
%
%  view = openWindow(viewType)
%
% djh, 6/2004

if ~exist('viewType','var')
    viewType = 'Volume';
end

mrGlobals;

view = newView(viewType);
% view is empty if it failed to initialize
if ~isempty(view)
    fig = mrLoadRetGUI('viewNum',view.viewNum);
    view = viewSet(view,'figure',fig);
else
    return
end
% set the location of the figure
figloc = mrGetFigLoc('mrLoadRetGUI');
if ~isempty(figloc)
    set(fig,'Position',figloc);
end

% Initialize the scan slider
nScans = viewGet(view,'nScans');
mlrGuiSet(view,'nScans',nScans);
mlrGuiSet(view,'scan',min(1,nScans));
% Initialize the slice slider
mlrGuiSet(view,'nSlices',0);

baseLoaded = 0;
if isfile('mrLastView.mat')
    mrLastView=load('mrLastView');
    % if the old one exists, then set up fields
    if isfield(mrLastView,'view')
        % open up base anatomy from last session
        if isfield(mrLastView.view,'baseVolumes')
            if length(mrLastView.view.baseVolumes) >= 1
                baseLoaded = 1;
                % Add it to the list of base volumes and select it
                view = viewSet(view,'newBase',mrLastView.view.baseVolumes);
            end
        end
        if baseLoaded && isfield(mrLastView,'viewSettings')
            % slice orientation from last run
            view = viewSet(view,'sliceOrientation',mrLastView.viewSettings.sliceOrientation);
            % rotate
            mlrGuiSet(view.viewNum,'rotate',mrLastView.viewSettings.rotate);
            % change group
            view = viewSet(view,'curGroup',mrLastView.viewSettings.curGroup);
            % change scan
            mlrGuiSet(view.viewNum,'scan',mrLastView.viewSettings.curScan);
            % change slice
            mlrGuiSet(view.viewNum,'slice',mrLastView.viewSettings.curSlice);
            % and refresh
            refreshMLRDisplay(view.viewNum);
        end
        % read ROIs into current view
        if isfield(mrLastView.view,'ROIs')
            for roinum = 1:length(mrLastView.view.ROIs)
                view = viewSet(view,'newROI',mrLastView.view.ROIs(roinum));
            end
        end

        % read analyses
        % if isfield(mrLastView.view,'analyses')
        %     for anum = 1:length(mrLastView.view.analyses)
        %         view = viewSet(view,'newAnalysis',mrLastView.view.analyses{anum});
        %     end
        % end
        % add here, to load more info...
    end

end

if ~baseLoaded
    % for when there is no mrLastView
    % open an anatomy, if there is one
    anatdir = dir('Anatomy/*.img');
    if ~isempty(anatdir)
        % load the first anatomy in the list
        view = loadAnat(view,anatdir(1).name);
        view = viewSet(view,'sliceOrientation','coronal');
        % rotate 270
        mlrGuiSet(view.viewNum,'rotate',270);
        % change group to last in list
        view = viewSet(view,'curGroup',viewGet(view,'numberOfGroups'));
        % and refresh
        refreshMLRDisplay(view.viewNum);
    end
end

