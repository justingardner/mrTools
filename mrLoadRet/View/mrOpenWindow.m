function view = mrOpenWindow(viewType)
%
%  view = openWindow(viewType)
%
% djh, 6/2004

if ~exist('viewType','var')
    viewType = 'Volume';
end

view = newView(viewType);
% view is empty if it failed to initialize
if ~isempty(view)
    fig = mrLoadRetGUI('viewNum',view.viewNum);
    view = viewSet(view,'figure',fig);
	
	% Initialize the scan slider
	nScans = viewGet(view,'nScans');
    mlrGuiSet(view,'nScans',nScans);
	mlrGuiSet(view,'scan',min(1,nScans));
    % Initialize the slice slider
    mlrGuiSet(view,'nSlices',0);
end