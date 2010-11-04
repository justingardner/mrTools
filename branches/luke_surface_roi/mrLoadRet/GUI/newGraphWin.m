function newGraphWin
%
% newGraphWin
%
% Opens a new window, sets the global MLR.graphFigure to be the handle to
% the new window. Sets closeRequestFcn to clean up properly (by calling
% closeGraphWin) when the window is closed.
%
% djh, 3/3/98
% djh, 9/2005 updated to MLR 4.0

mrGlobals
MLR.graphFigure = figure;
set(gcf,'CloseRequestFcn','closeGraphWin');
selectGraphWin;
figloc = mrGetFigLoc('graphFigure');
if ~isempty(figloc)
    set(MLR.graphFigure,'Position',figloc);
end