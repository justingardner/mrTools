function h = selectGraphWin(noClear)
%
% selectGraphWin
%
% Checks the global MLR.graphFigure.  If nonzero, selects it. Otherwise,
% makes a new one.
%
% djh, 3/3/98
% djh, 9/2005 updated to MLR 4.0
% jlg 9/2006 added option not to clear window

mrGlobals;

% if there is no mrLoadRet running, just return a figure
if isempty(MLR.views)
  h = figure;
  return
end

% get figure from MLR variable
h = MLR.graphFigure;
if (isempty(h) | h == 0) | ~ishandle(h)
  newGraphWin;
  h = MLR.graphFigure;
else
  set(0,'CurrentFigure',h);
  figure(h);
end

% Clear the figure
if ieNotDefined('noClear') | (noClear == 0)
  clf
end  

return;
