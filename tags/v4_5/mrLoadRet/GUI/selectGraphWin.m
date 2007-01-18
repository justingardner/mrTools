function selectGraphWin(noClear)
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

h = MLR.graphFigure;
if (isempty(h) | h == 0) | ~ishandle(h)
  newGraphWin;
else
  set(0,'CurrentFigure',h);
end

% Clear the figure
if ieNotDefined('noClear') | (noClear == 0)
  clf
end  
return;
