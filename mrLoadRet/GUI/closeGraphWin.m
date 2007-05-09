function closeGraphWin
%
% closeGraphWin
%
% Closes the current graphWin.  Sets the global GRAPHWIN=0
%
% djh, 3/3/98

mrGlobals

h = MLR.graphFigure;
if h
  if (h == get(0,'CurrentFigure'))
    MLR.graphFigure = []; 
    % save the position
    mrSetFigLoc('graphFigure',get(h,'Position'));
  end
end 

delete(get(0,'CurrentFigure'));


