function eventRelatedPlot(view,overlayNum,scan,x,y,s)
% eventRelatedPlot.m
%
%      usage: eventRelatedPlot()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%


% check arguments
if ~any(nargin == [1:6])
  help eventRelatedPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};

% get the estimated hemodynamic responses
[ehdr time] = gethdr(d,x,y,s);

% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'MenuBar','none');
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');

% and display
for i = 1:d.nhdr
  h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,'-')),'MarkerSize',8);
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,s,analysis.overlays.data{scan}(x,y,s)));
xaxis(0,d.hdrlen*d.tr);
