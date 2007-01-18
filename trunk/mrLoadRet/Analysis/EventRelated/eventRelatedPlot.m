function eventRelatedPlot(view,overlayNum,scan,x,y,z)
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
d = analysis.d;

% get the estimated hemodynamic responses
ehdr = gethdr(d,x,y,z);

% select the window to plot into
selectGraphWin;

% and display
time = d.tr/2:d.tr:(d.hdrlen*d.tr);
for i = 1:d.nhdr
  plot(time,ehdr(i,:),getcolor(i,'.-'));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,z,d.r2(x,y,z)));
xaxis(0,d.hdrlen*d.tr);
