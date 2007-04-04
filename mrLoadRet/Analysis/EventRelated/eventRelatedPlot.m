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

% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');

% get the number of ROIs, we first will 
% look to see if the user clicked on a 
% voxel that is contained in an ROI.
roi = {};
nROIs = viewGet(view,'numberOfROIs');
for roinum = 1:nROIs
  roicoords = getRoiCoordinates(view,roinum,scan);
  % see if this is a matching roi
  if ismember([x y s],roicoords','rows')
    % get the roi
    roi{end+1} = viewGet(view,'roi',roinum);
    % change the coordinates to our coordinates
    roi{end}.coords = roicoords;
  end
end

if isempty(d)
  disp('No analysis');
  reutrn
end
 
% get the estimated hemodynamic responses
[ehdr time ehdrste] = gethdr(d,x,y,s);

% plot the timecourse
subplot(2,2,1:2)
tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));
plot(tSeries);
xlabel('Volume number');
ylabel('MRI signal');
% and the stimulus times
hold on
axis tight;
for i = 1:d.nhdr
  vline(d.stimvol{i},getcolor(i));
end
subplot(2,2,3);
% and display ehdr
plotEhdr(time,ehdr,ehdrste);
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,s,analysis.overlays(1).data{scan}(x,y,s)));
xaxis(0,d.hdrlen*d.tr);
if isfield(d,'stimNames')
  legend(d.stimNames);
end

% if there is a roi, compute its average hemodynamic response
for roinum = 1:length(roi)
  subplot(2,2,4);
  ehdr = [];
  for voxnum = 1:size(roi{roinum}.coords,2)
    [ehdr(voxnum,:,:) time] = gethdr(d,roi{roinum}.coords(1:3,voxnum));
  end
  plotEhdr(time,squeeze(mean(ehdr)),squeeze(std(ehdr))/sqrt(size(roi{roinum}.coords,2)));
  title(sprintf('%s (n=%i)',roi{roinum}.name,size(roi{roinum}.coords,2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste)

% and display ehdr
for i = 1:size(ehdr,1)
  if nargin == 2
    h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,'-')),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,'-')),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');


