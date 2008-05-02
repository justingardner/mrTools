function eventRelatedPlot(view,overlayNum,scan,x,y,s,roi)
% eventRelatedPlot.m
%
%       $Id$	
%      usage: eventRelatedPlot()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%


% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(eventRelatedPlot) Event related not for scan %i',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};
% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{end}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

% get cutoff value
cutoffr2 = viewGet(view,'overlayMin');

if isempty(d)
  disp('No analysis');
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the hemodynamic response for voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
[ehdr time ehdrste] = gethdr(d,x,y,s);
% display ehdr with out lines if we have a fit
% since we also need to plot fit
if isfield(d,'peak') & isfield(d.peak,'fit') & ~any(isnan(d.peak.amp(x,y,s,:)))
  plotEhdr(time,ehdr,ehdrste,'');
  for r = 1:d.nhdr
    d.peak.fit{x,y,s,r}.smoothX = 1:.1:d.hdrlen;
    fitTime = d.tr*(d.peak.fit{x,y,s,r}.smoothX-0.5);
    plot(fitTime+d.tr/2,d.peak.fit{x,y,s,r}.smoothFit,getcolor(r,'-'));
  end
else
  plotEhdr(time,ehdr,ehdrste);
end
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,s,analysis.overlays(1).data{scan}(x,y,s)));
xaxis(0,max(time));
% add peaks if they exist to the legend
if isfield(d,'stimNames')
  stimNames = d.stimNames;
  if isfield(d,'peak')
    for i = 1:d.nhdr
      stimNames{i} = sprintf('%s: %s=%0.2f',stimNames{i},d.peak.params.method,d.peak.amp(x,y,s,i));
    end
  end
  legend(stimNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there is an roi at this voxel
% then plot mean response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for roinum = 1:length(roi)
  subplot(2,2,4);
  ehdr = [];
  roin = 0;
  for voxnum = 1:size(roi{roinum}.scanCoords,2)
    x = roi{roinum}.scanCoords(1,voxnum);
    y = roi{roinum}.scanCoords(2,voxnum);
    s = roi{roinum}.scanCoords(3,voxnum);
    if d.r2(x,y,s) >= cutoffr2
      roin = roin+1;
      [ehdr(roin,:,:) time] = gethdr(d,x,y,s);
      % if there is a peak field, calculate average peak
      if isfield(d,'peak')
	for i = 1:d.nhdr
	  amp(i,roin) = d.peak.amp(x,y,s,i);
	end
      end
    end
  end
  % plot the average of the ehdrs that beat the r2 cutoff
  if roin
    plotEhdr(time,shiftdim(mean(ehdr),1),shiftdim(std(ehdr),1)/sqrt(size(roi{roinum}.scanCoords,2)));
  end
  title(sprintf('%s (n=%i/%i)',roi{roinum}.name,roin,size(roi{roinum}.scanCoords,2)),'Interpreter','none');
  % create a legend (only if peaks exist) to display mean amplitudes
  if isfield(d,'peak')
    for i = 1:d.nhdr
      % get the stimulus name
      if isfield(d,'stimNames')
	stimNames{i} = d.stimNames{i};
      else
	stimNames{i} = '';
      end
      % and now append the peak info
      stimNames{i} = sprintf('%s: median=%0.2f',stimNames{i},median(amp(i,:)));
    end
    legend(stimNames);
  end
end

drawnow;

% save info in a global so that we can plot the time course later, since this takes
% a long time to load the tseries
global gEventRelatedPlot;
% make a lighter view
v = newView;
v = viewSet(v,'curGroup',viewGet(view,'curGroup'));
v = viewSet(v,'curScan',viewGet(view,'curScan'));
gEventRelatedPlot.v = v;

gEventRelatedPlot.scan = scan;
gEventRelatedPlot.vox = [x y s];
gEventRelatedPlot.d = d;
gEventRelatedPlot.d.ehdr = [];
gEventRelatedPlot.d.ehdrste = [];
gEventRelatedPlot.d.r2 = [];
gEventRelatedPlot.h = [];

% now put up a button that will call gEventRelatedPlotTSeries below when it is clicked
% but only if this is a long scan
if viewGet(view,'nFrames') > 500
  figpos = get(fignum,'position');
  gEventRelatedPlot.h = uicontrol('Parent',fignum,'Style','pushbutton','Callback',@eventRelatedPlotTSeries,'String','Plot the time series','Position',[figpos(3)/20 5*figpos(4)/8 9*figpos(3)/10 figpos(4)/4]);
else
  eventRelatedPlotTSeries;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to plot the time series for the voxel   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotTSeries(varargin)

global gEventRelatedPlot
subplot(2,2,1:2)
tSeries = squeeze(loadTSeries(gEventRelatedPlot.v,gEventRelatedPlot.scan,gEventRelatedPlot.vox(3),[],gEventRelatedPlot.vox(1),gEventRelatedPlot.vox(2)));
legendHandle(1) = plot(tSeries,'k.-');
legendStr{1} = 'TSeries';
xlabel('Volume number');
ylabel('MRI signal');
% and the stimulus times
hold on
axis tight;
d = gEventRelatedPlot.d;
if isfield(d, 'stimvol')
  for i = 1:d.nhdr
    vlineHandle = vline(d.stimvol{i},getcolor(i));
    legendHandle(i+1) = vlineHandle(1);
    nStimvol(i) = length(d.stimvol{i});
    if isfield(d,'stimNames')
      legendStr{i+1} = sprintf('%s (n=%i)',d.stimNames{i},nStimvol(i));
    else
      legendStr{i+1} = sprintf('%i (n=%i)',i,nStimvol(i));
    end
  end
end
legend(legendHandle,legendStr);
% get distribution of ISI
%diff(sort(cell2mat(d.stimvol)));

if ~isempty(gEventRelatedPlot.h)
  set(gEventRelatedPlot.h,'Visible','off');
end
clear global gEventRelatedPlot;
%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end

% and display ehdr
for i = 1:size(ehdr,1)
  if nargin == 2
    h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');


