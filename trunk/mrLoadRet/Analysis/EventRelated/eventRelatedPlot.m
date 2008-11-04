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
fignum = selectGraphWin;

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
  lhandle=legend(stimNames);
  set(lhandle,'Interpreter','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   make a global so that we can compute things that take   %%
%%   a long time only when user presses a button             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
gEventRelatedPlot.plotTSeriesHandle = [];
gEventRelatedPlot.computeErrorBarsHandle = [];
gEventRelatedPlot.roi = roi;
gEventRelatedPlot.time = time;
gEventRelatedPlot.cutoffr2 = cutoffr2;
gEventRelatedPlot.computingErrorBars = 0;
gEventRelatedPlot.loadingTimecourse = 0;

% set the delete function, so that we can delete the global we create
set(fignum,'DeleteFcn',@eventRelatedCloseWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there is an roi at this voxel
% then plot mean response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for roinum = 1:length(roi)
  subplot(2,2,4);
  ehdr = [];
  roin = 0;
  % first go for the quick and dirty way, which is
  % to load up the computed hemodynamic responses
  % and average them. 
  disppercent(-inf,'(eventRelatedPlot) Computing mean hdr');
  for voxnum = 1:size(roi{roinum}.scanCoords,2)
    disppercent(voxnum,size(roi{roinum}.scanCoords,2));
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
    disppercent(voxnum/size(roi{roinum}.scanCoords,2));
  end
  % plot the average of the ehdrs that beat the r2 cutoff
  if roin
    plotEhdr(time,shiftdim(mean(ehdr),1));
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
    lhandle = legend(stimNames);
    set(lhandle,'Interpreter','none');
  end
  % put up button whose call back will be to compute the error bars
  figpos = get(fignum,'position');
  gEventRelatedPlot.computeErrorBarsHandle = uicontrol('Parent',fignum,'Style','pushbutton','Callback',@eventRelatedPlotComputeErrorBars,'String','Compute error bars','Position',[figpos(3)/2+figpos(3)/20 figpos(4)/24 figpos(3)/2-figpos(3)/8 figpos(4)/14]);
  disppercent(inf);
end

drawnow;

% now put up a button that will call gEventRelatedPlotTSeries below
% when it is clicked but only if this is a long scan
if viewGet(view,'nFrames') > 500
  figpos = get(fignum,'position');
  gEventRelatedPlot.plotTSeriesHandle = uicontrol('Parent',fignum,'Style','pushbutton','Callback',@eventRelatedPlotTSeries,'String','Plot the time series','Position',[figpos(3)/20 5*figpos(4)/8 9*figpos(3)/10 figpos(4)/4]);
else
  eventRelatedPlotTSeries;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to plot the time series for the voxel   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotTSeries(varargin)

global gEventRelatedPlot
if gEventRelatedPlot.loadingTimecourse
  disp(sprintf('(eventRelatedPlot) Still loading timecourse. Please wait.'));
  return
end
gEventRelatedPlot.loadingTimecourse = 1;

disppercent(-inf,'(eventRelatedPlot) Plotting time series');
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
lhandle = legend(legendHandle,legendStr);
set(lhandle,'Interpreter','none');
% get distribution of ISI
%diff(sort(cell2mat(d.stimvol)));

if ~isempty(gEventRelatedPlot.plotTSeriesHandle)
  set(gEventRelatedPlot.plotTSeriesHandle,'Visible','off');
end

disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   eventRelatedPlotComputeErrorBars   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotComputeErrorBars(varargin)

global gEventRelatedPlot
if gEventRelatedPlot.computingErrorBars
  disp(sprintf('(eventRelatedPlot) Still computing error bars. Please wait.'));
  return
end

disp(sprintf('(eventRelatedPlot) Computing error bars over stimulus repetitions (i.e. averaging together all voxels that meet the r2 cutoff to form a single timecourse and then computing error bars using the inverse of the design covariance matrix)'));
gEventRelatedPlot.computingErrorBars = 1;
v = gEventRelatedPlot.v;
roi = gEventRelatedPlot.roi;
d = gEventRelatedPlot.d;
cutoffr2 = gEventRelatedPlot.cutoffr2;

subplot(2,2,4)

% get the time series
roi = loadROITSeries(v,roi{1});

n = 0;
for voxnum = 1:roi.n
  % get coordinates
  x = roi.scanCoords(1,voxnum);
  y = roi.scanCoords(2,voxnum);
  s = roi.scanCoords(3,voxnum);
  if d.r2(x,y,s) > cutoffr2
    n = n+1;
    meanTimecourse(n,:) = roi.tSeries(voxnum,:);
  end
end

yowsa = meanTimecourse;
if n == 0
  disp(sprintf('(eventRelatedPlot) No voxels met r2 > %0.3f',cutoffr2));
elseif n > 1
  meanTimecourse = mean(meanTimecourse);
end

% compute the event related analysis and the error bars
er = getr2timecourse(meanTimecourse,d.nhdr,d.hdrlen,d.scm,d.tr);

% plot them
plotEhdr(er.time,er.ehdr,er.ehdrste);

if ~isempty(gEventRelatedPlot.computeErrorBarsHandle)
  set(gEventRelatedPlot.computeErrorBarsHandle,'Visible','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end

% and display ehdr
for i = 1:size(ehdr,1)
  if ieNotDefined('ehdrste')
    h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   eventRelateCloseWindow   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedCloseWindow(varargin)

clear global gEventRelatedPlot;
