function eventRelatedPlot(view,overlayNum,scan,x,y,s,roi)
% eventRelatedPlot.m
%
%       $Id$ 
%      usage: eventRelatedPlot()
%         by: justin gardner, modified by julien besle
%       date: 10/20/06
%    purpose: 
%
%     modifications: plots TSeries for ROI


lineWidth = 3;
fontSize = 15; 

% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data
% get the analysis structure
analysis = viewGet(view,'analysis');
if ~isfield(analysis,'d') || (length(analysis.d) < scan) || isempty(analysis.d)
  disp(sprintf('(eventRelatedPlot) Event related not for scan %i',scan));
  return
end
d = analysis.d{scan};
if isempty(d)
  mrWarnDlg(sprintf('(eventRelatedPlot) Could not find d structure for scan %i. Has eventRelated been run for this scan?',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};

% make a lighter view
v = newView;
v = viewSet(v,'curGroup',viewGet(view,'curGroup'));
v = viewSet(v,'curScan',viewGet(view,'curScan'));
plotTSeriesData.v = v;
errorBarData.v = v;

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{end}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

% get cutoff value
cutoffr2 = viewGet(view,'overlayMin');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set graph constants 
% select the window to plot into
fignum = selectGraphWin;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');
%set default values for figure aspect
set(fignum,'DefaultLineLineWidth',lineWidth);
set(fignum,'DefaultAxesFontSize',fontSize);
%set the colors to digit colors
colors = color2RGB;
colors = colors([7 5 6 8 4 3 2 1 10]);
for i_color = 1:length(colors)
   colorOrder(i_color,:) = color2RGB(colors{i_color});
end
set(fignum,'DefaultAxesColorOrder',colorOrder);
%set plotting dimension
subplot2marginRatio = 4;
legend2marginRatio = 1.5;
numberHorizontalPlots = (length(roi)+1);
numberVerticalPlots = (length(roi)+2);
margin = 1/((subplot2marginRatio+1)*numberHorizontalPlots+legend2marginRatio+1+1);
legendWidth = legend2marginRatio*margin;
totalPlotWidth = 1-2*margin;
subPlotWidth = subplot2marginRatio*margin;
subPlotHeight = (1-(numberVerticalPlots+2)*margin)/numberVerticalPlots;
subplotEhdrPosition = [margin margin*2+(numberVerticalPlots-1)*(subPlotHeight+margin) subPlotWidth subPlotHeight ];
subplotTSeriesPosition = [margin margin+(length(roi))*(subPlotHeight+margin) totalPlotWidth subPlotHeight];
for iRoi = 1:length(roi)
   subplotRoiEhdrPosition(iRoi,:) = [margin+iRoi*(subPlotWidth+margin) margin*2+(length(roi)+1)*(subPlotHeight+margin) subPlotWidth subPlotHeight ];
   subplotRoiEhdrButtonPosition(iRoi,:)  = [margin+iRoi*(subPlotWidth+margin) margin*.5+(length(roi)+1)*(subPlotHeight+margin) subPlotWidth margin/2];
   subplotRoiTSeriesPosition(iRoi,:) = [margin margin+(length(roi)-iRoi)*(subPlotHeight+margin) totalPlotWidth subPlotHeight];
end
legendPosition = [1-margin-legendWidth margin*2+(numberVerticalPlots-1)*(subPlotHeight+margin) legendWidth subPlotHeight];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the hemodynamic response for voxel
subplot('position',subplotEhdrPosition);

[voxelEhdr time ehdrste] = gethdr(d,x,y,s);
% display ehdr with out lines if we have a fit
% since we also need to plot fit
if isfield(d,'peak') & isfield(d.peak,'fit') & ~any(isnan(d.peak.amp(x,y,s,:)))
  plotEhdr(time,voxelEhdr,ehdrste,'');
  for r = 1:d.nhdr
    d.peak.fit{x,y,s,r}.smoothX = 1:.1:d.hdrlen;
    fitTime = d.tr*(d.peak.fit{x,y,s,r}.smoothX-0.5);
    plot(fitTime+d.tr/2,d.peak.fit{x,y,s,r}.smoothFit,getcolor(r,'-'),'lineWidth',lineWidth);
  end
else
  plotEhdr(time,voxelEhdr,ehdrste);
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
  lhandle=legend(stimNames,'position',legendPosition);
  set(lhandle,'Interpreter','none','box','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the time series for voxel

% now put up a button that will call EventRelatedPlotTSeries below
% when it is clicked 
plotTSeriesData.scan = scan;
plotTSeriesData.vox = [x y s];
plotTSeriesData.d = d;
plotTSeriesData.plotTSeriesHandle = [];
plotTSeriesData.computeErrorBarsHandle = [];
plotTSeriesData.time = time;
plotTSeriesData.cutoffr2 = cutoffr2;
plotTSeriesData.computingErrorBars = 0;
plotTSeriesData.loadingTimecourse = 0;
plotTSeriesData.ehdr = voxelEhdr;

uicontrol('Parent',fignum,...
        'units','normalized',...
        'Style','pushbutton',...
        'Callback',{@eventRelatedPlotTSeries,plotTSeriesData,'voxel'},...
        'String','Plot the time series for voxel',...
        'Position',subplotTSeriesPosition);
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if there is an roi at this voxel
%                                                      then plot mean response an time series
errorBarData.scan = scan;
errorBarData.vox = [x y s];
errorBarData.d = d;
errorBarData.voxelEhdr = voxelEhdr;
%errorBarData.ehdrste = [];
errorBarData.computeErrorBarsHandle = [];
errorBarData.time = time;
errorBarData.cutoffr2 = cutoffr2;
errorBarData.computingErrorBars = 0;

for iRoi = 1:length(roi)
   hSubplot = subplot('position',subplotRoiEhdrPosition(iRoi,:));
   roiEhdr = [];
  roin = 0;
  % first go for the quick and dirty way, which is
  % to load up the computed hemodynamic responses
  % and average them. 
  for voxnum = 1:size(roi{iRoi}.scanCoords,2)
    x = roi{iRoi}.scanCoords(1,voxnum);
    y = roi{iRoi}.scanCoords(2,voxnum);
    s = roi{iRoi}.scanCoords(3,voxnum);
    if d.r2(x,y,s) >= cutoffr2
      roin = roin+1;
      [roiEhdr(roin,:,:) time] = gethdr(d,x,y,s);
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
    plotEhdr(time,shiftdim(mean(roiEhdr),1));
  end
  title(sprintf('ROI %s (n=%i/%i)',roi{iRoi}.name,roin,size(roi{iRoi}.scanCoords,2)),'Interpreter','none');
  errorBarData.ehdr = squeeze(mean(roiEhdr));
  errorBarData.roi = roi{iRoi};
  errorBarData.hSubplot = hSubplot;
   uicontrol('Parent',fignum,...
      'units','normalized',...
     'Style','pushbutton',...
     'Callback',{@eventRelatedPlotComputeErrorBars,errorBarData},...
     'String','Compute error bars',...
     'Position',subplotRoiEhdrButtonPosition(iRoi,:));
  
  % create a legend (only if peaks exist) to display mean amplitudes
  % put up button whose call back will be to compute the error bars
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
    lhandle = legend(stimNames,'position',legendPosition);
    set(lhandle,'Interpreter','none');
  end
   
  %put a button to plot the time series of this roi
   plotTSeriesData.roi = roi{iRoi};
   plotTSeriesData.ehdr = squeeze(mean(roiEhdr,1));
   uicontrol('Parent',fignum,...
      'units','normalized',...
     'Style','pushbutton',...
     'Callback',{@eventRelatedPlotTSeries,plotTSeriesData,'roi'},...
     'String',['Plot the time series for ROI ' roi{iRoi}.name],...
     'Position',subplotRoiTSeriesPosition(iRoi,:));
end



drawnow;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function to plot the time series for the voxel and rois   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotTSeries(handle,eventData,plotTSeriesData,tSeriesType)

if plotTSeriesData.loadingTimecourse
  disp(sprintf('(eventRelatedPlot) Still loading timecourse. Please wait.'));
  return
end
plotTSeriesData.loadingTimecourse = 1;
d = plotTSeriesData.d;
scm = plotTSeriesData.d.scm;
position = get(handle,'position');

disppercent(-inf,'(eventRelatedPlot) Plotting time series');

subplot('position',position);
hold on
switch(tSeriesType)
   case 'voxel'
      tSeries = squeeze(loadTSeries(plotTSeriesData.v,plotTSeriesData.scan,plotTSeriesData.vox(3),[],plotTSeriesData.vox(1),plotTSeriesData.vox(2)));
      title('Voxel Time Series');
      modelTSeries = scm*reshape(plotTSeriesData.ehdr',numel(plotTSeriesData.ehdr),1)/100+1;
   case 'roi'
      fprintf(1,'\n');
      roi = loadROITSeries(plotTSeriesData.v,plotTSeriesData.roi);
      tSeries = mean(roi.tSeries,1);
      title(['ROI' plotTSeriesData.roi.name ' Time Series']);
      modelTSeries = scm*reshape(plotTSeriesData.ehdr',numel(plotTSeriesData.ehdr),1)/100+1;
end
junkFrames = viewGet(plotTSeriesData.v, 'junkFrames', plotTSeriesData.scan);
nFrames = viewGet(plotTSeriesData.v,'nFrames',plotTSeriesData.scan);
tSeries = tSeries(junkFrames+1:junkFrames+nFrames);
legendHandle(1) = plot(tSeries,'k.-');
legendHandle(2) = plot(modelTSeries,'--','Color',[.3 .3 .3]);

legendStr{1} = 'Actual TSeries';
legendStr{2} = 'Model TSeries';
xlabel('Volume number');
ylabel('MRI signal');
% and the stimulus times

axis tight;
if isfield(d, 'stimvol')
  colorOrder = get(gca,'colorOrder');
  scale = axis; 
  for iStim = 1:d.nhdr
    vlineHandle = plot([d.stimvol{iStim};d.stimvol{iStim}],repmat(scale(3:4)',1,size(d.stimvol{iStim},2)),'color',colorOrder(iStim,:));
    legendHandle(iStim+2) = vlineHandle(1);
    if isfield(d,'stimNames') && (length(d.stimNames) >= iStim)
      legendStr{iStim+2} = sprintf('%s (n=%i)',d.stimNames{iStim},length(d.stimvol{iStim}));
    else
      legendStr{iStim+2} = sprintf('%iStim (n=%i)',iStim,length(d.stimvol{iStim}));
    end
  end
end
lhandle = legend(legendHandle,legendStr);
set(lhandle,'Interpreter','none');

delete(handle);

disppercent(inf);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   eventRelatedPlotComputeErrorBars   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotComputeErrorBars(handle,eventData,errorBarData)

if errorBarData.computingErrorBars
  disp(sprintf('(eventRelatedPlot) Still computing error bars. Please wait.'));
  return
end

axes(errorBarData.hSubplot);

disp(sprintf('(eventRelatedPlot) Computing error bars over stimulus repetitions (i.e. averaging together all voxels that meet the r2 cutoff to form a single timecourse and then computing error bars using the inverse of the design covariance matrix)'));
errorBarData.computingErrorBars = 1;
v = errorBarData.v;
roi = errorBarData.roi;
d = errorBarData.d;
cutoffr2 = errorBarData.cutoffr2;

% get the time series
roi = loadROITSeries(v,roi);

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

if n == 0
  disp(sprintf('(eventRelatedPlot) No voxels met r2 > %0.3f',cutoffr2));
elseif n > 1
  meanTimecourse = mean(meanTimecourse);
end

% compute the event related analysis and the error bars
er = getr2timecourse(meanTimecourse,d.nhdr,d.hdrlen,d.scm,d.tr);

% plot them
titleString=get(get(gca,'title'),'string');
plotEhdr(er.time,er.ehdr,er.ehdrste);

title(titleString,'interpreter','none');

delete(handle)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function to plot ehdr  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end

% and display ehdr
  if ieNotDefined('ehdrste')
    h = plot(repmat(time,size(ehdr,1),1)',ehdr');
    minScale = min(ehdr(:));
    maxScale = max(ehdr(:));
  else
    h =errorbar(repmat(time,size(ehdr,1),1)',ehdr',ehdrste');
    minScale = min(ehdr(:)-ehdrste(:)/2);
    maxScale = max(ehdr(:)+ehdrste(:)/2);
  end

% for i = 1:size(ehdr,1)
%   if ieNotDefined('ehdrste')
%     h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
%   else
%     h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
%   end
%   set(h,'MarkerFaceColor',getcolor(i));
%   hold on
% end
xlabel('Time (sec)');
ylabel('% Signal change');
set(gca,'box','off')
scale(1:2) = [time(1)-1 time(end)+1];
scale(3) = floor(minScale);
scale(4) = ceil(maxScale);
axis(scale);
pos = get(gca,'position');

uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.6*pos(4) .02 .03],'string',num2str(scale(4)),'callback',{@changeScale,gca,'max'});
uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.4*pos(4) .02 .03 ],'string',num2str(scale(3)),'callback',{@changeScale,gca,'min'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%         changeScale        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeScale(handle,eventData,axisHandle,whichScale)
   
axes(axisHandle);
scale = axis;
switch(whichScale)
   case 'min'
      scale(3) = str2num(get(handle,'String'));
   case 'max'
      scale(4) = str2num(get(handle,'String'));
end
axis(scale);


