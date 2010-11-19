function glmContrastPlot(thisView,overlayNum,scan,x,y,s,roi)
% glmContrastPlot.m
%
%        $Id$
%      usage: glmContrastPlot()
%         by: farshad moradi, modified by julien besle
%       date: 09/14/07, 12/02/2010
%    purpose: 

lineWidth = 2;
fontSize = 15; 

% check arguments
if ~any(nargin == [1:7])
  help glmContrastPlot
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data
% get the analysis structure
analysis = viewGet(thisView,'analysis');
if ~ismember(analysis.type,{'glmAnalStats','glmAnal','erAnal','deconvAnal'})
   disp(['(glmContrastPlot) Wrong type of analysis (' analysis.type ')']);
   return;
end
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(glmContrastPlot) No GLM analysis for scan %i',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};
numberOfConditions = d.nhdr;
if isfield(analysis.params, 'contrasts') && ~isempty(analysis.params.contrasts)
   d.contrasts = analysis.params.contrasts;
   numberOfConditions = size(d.contrasts,1);
end
if isfield(analysis.params, 'fTests') && ~isempty(analysis.params.fTests)
   d.fTests = analysis.params.fTests;
end


% check to see if there is a regular event related analysis
erAnalyses = [];
for anum = 1:viewGet(thisView,'nAnalyses')
  if ismember(viewGet(thisView,'analysisType',anum),{'erAnal','deconvAnal'})
    erAnalyses = [erAnalyses anum];
  end
end
if ~isempty(erAnalyses)
  if length(erAnalyses)==1
    erAnalNum = erAnalyses;
  else
    erAnalNum = 1:length(erAnalyses);
    while length(erAnalNum)>1
      erAnalNames = viewGet(thisView,'analysisNames');
      erAnalNum = find(buttondlg('Choose a deconvolution analysis or press Cancel',erAnalNames(erAnalyses)));
    end
    erAnalNum = erAnalyses(erAnalNum);
  end
  if ~isempty(erAnalNum)
    % get the event related heaerDatas
    erData = viewGet(thisView,'d',scan,erAnalNum);
  end
end

% set roi coords
for iRoi = 1:length(roi)
  % get scan coordinates
  roi{iRoi}.scanCoords = getROICoordinates(thisView,roi{iRoi},scan);
end

% get cutoff value
cutoffr2 = viewGet(thisView,'overlayMin');
framePeriod = viewGet(thisView,'framePeriod',scan);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set graph constants 
% select the window to plot into
fignum = selectGraphWin;


%set default values for figure aspect
set(fignum,'DefaultLineLineWidth',lineWidth);
set(fignum,'DefaultAxesFontSize',fontSize);
%set the colors
colors = color2RGB;
colors = colors([7 5 6 8 4 3 2 1 10]); %remove white and re-orerData
for i_color = 1:length(colors)
   colorOrder(i_color,:) = color2RGB(colors{i_color});
end
if numberOfConditions>size(colorOrder,1)
   colorOrder = repmat(colorOrder,ceil(numberOfConditions/size(colorOrder,1)),1);
end
colorOrder = colorOrder(1:numberOfConditions,:);

      
set(fignum,'DefaultAxesColorOrder',colorOrder);
%for bars, need to set the colormap
set(fignum,'colormap',colorOrder);

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','glmContrastPlot');

%set plotting dimension
subplot2marginRatio = 4;
legend2marginRatio = 1.5;
numberHorizontalPlots = (length(roi)+1);
margin = 1/((subplot2marginRatio+1)*numberHorizontalPlots+legend2marginRatio+1+1);
legendWidth = legend2marginRatio*margin;
subPlotWidth = subplot2marginRatio*margin;
subPlotHeight = (1-4*margin)/2;
subplotBetaPosition = [margin margin*3+subPlotHeight subPlotWidth subPlotHeight ];                                                                                                                                                 
subplotEhdrPosition = [margin margin*2 subPlotWidth subPlotHeight];
for iRoi = 1:length(roi)
   subplotRoiBetaPosition(iRoi,:) = [margin+iRoi*(subPlotWidth+margin) margin*3+subPlotHeight subPlotWidth subPlotHeight ];
   subplotRoiEhdrPosition(iRoi,:) = [margin+iRoi*(subPlotWidth+margin) margin*2 subPlotWidth subPlotHeight];
end
subplotRoiEhdrButtonPosition  = [margin .01 subPlotWidth margin/2];
legendBetaPosition = [1-margin-legendWidth margin*3+subPlotHeight legendWidth subPlotHeight];
legendEhdrPosition = [1-margin-legendWidth margin*2 legendWidth subPlotHeight];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the (contrast) betas for voxel
subplot('position',subplotBetaPosition);
betas = shiftdim(d.ehdr(x,y,s,:,:), 3);
if (isfield(d,'contrasts') && ~isempty(d.contrasts))
   betas = d.contrasts*betas;
   beta_errors = [];
   name = 'Contrasts';
   for iContrast = 1:size(d.contrasts,1)
      names{iContrast} = num2str(d.contrasts(iContrast,:));
   end
else
   beta_errors = shiftdim(d.ehdrste(x,y,s,:,:), 3);
   name = 'Explanatory Variables';
   for i_beta = 1:size(betas,1)
      names{i_beta} = [num2str(i_beta) ': ' d.stimNames{i_beta}];
   end
   %uicontrol('Parent',fignum, 'String',names,'style','text','unit','normalized','Position',[.91 .05 .08 .9]);
end
plotEcontrast(betas,beta_errors,name);
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,s,analysis.overlays(1).data{scan}(x,y,s)));
% subplot('position',legendBetaPosition);
% h=plotEcontrast(zeros(size(betas)),[],'');
% set(gca,'visible','off');
% set(h,'visible','off');
% lhandle = legend(names);
lhandle = legend(names,'position',legendBetaPosition);
set(lhandle,'Interpreter','none','box','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the hemodynamic response for voxel
subplot('position',subplotEhdrPosition);

% reconstruct what the glm produces
if ~isempty(analysis.params.scanParams{scan}) && isfield(analysis.params.scanParams{scan},'segmentBegin')
  % first want to get stimvols for just the beginning of the 
  % stimulus so that we can run an event-related analysis
  % on the glm modeled timecourse
  %
  % so, construct an appropriate params variable
  params = analysis.params.scanParams{scan};
  params.segmentNum = params.segmentBegin;
  params = rmfield(params,'segmentEnd');
  % get the stimvols 
  dglm = loadScan(view,scan,viewGet(view,'curGroup'),0); 
  dglm = getStimvol(dglm,params);
  % make stimulus convolution matrix
  dglm = makescm(dglm);
  % check to make sure dimensions match
  if size(d.scm,2) ~= size(d.ehdr,4)
    % no match, give up
    [ehdr time ehdrste] = gethdr(d,x,y,s);
    ehdrTitle = 'Scaled HRF';
  else
    % now we want to create a time course that is the glm design
    % multiplied by the beta weights (this is the estimated time
    % course)
    estTimecourse = d.scm*squeeze(d.ehdr(x,y,s,:));
    % now deconvolve to get the stimulus triggered averages
    % of the estimated timecourse
    estResponses = pinv(dglm.scm)*(estTimecourse-mean(estTimecourse));
    % and reshape
    ehdr = reshape(estResponses,dglm.hdrlen,dglm.nhdr)';
    time = dglm.tr/2:dglm.tr:(dglm.hdrlen*dglm.tr);
    ehdrste = [];
    ehdrTitle = 'GLM modeled response';
  end
else
  % get the estimated hdr from the beta weight times
  % the hdr function
  [ehdr time ehdrste] = gethdr(d,x,y,s);
  ehdrTitle = 'Scaled HRF';
end

if isfield(d,'contrasts') && ~isempty(d.contrasts)
   ehdr = d.contrasts*ehdr;
   ehdrste = [];
end
% display ehdr with out lines if we have a fit
% since we also need to plot fit
if isfield(d,'peak') & isfield(d.peak,'fit') & ~any(isnan(d.peak.amp(x,y,s,:)))
  plotEhdr(time,ehdr,ehdrste,'');
  for r = 1:d.nhdr
    d.peak.fit{x,y,s,r}.smoothX = 1:.1:d.hdrlen;
    fitTime = d.tr*(d.peak.fit{x,y,s,r}.smoothX-0.5);
    plot(fitTime+d.tr/2,d.peak.fit{x,y,s,r}.smoothFit,getcolor(r,'-'));
  end
    xaxis(0,time(end)+framePeriod/2);
else
  % if there is deconvolution data, display that too
  if exist('erData','var') %&& (length(erData.stimvol) == length(d.stimvol)) Why would we need to have the same length 
    [deconvEhdr deconvTime deconvEhdrste] = gethdr(erData,x,y,s);

    % note that we need to subtract the mean of the glm ehdr            %
    % to match the event related and the glm data. This is              %
    % because the glm has the mean subtracted from the columns          %
    % i.e. the glm gives an estimate of the response *before*           %
    % mean subtraction
    %deconvEhdr = deconvEhdr-repmat(mean(ehdr,2),1,size(deconvEhdr,2)); %NOT SO SURE... JB
    %what about this: 
    %centered_ehdr = ehdr-repmat(mean(ehdr,2),1,size(ehdr,2)); %Actually, let's just do nothing. 
    plotEhdr(time,ehdr,ehdrste,'-',0);
    hold on;
    h_deconv = plotEhdr(deconvTime,deconvEhdr,deconvEhdrste,'');
    set(h_deconv,'visible','off');
    uicontrol('Parent',fignum,...
       'units','normalized',...
       'Style','pushbutton',...
       'Callback',{@makeVisible,h_deconv},...
       'String','Plot estimated HDR',...
       'Position',subplotRoiEhdrButtonPosition,...
       'userdata','invisible');
    xaxis(0,deconvTime(end)+framePeriod/2);
  else
    plotEhdr(time,ehdr,ehdrste);
    xaxis(0,time(end)+framePeriod/2);
  end
end
title(ehdrTitle);

% add peaks if they exist to the legend
if isfield(d,'peak')
 for i = 1:d.nhdr
   names{i} = sprintf('%s: %s=%0.2f',names{i},d.peak.params.method,d.peak.amp(x,y,s,i));
 end
end
% subplot('position',legendEhdrPosition);
% plot([]);
% set(gca,'visible','off');
% lhandle = legend(names);
lhandle = legend(names,'position',legendEhdrPosition);
set(lhandle,'Interpreter','none','box','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if there is an roi at this voxel
% then plot mean beta values and Ehdr for this ROI
volumeEhdr = reshape(d.ehdr,[numel(d.r2) size(d.ehdr,4) size(d.ehdr,5)]);
volumeEhdrste = reshape(d.ehdrste,[numel(d.r2) size(d.ehdrste,4) size(d.ehdrste,5)]);
for iRoi = 1:length(roi)
  ehdr = [];
  econt = [];
  
  %array version
  %get ROI betas (econt)
  voxelIndices = sub2ind(size(d.r2),roi{iRoi}.scanCoords(1,:),roi{iRoi}.scanCoords(2,:),roi{iRoi}.scanCoords(3,:));
  voxelIndices = voxelIndices(d.r2(voxelIndices)>cutoffr2);
  voxelIndices = voxelIndices(~isnan(volumeEhdr(voxelIndices,1,1)));
  roiEhdr = volumeEhdr(voxelIndices,:,:);
  roiEhdrste = volumeEhdrste(voxelIndices,:,:);
  if ~isempty(roiEhdr)
     %if there is a contrast, compute the contrast values
     if isfield(d,'contrasts') && ~isempty(d.contrasts)
        for iHdr = 1:size(roiEhdr,3)
           roiEhdr(:,:,iHdr) = (d.contrasts*roiEhdr(:,:,iHdr)')';
        end
     end
     %put voxels on last dimension
     roiEhdr = permute(roiEhdr,[2 3 1]);
     roiEhdrste = permute(roiEhdrste,[2 3 1]);
     %scale the canonical hrf with contrast/beta value
     
     %construct a data structure with the mean beta estimates and std across voxels
     %although, ideally we would like to re compute the analysis like in eventRelatedPlot (but later...)
     roiD = d;
     roiD.dim = [1 1 1 d.dim(4)];
     roiD.ehdr = permute(mean(roiEhdr,3),[3 4 5 1 2]);
     %we could either compute the std error across voxels
     %roiD.ehdrste = permute(std(roiEhdr,0,3),[3 4 5 1 2])/sqrt(size(roi{iRoi}.scanCoords,2));
     %or take the mean std error across voxels
     roiD.ehdrste = permute(mean(roiEhdrste,3),[3 4 5 1 2]);
     %but what we should do is recompute...
     
     %get the scaled hdr and std errors
     [ehdr time ehdrste] = gethdr(roiD,1,1,1);

     if isfield(d,'peak')
        amp = reshape(d.peak.amp,numel(d.r2),d.nhdr);
        amp = amp(voxelIndices,:)';
     end

  % plot the average of the ehdrs that beat the r2 cutoff
    subplot('position',subplotRoiEhdrPosition(iRoi,:));
    plotEhdr(time,ehdr,ehdrste);
    xaxis(0,time(end)+framePeriod/2);
    subplot('position',subplotRoiBetaPosition(iRoi,:));
    %std error across voxels
    %plotEcontrast(mean(roiEhdr,3),mean(std(roiEhdr,0,3)/sqrt(size(roi{iRoi}.scanCoords,2)),name);
    %mean of std error across voxels
    plotEcontrast(mean(roiEhdr,3),mean(roiEhdrste,3),name);
  end
      
%   %loop version
%   roin = 0;
%   for voxnum = 1:size(roi{iRoi}.scanCoords,2)
%     x = roi{iRoi}.scanCoords(1,voxnum);
%     y = roi{iRoi}.scanCoords(2,voxnum);
%     s = roi{iRoi}.scanCoords(3,voxnum);
%     if d.r2(x,y,s) >= cutoffr2
%       roin = roin+1;
%       [ehdr(roin,:,:) time] = gethdr(d,x,y,s);
%       econt(roin,:,:)  = shiftdim(d.ehdr(x,y,s,:,:), 3);
%       % if there is a peak field, calculate average peak
%       if isfield(d,'peak')
% 	for i = 1:d.nhdr
% 	  amp(i,roin) = d.peak.amp(x,y,s,i);
% 	end
%       end
%     end
%   end
%   % plot the average of the ehdrs that beat the r2 cutoff
%   if ~isempty(ehdr)
%      if isfield(d,'contrasts') && ~isempty(d.contrasts)
%         for i_sample = 1:size(ehdr,3)
%            ehdr(:,:,i_sample) = (d.contrasts*ehdr(:,:,i_sample)')';
%         end
%          econt = (d.contrasts*econt')';
%      end
% 
%     subplot('position',subplotRoiEhdrPosition(iRoi,:));
%     plotEhdr(time,shiftdim(mean(ehdr),1),shiftdim(std(ehdr),1)/sqrt(size(roi{iRoi}.scanCoords,2)));
%     xaxis(0,time(end)+framePeriod/2);
%     subplot('position',subplotRoiBetaPosition(iRoi,:));
%     plotEcontrast(shiftdim(mean(econt),1),shiftdim(std(econt),1)/sqrt(size(roi{iRoi}.scanCoords,2)),name);
%   end

   title(sprintf('ROI %s (n=%i/%i)',roi{iRoi}.name,size(roiEhdr,3),size(roi{iRoi}.scanCoords,2)),'Interpreter','none');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function to plot ehdr  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotEhdr(time,ehdr,ehdrste,lineSymbol,drawSymbols)

colorOrder = get(gca,'colorOrder');
% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end
if ~exist('drawSymbols','var'),drawSymbols = 1;,end

% and display ehdr
if ieNotDefined('ehdrste')
   h=plot(repmat(time,size(ehdr,1),1)',ehdr',lineSymbol);
else
   h=errorbar(repmat(time,size(ehdr,1),1)',ehdr',ehdrste',ehdrste',lineSymbol);
end
 
if drawSymbols
   for iEv = 1:size(ehdr,1)
      set(h(iEv),'Marker',getsymbol(iEv),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colorOrder(iEv,:));
   end
end

xlabel('Time (sec)');
ylabel('% Signal change');

pos = get(gca,'position');
scale = axis;
uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.6*pos(4) .02 .03],'string',num2str(scale(4)),'callback',{@changeScale,gca,'max'});
uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.4*pos(4) .02 .03 ],'string',num2str(scale(3)),'callback',{@changeScale,gca,'min'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% function to plot contrasts  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotEcontrast(econt,econtste, name)

cla;
colorOrder = get(gca,'colorOrder');

% display econt
if size(econt,1)==1
    econt = econt';
    econtste = econtste';
end

if size(econt,2)==1
    set(gca,'nextPlot','add');
    for iEv = 1:size(econt,1)
       h(iEv) = bar(iEv,econt(iEv),'faceColor',colorOrder(iEv,:),'edgecolor','none');
    end
    set(gca,'xTickLabel',{})
    set(gca,'xTick',[])
else
   h = bar(econt','grouped','edgecolor','none');
   set(gca,'xtick',1:size(econt,2))
   xlabel('EV components');
end


if ~ieNotDefined('econtste')
    hold on;
    for i=1:length(h)
        % location of the bar
        x = get(get(h(i),'Children'),'XData');
        % find the center of the bar
        x = (x(2,:)+x(3,:))/2;
        
        herr = errorbar(x, econt(i,:), econtste(i,:), 'k');
        temp = get(herr, 'Children');
        set(temp(1), 'visible', 'off');
    end
end
xlabel(name);
switch(name)
   case 'Contrasts'
      yLabelString = 'Contrast values';
   case 'Explanatory Variables'
      yLabelString = 'Beta Values';
end
ylabel(yLabelString);

pos = get(gca,'position');
scale = axis;
uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.6*pos(4) .02 .03],'string',num2str(scale(4)),'callback',{@changeScale,gca,'max'});
uicontrol('style','edit','units','normalized','position',[pos(1)+pos(3)+.01 pos(2)+.4*pos(4) .02 .03 ],'string',num2str(scale(3)),'callback',{@changeScale,gca,'min'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% function to make lineseries visible  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeVisible(handle,eventdata,h_axis)

if strcmp(get(handle,'userdata'),'invisible')
   set(h_axis,'visible','on');
   set(handle,'String','Hide estimated HDR');
   set(handle,'userdata','visible');
else
   set(h_axis,'visible','off');
   set(handle,'String','Show estimated HDR');
   set(handle,'userdata','invisible');
end


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


