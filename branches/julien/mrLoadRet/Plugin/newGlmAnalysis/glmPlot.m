function glmPlot(thisView,overlayNum,scanNum,x,y,s,roi)
% glmPlot.m
%
%        $Id$
%      usage: glmPlot() is an interrogator function
%         by: julien besle, modified from eventRelatedPlot and glmPlot
%       date: 09/14/07, 12/02/2010
%    purpose: plot GLM beta weights, contrasts, estimated HDR and time-series from GLM analysis


% check arguments
if ~any(nargin == [1:7])
  help glmPlot
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data
% get the analysis structure
analysisType = viewGet(thisView,'analysisType');
analysisParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));


if ~ismember(analysisType,{'glmAnalStats','glmAnal','glmcAnal','erAnal','deconvAnal'})
   disp(['(glmPlot) Wrong type of analysis (' analysisType ')']);
   return;
end
glmData = viewGet(thisView,'d');
if isempty(glmData)
  disp(sprintf('(glmPlot) No GLM analysis for scanNum %i',scanNum));
  return
end
r2data = viewGet(thisView,'overlaydata',scanNum,1);
r2clip = viewGet(thisView,'overlayClip',1);
numberEVs = glmData.nhdr;
framePeriod = viewGet(thisView,'framePeriod',scanNum);      

    
if isfield(glmData, 'contrasts') && ~isempty(glmData.contrasts) 
  numberContrasts = size(glmData.contrasts,1);
  if isfield(glmData,'EVnames')    
    contrastNames = makeContrastNames(glmData.contrasts,glmData.EVnames);
  else
    for iContrast = 1:numberContrasts
      contrastNames{iContrast} = num2str(glmData.contrasts(iContrast,:));
    end
  end
else
  numberContrasts=0;
end

if isfield(glmData,'EVnames')    
  EVnames = glmData.EVnames;
else
  for i_beta = 1:length(glmData.stimNames)
    EVnames{i_beta} = [num2str(i_beta) ': ' glmData.stimNames{i_beta}];
  end
end

if ismember(analysisType,{'glmAnalStats','glmAnal','glmcAnal'})
  plotBetaWeights = 1;
  plotDeconvolution=0;

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
        erAnalNum = find(buttondlg('Choose a deconvolution analysis or press OK if none is required',erAnalNames(erAnalyses)));
      end
      if all(~size(erAnalNum))
        return
      end
      erAnalNum = erAnalyses(erAnalNum);
    end
    if ~isempty(erAnalNum)
      plotDeconvolution=1;
      % get the event related data
      deconvData = viewGet(thisView,'d',scanNum,erAnalNum);
    end
    %
  end
else
  plotDeconvolution=0;
  plotBetaWeights = 0;
end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Set graph constants 

% select the window to plot into
fignum = selectGraphWin;
initializeFigure(fignum,max(numberEVs,numberContrasts))
set(fignum,'Name','glmPlot');

%set plotting dimension
maxNumberSte = 3;
subplotXgrid = [1.1 ones(1,length(roi)) .1 .4];
subplotYgrid = [.8*plotBetaWeights .6*logical(numberContrasts)*plotBetaWeights 1 logical(numberContrasts) .1 .1];
xMargin = .05;
yMargin = .01;
for iPlot = 1:length(roi)+1
  betaSubplotPosition(iPlot,:) = getSubplotPosition(iPlot,1,subplotXgrid,subplotYgrid,xMargin,yMargin);
  contrastSubplotPosition(iPlot,:) = getSubplotPosition(iPlot,2,subplotXgrid,subplotYgrid,xMargin,yMargin);
  ehdrSubplotPosition(iPlot,:) = getSubplotPosition(iPlot,3,subplotXgrid,subplotYgrid,xMargin,yMargin);
  contrastEhdrSubplotPosition(iPlot,:) = getSubplotPosition(iPlot,4,subplotXgrid,subplotYgrid,xMargin,yMargin);
  stePopupPosition(iPlot,:)  = getSubplotPosition(iPlot,5,subplotXgrid,subplotYgrid,xMargin,yMargin);
  tSeriesButtonPosition(iPlot,:) = getSubplotPosition(iPlot,6,subplotXgrid,subplotYgrid,xMargin,yMargin);
end
deconvButtonPosition  = getSubplotPosition(3+length(roi),5,subplotXgrid,subplotYgrid,xMargin,yMargin);
ehdrButtonPosition  = getSubplotPosition(3+length(roi),6,subplotXgrid,subplotYgrid,xMargin,yMargin);
legendBetaPosition = getSubplotPosition(3+length(roi),1,subplotXgrid,subplotYgrid,xMargin,yMargin);
legendContrastsPosition = getSubplotPosition(3+length(roi),2,subplotXgrid,subplotYgrid,xMargin,yMargin);
legendEhdrPosition = getSubplotPosition(3+length(roi),3,subplotXgrid,subplotYgrid,xMargin,yMargin);
legendHdrContrastPosition = getSubplotPosition(3+length(roi),4,subplotXgrid,subplotYgrid,xMargin,yMargin);


%&&&&&&&&&&&&&&&&&&&&&& LOOP on voxel + ROIs (if any) %&&&&&&&&&&&&&&&&&&&&&&

if ~isempty(roi)
  volumeBetas = reshape(glmData.ehdr,[numel(r2data) size(glmData.ehdr,4) size(glmData.ehdr,5)]);
%   volumeBetaSte = reshape(glmData.ehdrste,[numel(r2data) size(glmData.ehdrste,4) size(glmData.ehdrste,5)]);
%   if numberContrasts
%     volumeRSS = glmData.rss;
%   end
%   if plotDeconvolution
%     volumeDeconv = reshape(deconvData.ehdr,[numel(r2data) size(deconvData.ehdr,4) size(deconvData.ehdr,5)]);
%   end
end

hEhdr = [];
hDeconv = [];
for iPlot = 1:length(roi)+1
  hEhdrSte = zeros(numberEVs+numberContrasts,plotBetaWeights+1,maxNumberSte);
  if iPlot==1 %this is the voxel data
    titleString{1}=sprintf('Voxel (%i,%i,%i)',x,y,s);
    titleString{2}=sprintf('r2=%0.3f',r2data(x,y,s));
    volumeIndices = [x y s];
    e = getEstimates(glmData,analysisParams,volumeIndices);
    buttonString{1} = 'estimate std error';
    
  else  %this is an ROI
    
    hWait = uicontrol('parent',fignum,'style','text','unit','normalized',...
        'string','Computing Standard Errors for ROI. Please wait...','position',ehdrSubplotPosition(iPlot,:),'backgroundcolor',get(fignum,'color'));
    drawnow;
    roiNum = iPlot-1;
    % get roi scan coords
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum},scanNum);
    %get ROI estimates 
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    roiIndices = (r2data(volumeIndices)>r2clip(1)) & (r2data(volumeIndices)<r2clip(2));% & (~isnan(volumeBetas(volumeIndices,1,1)))';
    volumeIndices = volumeIndices(roiIndices);
    [e,volumeIndices] = getEstimates(glmData,analysisParams,volumeIndices);
    nVoxels = length(volumeIndices);
    nTotalVoxels = length(roiIndices);
    
    if nVoxels
      f=e;
      titleString{1}=sprintf('ROI %s (n=%i/%i)',roi{roiNum}.name,nVoxels,nTotalVoxels);
      titleString{2}=sprintf('%f<r2<%f',r2clip(1),r2clip(2));
      %compute mean across voxels
      e.betas = mean(f.betas,3);
      e.contrastBetas = mean(f.contrastBetas,3);
      e.hdr = mean(f.hdr,3);
      e.contrastHdr = mean(f.contrastHdr,3);
      
      %there are several possible ways of computing the ROI std error:
      % 1) as the mean of the std error across voxels (=treating the standard error as a measure rather than an estimate)
      buttonString{1} = 'average (across voxels) of the voxel-wise estimate standard errors';
      e.betaSte = mean(f.betaSte,3);
      e.contrastBetaSte = mean(f.contrastBetaSte,3);
      e.hdrSte = mean(f.hdrSte,3);
      e.contrastHdrSte = mean(f.contrastHdrSte,3);
      
      % 2) as the std error of the estimate across voxels (=ignoring the intra-voxels variability = treating the estimate as a measure)
      buttonString{2} = 'standard error (across voxels) of the voxel-wise estimates';
      e.betaSte(:,:,2) = std(f.betas,0,3)/sqrt(nVoxels);
      e.contrastBetaSte(:,:,2) = std(f.contrastBetas,0,3)/sqrt(nVoxels);
      e.hdrSte(:,:,2) = std(f.hdr,0,3)/sqrt(nVoxels);
      e.contrastHdrSte(:,:,2) = std(f.contrastHdr,0,3)/sqrt(nVoxels);
      
      % 3) as the std error of a sum (mean) of estimates: (=treating the ROI beta as a random variable that is a weigthed sum of estimates)
      buttonString{3} = 'ROI estimate standard error (assuming no spatial correlation)';
      e.betaSte(:,:,3) = sqrt(mean(f.betaSte.^2,3));
      e.contrastBetaSte(:,:,3) = sqrt(mean(f.contrastBetaSte.^2,3));
      e.hdrSte(:,:,3) = sqrt(mean(f.hdrSte.^2,3));
      e.contrastHdrSte(:,:,3) = sqrt(mean(f.contrastHdrSte.^2,3));
      
      % 4) as the std error of an estimate from the mean time-series (by rerunning the glm analysis)
      %buttonString{4} = 'ROI estimate standard error from ROI time-series';
      %later...
      
      
            % % %       % create a legend (only if peaks exist) to display mean amplitudes
            % % %       if isfield(glmData,'peak')
            % % %         for i = 1:glmData.nhdr
            % % %           % get the stimulus name
            % % %           if isfield(glmData,'stimNames')
            % % %             stimNames{i} = glmData.stimNames{i};
            % % %           else
            % % %             stimNames{i} = '';
            % % %           end
            % % %           % and now append the peak info
            % % %           stimNames{i} = sprintf('%s: median=%0.2f',stimNames{i},median(amp(i,:)));
            % % %         end
            % % %         legend(stimNames);
            % % %       end
            
%       if plotDeconvolution
%         deconvData.ehdr = permute(mean(volumeDeconv(volumeIndices,:,:),1),[4 5 1 2 3]);
%       end
    end
    delete(hWait);
  end
  
  %mean confidence intervals across voxels 
  if isfield(glmData,'ehdrBootstrapCIs')
    e.betaSte(:,:,end+1) = mean(e.bootstrapBetaCIs,3);
    if iPlot==1
      buttonString{end+1} = ['Bootstrap Confidence Intervals (alpha = ' num2str(analysisParams.alphaConfidenceIntervals) ')'];
    else
      buttonString{end+1} = ['average (across voxels) of voxel-wise Bootstrap Confidence Intervals (alpha = ' num2str(analysisParams.alphaConfidenceIntervals) ')'];
    end
    e.contrastBetaSte(:,:,end+1) = mean(e.bootstrapContrastCIs,3);
    %no value for HDR
    e.hdrSte(:,:,end+1) = NaN(size(e.hdrSte(:,:,1)));
    e.contrastHdrSte(:,:,end+1) = NaN(size(e.contrastHdrSte(:,:,1)));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare axes
  if plotBetaWeights
    betaAxes = axes('parent',fignum,'outerposition',betaSubplotPosition(iPlot,:),'DrawMode','fast');
    hold on
    %hold(betaAxes);
    title(titleString,'Interpreter','none');
    if numberContrasts
      contrastAxes = axes('parent',fignum,'outerposition',contrastSubplotPosition(iPlot,:),'DrawMode','fast');
      hold on
      %hold(contrastAxes);
    end
  end
  ehdrAxes = axes('parent',fignum,'outerposition',ehdrSubplotPosition(iPlot,:),'DrawMode','fast');
  hold on
  %hold(ehdrAxes);
  %plot baseline
  plot(ehdrAxes,[0 (glmData.hdrlen+1)*framePeriod],[0 0],'--k','lineWidth',1);
  if numberContrasts
    hdrContrastAxes = axes('parent',fignum,'outerposition',contrastEhdrSubplotPosition(iPlot,:),'DrawMode','fast');
    hold on
    %hold(hdrContrastAxes);
    %plot baseline
    plot(hdrContrastAxes,[0 (glmData.hdrlen+1)*framePeriod],[0 0],'--k','lineWidth',1);
  end
  
      
  %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& PLOT DATA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  for iSte = 1:size(e.betaSte,3)
    if plotBetaWeights
      % plot the beta weights
      [h,hEhdrSte(1:numberEVs,1,iSte)] = plotBetas(betaAxes,e.betas,e.betaSte(:,:,iSte),iSte~=1);
      if iSte==1 && iPlot==1,hBeta=h;end;

      if numberContrasts
        % plot the contrast estimates
        [h,hEhdrSte(numberEVs+1:numberEVs+numberContrasts,1,iSte)] = plotBetas(contrastAxes,e.contrastBetas,e.contrastBetaSte(:,:,iSte),iSte~=1);
        if iSte==1 && iPlot==1,hContrastBeta=h;end;
      end
    end  
    % plot the hemodynamic response for voxel
    [h,hEhdrSte(1:numberEVs,plotBetaWeights+1,iSte)]=plotEhdr(ehdrAxes,e.time,e.hdr,e.hdrSte(:,:,iSte),[],[],iSte~=1);
    if iSte==1, hHdr = h; hEhdr = [hEhdr;h];end
    if numberContrasts
      [h,hEhdrSte(numberEVs+1:numberEVs+numberContrasts,plotBetaWeights+1,iSte)] = plotEhdr(hdrContrastAxes,e.time,e.contrastHdr, e.contrastHdrSte(:,:,iSte),'','',iSte~=1);
      if iSte==1, hContrastHdr =h; hEhdr = [hEhdr;h];end;
    end
    
    if iSte~=1
      set(hEhdrSte(:,:,iSte),'visible','off');
    end

  end
  uicontrol('Parent',fignum,...
     'units','normalized',...
     'Style','popupmenu',...
     'Callback',{@makeVisible,hEhdrSte},...
     'String', [buttonString {'No error bars'}],...
     'Position',stePopupPosition(iPlot,:));
                                                                 % % % % display ehdr with out lines if we have a fit
                                                                  % % % % since we also need to plot fit
                                                                  % % % if isfield(glmData,'peak') & isfield(glmData.peak,'fit') & ~any(isnan(glmData.peak.amp(x,y,s,:)))
                                                                  % % %   h = plotEhdr(e.time,e.hdr,e.hdrSte,'');
                                                                  % % %   for r = 1:glmData.nhdr
                                                                  % % %     glmData.peak.fit{x,y,s,r}.smoothX = 1:.1:glmData.hdrlen;
                                                                  % % %     fitTime = glmData.tr*(glmData.peak.fit{x,y,s,r}.smoothX-0.5);
                                                                  % % %     plot(fitTime+glmData.tr/2,glmData.peak.fit{x,y,s,r}.smoothFit,getcolor(r,'-'));
                                                                  % % %   end
                                                                  % % %     xaxis(0,e.time(end)+framePeriod/2);
                                                                  % % % end
                                                                  % % % % add peaks if they exist to the legend
                                                                  % % % if isfield(glmData,'peak')
                                                                  % % %  for i = 1:glmData.nhdr
                                                                  % % %    names{i} = sprintf('%s: %s=%0.2f',names{i},glmData.peak.params.method,glmData.peak.amp(x,y,s,i));
                                                                  % % %  end
                                                                  % % % end

  % if there is deconvolution data, display that too
  if plotDeconvolution 
    keyboard %I need to finish this
    if numberContrasts
      eDeconv = getEstimates(deconvData,analysisParams,volumeIndices);
      [deconvEhdr deconvTime deconvEhdrste deconvContrastHdr] = getEstimates(deconvData,x,y,s,glmData.contrasts);
    else
      [deconvEhdr deconvTime] = getEhdr(deconvData,x,y,s);
    end
    if any(any(deconvEhdr))
      set(hEhdr,'MarkerEdgeColor','none','MarkerFaceColor','none');
      hDeconv = [hDeconv; plotEhdr(ehdrAxes,deconvTime,deconvEhdr)];
      if numberContrasts
        hDeconv = [hDeconv; plotEhdr(hdrContrastAxes,deconvTime,deconvContrastHdr)];
      end
    end
  else
    deconvTime=0;
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finalize axes
  %plot baselines of histograms
  if plotBetaWeights
    plot(betaAxes,get(betaAxes,'Xlim'),[0 0],'--k','lineWidth',1);
      maxSte = max(e.betaSte,[],3);
      makeScaleEditButton(fignum,betaAxes,...
      [min(min((e.betas-maxSte))),max(max((e.betas+maxSte)))]);
    if iPlot==1
      ylabel(betaAxes,{'Beta' 'Estimates'});
      lhandle = legend(hBeta,EVnames,'position',legendBetaPosition);
      set(lhandle,'Interpreter','none','box','off');
    end
    if numberContrasts
      %plot baseline
      plot(contrastAxes,get(contrastAxes,'Xlim'),[0 0],'--k','lineWidth',1);
      maxSte = max(e.contrastBetaSte,[],3);
      if isnan(maxSte)
        maxSte = 0;
      end
      makeScaleEditButton(fignum,contrastAxes,...
        [min(min(e.contrastBetas-maxSte)),max(max(e.contrastBetas+maxSte))]);
      if iPlot==1
        ylabel(contrastAxes,{'Contrast' 'Estimates'});
        lhandle = legend(hContrastBeta,contrastNames,'position',legendContrastsPosition);
        set(lhandle,'Interpreter','none','box','off');
      end
    end
  end
  set(ehdrAxes,'xLim',[0,max(deconvTime(end),e.time(end))+framePeriod/2]);
  maxSte = abs(max(e.hdrSte,[],3));
  makeScaleEditButton(fignum,ehdrAxes,...
    [min(min((e.hdr-maxSte))),max(max((e.hdr+maxSte)))]);
  if iPlot==1
    ylabel(ehdrAxes,{'Scaled HRF','% Signal change'});
    lhandle = legend(hHdr,EVnames,'position',legendEhdrPosition);
    set(lhandle,'Interpreter','none','box','off');
  end
  if numberContrasts
    set(hdrContrastAxes,'xLim',[0,max(deconvTime(end),e.time(end))+framePeriod/2]);
    maxSte = max(e.contrastHdrSte,[],3);
    if isnan(maxSte)
      maxSte = 0;
    end
    makeScaleEditButton(fignum,hdrContrastAxes,...
      [min(min(e.contrastHdr-maxSte)),max(max(e.contrastHdr+maxSte))]);
    if iPlot==1
      ylabel(hdrContrastAxes,{'Scaled Contrast HRF','% Signal change'});
      lhandle = legend(hContrastHdr,contrastNames,'position',legendHdrContrastPosition);
      set(lhandle,'Interpreter','none','box','off');
    end
    xlabel(hdrContrastAxes,'Time (sec)');
  else
    xlabel(ehdrAxes,'Time (sec)');
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the time series
    %put a button to plot the time series of this roi
  if iPlot==1
    thisRoi = [x y s];
  else
    thisRoi = roi{roiNum};
  end
  uicontrol('Parent',fignum,...
    'units','normalized',...
   'Style','pushbutton',...
   'Callback',{@eventRelatedPlotTSeries, thisView, analysisParams, glmData, thisRoi},...
   'String',['Plot the time series for ' titleString{1}],...
   'Position',tSeriesButtonPosition(iPlot,:));


end


if plotDeconvolution && ~isempty(hDeconv)
  set(hDeconv,'visible','off');
  uicontrol('Parent',fignum,...
     'units','normalized',...
     'Style','pushbutton',...
     'Callback',{@makeVisible,hDeconv},...
     'String','Show deconvoluted HDR',...
     'Position',deconvButtonPosition);
  uicontrol('Parent',fignum,...
   'units','normalized',...
   'Style','pushbutton',...
   'Callback',{@makeVisible,hEhdr},...
   'String','Hide estimated HDR',...
   'Position',ehdrButtonPosition);
end

drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function to plot the time series for the voxel and rois   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventRelatedPlotTSeries(handle,eventData,thisView,analysisParams, d, roi)

fignum = selectGraphWin(0,'Make new');
initializeFigure(fignum,d.nhdr);
monitorPositions = getMonitorPositions;
figurePosition = get(fignum,'position');
[whichMonitor,figurePosition]=getMonitorNumber(figurePosition,monitorPositions);
screenSize = monitorPositions(whichMonitor,:); % find which monitor the figure is displayed in
figurePosition(3) = screenSize(3);
figurePosition(4) = screenSize(3)/5;
set(fignum,'position',figurePosition);
h=uicontrol('parent',fignum,'unit','normalized','style','text',...
  'string','Loading timecourse. Please wait...',...
  'position',[0.1 0.45 .8 .1],'backgroundcolor',get(fignum,'color'));
drawnow;
disppercent(-inf,'(eventRelatedPlot) Plotting time series');

if isnumeric(roi)
  tSeries = squeeze(loadTSeries(thisView,[],roi(3),[],roi(1),roi(2)));
  titleString = sprintf('Voxel %i,%i,%i Time Series',roi(1),roi(2),roi(3));
  ehdr = shiftdim(d.ehdr(roi(1),roi(2),roi(3),:,:),3);
else
  fprintf(1,'\n');
  roi = loadROITSeries(thisView,roi);
  tSeries = mean(roi.tSeries,1);
  titleString = ['ROI ' roi.name ' Time Series'];
  ehdr = permute(d.ehdr,[4 5 1 2 3]);
  ehdr = ehdr(:,:,sub2ind(d.dim(1:3),roi.scanCoords(1,:)',roi.scanCoords(2,:)',roi.scanCoords(3,:)'));
  ehdr = mean(ehdr,3);
end
%convert to percent signal change the way it's done in getGlmStatistics
tSeries = (tSeries - mean(tSeries))'/mean(tSeries)*100;
junkFrames = viewGet(thisView, 'junkFrames');
nFrames = viewGet(thisView,'nFrames');
tSeries = tSeries(junkFrames+1:junkFrames+nFrames);

if isfield(analysisParams,'scanParams') && isfield(analysisParams.scanParams{thisView.curScan},'acquisitionSubsample')...
    && ~isempty(analysisParams.scanParams{thisView.curScan}.acquisitionSubsample)
  acquisitionSubsample = analysisParams.scanParams{thisView.curScan}.acquisitionSubsample;
else
  acquisitionSubsample = 1;
end
if ~isfield(d,'estimationSupersampling')
  d.estimationSupersampling=1;
end
time = ((1:length(tSeries))+(acquisitionSubsample-.5)/d.estimationSupersampling)*d.tr;


delete(h);
set(fignum,'name',titleString)
tSeriesAxes = axes('parent',fignum,'outerposition',getSubplotPosition(1,1,[7 1],1,0,0));
hold on
%hold(tSeriesAxes);


%Plot the stimulus times
set(tSeriesAxes,'Ylim',[min(tSeries);max(tSeries)])
if ~isfield(d,'designSupersampling')
  d.designSupersampling=1;
end

colorOrder = get(tSeriesAxes,'colorOrder');
if isfield(d,'EVmatrix') && isfield(d,'EVnames')
  stimOnsets = d.EVmatrix;
  stimDurations = [];
  legendString = d.EVnames;
elseif isfield(d,'stimDurations') && isfield(d, 'stimvol')
  stimOnsets = d.stimvol;
  stimDurations = d.stimDurations;
  legendString = d.stimNames;
elseif isfield(d, 'stimvol')
  stimOnsets = d.stimvol;
  stimDurations = [];
  legendString = d.stimNames;
end
if isfield(d,'runTransitions')
  runTransitions = d.runTransitions;
else
  runTransitions = [];
end

[h,hTransitions] = plotStims(stimOnsets, stimDurations, d.tr/d.designSupersampling, colorOrder, tSeriesAxes, runTransitions);
legendString = legendString(h>0);
h = h(h>0);

%and the time-series
h(end+1) = plot(tSeriesAxes,time,tSeries,'k.-');
if size(d.scm,2)==numel(ehdr)
  %compute model time series
  modelTSeries = d.scm*reshape(ehdr',numel(ehdr),1);
  h(end+1) = plot(tSeriesAxes,time,modelTSeries,'--r');
end
legendString{end+1} = 'Actual TSeries';
legendString{end+1} = 'Model TSeries';
if ~isempty(hTransitions)
  h = [h hTransitions];
  legendString{end+1} = 'Run transitions';
end
ylabel('Percent Signal Change');
xlim([0 ceil(time(end)+1)]);
%legend
lhandle = legend(h,legendString,'position',getSubplotPosition(2,1,[7 1],1,0,.2));
set(lhandle,'Interpreter','none','box','off');

disppercent(inf);
%delete(handle); %don't delete the button to plot the time-series



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function to initialize figure%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeFigure(fignum,numberColors)

lineWidth = 2;
fontSize = 15; 

%set default values for figure aspect
set(fignum,'DefaultLineLineWidth',lineWidth);
set(fignum,'DefaultAxesFontSize',fontSize);
%set the colors
colors = color2RGB;
colors = colors([7 5 6 8 4 3 2 1]); %remove white and black and re-orerData
for i_color = 1:length(colors)
   colorOrder(i_color,:) = color2RGB(colors{i_color});
end
if numberColors>size(colorOrder,1)
   colorOrder(end+1:numberColors,:) = randomColors(numberColors-size(colorOrder,1));
end
colorOrder = colorOrder(1:numberColors,:);

      
set(fignum,'DefaultAxesColorOrder',colorOrder);
%for bars, need to set the colormap
set(fignum,'colormap',colorOrder);

% turn off menu/title etc.
set(fignum,'NumberTitle','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function to plot ehdr  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,hSte] = plotEhdr(hAxes,time,ehdr,ehdrSte,lineSymbol,drawSymbols, steOnly)

colorOrder = get(hAxes,'colorOrder');
% whether to plot the line inbetween points or not
if ieNotDefined('lineSymbol'),lineSymbol = '-';,end
if ieNotDefined('drawSymbols'),drawSymbols = 1;,end
if ieNotDefined('steOnly'),steOnly = 0;,end

% and display ehdr
if steOnly
  h=[];
else
  h=plot(hAxes,repmat(time,1,size(ehdr,2)),ehdr,lineSymbol);
  if drawSymbols
     for iEv = 1:size(ehdr,2)
        set(h(iEv),'Marker',getsymbol(iEv),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colorOrder(iEv,:));
     end
  end
end

if ~ieNotDefined('ehdrSte')
  hold on
  %if ~ishold(hAxes),hold(hAxes);end;
  hSte=errorbar(hAxes,repmat(time,1,size(ehdr,2)),ehdr,ehdrSte,ehdrSte,'lineStyle','none')';
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% function to plot contrasts  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h hSte]= plotBetas(hAxes,econt,econtste, steOnly)

colorOrder = get(hAxes,'colorOrder');
if ieNotDefined('steOnly'),steOnly = 0;,end

% display econt
if size(econt,1)==1
  econt = econt';
  econtste = econtste';
end

if size(econt,2)==1
  set(hAxes,'nextPlot','add');
  h=zeros(size(econt,1),1);
  for iEv = 1:size(econt,1)
     h(iEv) = bar(hAxes,iEv,econt(iEv),'faceColor',colorOrder(iEv,:),'edgecolor','none');
  end
  %delete baseline
  delete(get(h(iEv),'baseline'));
  set(hAxes,'xTickLabel',{})
  set(hAxes,'xTick',[])
else
  h = bar(hAxes,econt','grouped','edgecolor','none');
  set(hAxes,'xtick',1:size(econt,2))
  xlabel('EV components');
end
if steOnly
  set(h,'visible','off');
end

if ~ieNotDefined('econtste')
  hold on
  %if ~ishold(hAxes),hold(hAxes);end;
  hSte = zeros(size(h));
  for i=1:length(h)
    % location of the bar
    x = get(get(h(i),'Children'),'XData');
    % find the center of the bar
    x = (x(2,:)+x(3,:))/2;
    hSte(i) = errorbar(hAxes,x, econt(i,:), econtste(i,:), 'k','lineStyle','none');
%         temp = get(hSte(i), 'Children');
%         set(temp(1), 'visible', 'off');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% function to make scale edit boxes  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeScaleEditButton(fignum,axisHandle,minMaxData)

pos = get(axisHandle,'position');
if ieNotDefined('minMaxData')
  minMaxData = get(axisHandle,'YLim');
else
  minMaxData(1) = minMaxData(1)-.02*diff(minMaxData);
  minMaxData(2) = minMaxData(2)+.02*diff(minMaxData);
end
set(axisHandle,'YLim',minMaxData);
uicontrol(fignum,'style','edit','units','normalized',...
  'position',[pos(1)+pos(3) pos(2)+.6*pos(4) .03 .03],...
  'string',num2str(minMaxData(2)),'callback',{@changeScale,axisHandle,'max'});
uicontrol(fignum,'style','edit','units','normalized',...
  'position',[pos(1)+pos(3) pos(2)+.4*pos(4) .03 .03 ],...
  'string',num2str(minMaxData(1)),'callback',{@changeScale,axisHandle,'min'});
set(axisHandle,'YLimMode','manual');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% function to make lineseries visible  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeVisible(handle,eventdata,hAxes)

switch(get(handle,'style'))
  case 'pushbutton'
    string = get(handle,'String');
    if strcmp(get(hAxes,'visible'),'off')
       set(hAxes,'visible','on');
       set(handle,'String',['Hide' string(5:end)]);
       set(handle,'userdata','visible');
    else
       set(hAxes,'visible','off');
       set(handle,'String',['Show' string(5:end)]);
       set(handle,'userdata','invisible');
    end
    
  case 'popupmenu'
    set(hAxes,'visible','off')
    handleNum = get(handle,'value');
    if handleNum ~= length(get(handle,'string'));
      set(hAxes(:,:,handleNum),'visible','on');
    end
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


