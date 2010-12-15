% getGlmEVParamsGUI.m
%
%        $Id$
%      usage: params = getGlmEVParamsGUI(thisView,params,useDefault)
%         by: julien besle
%       date: 03/12/2010
%    purpose: return EV specific parameters for GLM analysis
%

function [scanParams, params] = getGlmScanParamsGUI(thisView,params,useDefault)

keepAsking = 1;
groupNum = viewGet(thisView,'groupNum',params.groupName);
nScans = viewGet(thisView,'nScans',groupNum);
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = 1:nScans;
end
if isfield(params,'scanParams') && length(params.scanParams)==nScans
   scanParams = params.scanParams;
else
   % make the output as long as the number of scans
   scanParams = cell(1,nScans);
end
if ~isfield(params, 'EVnames') || ~isequal(length(params.EVnames),params.numberEVs) 
  EVnames = {};
else
  EVnames = params.EVnames;
end

while keepAsking
  % check for stimfile, and if it is mgl/type then ask the
  % user which variable they want to do the anlysis on
  for iScan = params.scanNum
    
    %get the number of events after running the pre-processing function
    disp(sprintf('Getting number of conditions from scan %d', iScan)); 
    d = loadScan(thisView, iScan, groupNum, 0);
    d = getStimvol(d,scanParams{iScan});
    d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
    params.stimNames = d.stimNames;
    nStims = length(d.stimvol);
    disp(sprintf('%d conditions found', nStims));
    
    if ~params.numberEVs
      params.numberEVs = nStims;
    end
    if ~isfield(scanParams{iScan}, 'stimToEVmatrix') || isempty(scanParams{iScan}.stimToEVmatrix) || ...
      ~isequal(size(scanParams{iScan}.stimToEVmatrix,1),nStims) || ~isequal(size(scanParams{iScan}.stimToEVmatrix,2),params.numberEVs)
        scanParams{iScan}.stimToEVmatrix = eye(nStims,params.numberEVs);
    end
    if isempty(EVnames) 
      EVnames = repmat({' '},1,params.numberEVs);
      if params.numberEVs < nStims
        EVnames = params.stimNames(1:params.numberEVs);
      else
        EVnames(1:params.numberEVs) = params.stimNames;
      end
    end

    paramsInfo = cell(1,nStims+4);
    paramsInfo{1}= {'scan', scanParams{iScan}.scan, 'editable=0','description of the scan'};
    for iEvent = 1:nStims
      paramsInfo{iEvent+1} = {fixBadChars(params.stimNames{iEvent}), scanParams{iScan}.stimToEVmatrix(iEvent,:),...
        'incdec=[-1 1]','incdecType=plusMinus','minmax=[0 inf]',...
        ['How much stimulus ' fixBadChars(params.stimNames{iEvent}) ' contributes to each EV']};
    end
    paramsInfo{nStims+2}= {'EVnames', EVnames, 'type=stringarray','self explanatory'};
    paramsInfo{nStims+3}= {'showDesign', 0, 'type=pushbutton','buttonString=Show Experimental Design','Plots the experimental design before convolution',...
            'callback',{@plotExperimentalDesign,scanParams,params,iScan,thisView},'passParams=1'};

    % give the option to use the same variable for remaining scans
    if (iScan ~= params.scanNum(end)) && (length(params.scanNum)>1)
      paramsInfo{nStims+4} = {'sameForNextScans',iScan == params.scanNum(1),'type=checkbox','Use the same parameters for all scans'};
    else
      paramsInfo(nStims+4) = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % now we have all the dialog information, ask the user to set parameters
    if useDefault
       tempParams = mrParamsDefault(paramsInfo);
    else
       tempParams = mrParamsDialog(paramsInfo,'Set Stimulus/EV matrix');
    end

    % user hit cancel
    if isempty(tempParams)
       scanParams = [];
       return
    end
    
    tempParams = mrParamsRemoveField(tempParams,'showDesign');
    
    % form stimToEV matrix from fields
    for iEvent = 1:nStims
      scanParams{iScan}.stimToEVmatrix(iEvent,:) = tempParams.(fixBadChars(params.stimNames{iEvent}));
      tempParams = mrParamsRemoveField(tempParams,fixBadChars(params.stimNames{iEvent}));
    end
    EVnames = tempParams.EVnames;
    tempParams = mrParamsRemoveField(tempParams,'EVnames');
    
    scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{iScan});

    % if sameForNextScans is set, copy all parameters into all remaining scans and break out of loop
    if isfield(scanParams{iScan},'sameForNextScans') && ...
       scanParams{iScan}.sameForNextScans
       for jScan = params.scanNum(find(params.scanNum>iScan,1,'first'):end)
          % set the other scans params to the same as this one
          scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{jScan});
       end
       break
    end
  %          paramsInfo = {};
    if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
      scanParams = [];
      return;
    end
  end
  params.EVnames = EVnames;
  keepAsking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plotStimDesign %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dummy = plotExperimentalDesign(thisScanParams,scanParams,params,scanNum,thisView)

dummy = [];
%we need to convert the stimToEVmatrix
for iEvent = 1:length(params.stimNames)
  thisScanParams.stimToEVmatrix(iEvent,:) = thisScanParams.(fixBadChars(params.stimNames{iEvent}));
end

if ~isfield(thisScanParams,'sameForNextScans') || thisScanParams.sameForNextScans
  %if it's the last scan or sameForNextScans is selected, we show all scans
  %copy params to remaining scans
  for jScan = params.scanNum(find(params.scanNum==scanNum,1,'first'):end)
    scanParams{jScan} = thisScanParams;
  end
  scanList = params.scanNum;
else
  %otherwise we only show the current scan
  scanParams = cell(1,scanNum);
  scanParams{scanNum} = thisScanParams;
  scanList = scanNum;
end

fignum = selectGraphWin(0,'Make new');
set(fignum,'name','plotExperimentalDesign');
nScans = length(scanList);
screenSize = get(0,'MonitorPositions');
screenSize = screenSize(1,:); % multiple screens
position = get(fignum,'position');
position(3) = screenSize(3);
position(4) = min(screenSize(3)/8*nScans,screenSize(4));
set(fignum,'position',position);

cScan=0;
axisLength = zeros(1,length(params.scanNum));
tSeriesAxes = zeros(1,length(params.scanNum));
for iScan = scanList
  params.scanParams{iScan} = mrParamsCopyFields(scanParams{iScan},params.scanParams{iScan});
  %replace all unused stimuli by one EV
  params.scanParams{iScan}.stimToEVmatrix(:,end+1) = ~any(params.scanParams{iScan}.stimToEVmatrix,2);
  params.scanParams{iScan}.EVnames{end+1} = 'Not used';
  
  d = loadScan(thisView, iScan, viewGet(thisView,'groupNum',params.groupName), 0);
  d = getStimvol(d,params.scanParams{iScan});
  [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,0,1);
  d = eventRelatedPreProcess(d,params.scanParams{iScan}.preprocess);
  d = makeDesignMatrix(d,params,1,iScan);
  cScan=cScan+1;
  
  if ~isempty(d.scm) && size(d.scm,2)==length(params.scanParams{iScan}.EVnames)*d.nHrfComponents
    %the length for the axis depends on the number of volumes times the frame period
    axisLength(cScan) = size(d.EVmatrix,1)*d.tr;
    colors = randomColors(length(params.scanParams{iScan}.EVnames));
    colors(end,:) = [.85 .85 .85]; %last color for unused stims

    tSeriesAxes(cScan) = axes('parent',fignum,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));
    hold on
    [h,hTransition] = plotStims(d.EVmatrix, d.stimDurations, d.tr/d.designSupersampling, colors, tSeriesAxes(cScan),d.runTransitions);

    title(viewGet(thisView,'description',iScan),'interpreter','none');
    %legend
    legendStrings = params.scanParams{iScan}.EVnames(h>0);
    h = h(h>0);
    if ~isempty(hTransition)
      h = [h hTransition];
      legendStrings = [legendStrings {'Run Transitions'}];
    end
    lhandle = legend(h,legendStrings,'position',getSubplotPosition(2,cScan,[7 1],ones(nScans,1),0,.05));
    set(lhandle,'Interpreter','none','box','off');
  else
    h = axes('parent',fignum,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));   
    set(h,'visible','off');
    text(.5,.5,{viewGet(thisView,'description',iScan), 'Number of EVs does not match', 'Change the stimToEV matrix for this scan'},...
      'parent',h,'units','normalized','HorizontalAlignment','center');
  end
end

%resize axes according to length of each sequence
maxAxisLength = max(axisLength);
for iScan = 1:cScan
  if axisLength(iScan)
    position = get(tSeriesAxes(iScan),'position');
    position(3) = position(3)*axisLength(iScan)/maxAxisLength;
    set(tSeriesAxes(iScan),'position',position);
  end
end
