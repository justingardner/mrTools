% searchlightClassification.m
%
%        $Id: searchlightClassification.m 1839 2010-11-14 17:45:36Z julien $
%      usage: view = eventRelated(view,params)
%         by: alex beckett
%       date: 10/20/06
%    purpose: event related volumetric searchlight
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = searchlightClassification(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = searchlightClassification(v,[],'justGetParams=1');
%
function [thisView, params] = searchlightClassification(thisView,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help searchlightClassification
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin));%,[],[],{'justGetParams','defaultParams','scanList'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end
if ieNotDefined('params'),params = [];end

% First get parameters
if isempty(params) || justGetParams
    params = searchlightClassGUI('thisView',thisView,'params',params,'defaultParams',defaultParams,'scanList',scanList);
%   else
%     params = searchlightClassGUI('groupName',viewGet(view,'groupName'),'scanList',scanList,'roilist',viewGet(view,'roiNames'));
%   end
end

% Abort if params empty
if ieNotDefined('params')
  disp('(searchlightAnalysis) Searchlight Analysis cancelled');
  return
% just return parameters
elseif justGetParams
  return
end


% set the group
thisView = viewSet(thisView,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

params.fweMethod='Adaptive Step-down';
params.fdrAssumption= 'Independence/Positive dependence';
params.fdrMethod= 'Adaptive Step-up';
% params.fweAdjustment= 1;
% params.fdrAdjustment =  1;
params.testOutput= 'P value';
params.trueNullsEstimationMethod= 'Least Squares';
params.trueNullsEstimationThreshold= 0.0500;

params.offsets=make_sphere(params.radius);
scanParams = params.scanParams;
tic



%-------------------Main Loop Over Scans----------
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
%initialize the data we're keeping for output overlays
precision = mrGetPref('defaultPrecision');
accDiag_Class = cell(1,params.scanNum(end));
accDiag = cell(1,params.scanNum(end));
if params.sigTest
    pDiag= cell(1,params.scanNum(end));
    pDiag_Class = cell(1,params.scanNum(end));
end

%scanloop
for scanNum = params.scanNum
    scanDims = viewGet(thisView,'dims',scanNum);
    [d, d.roiVoxelIndices, d.roiCoords] = loadScanRois(thisView,scanNum,viewGet(thisView,'roinum',params.roiMask));
    
    d = getStimvol(d,scanParams{scanNum});
    if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    
    % do any call for preprocessing
    if ~isempty(scanParams{scanNum}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
    end

    for i=1:size(scanParams{scanNum}.stimToEVmatrix,2)
        tmp.stimvol{i}=[d.stimvol{[find(scanParams{scanNum}.stimToEVmatrix(:,i))]}];
    end
    d.stimvol=tmp.stimvol;
    d.stimNames=scanParams{scanNum}.EVnames;
    
    lab = [];
    for i=1:length(d.stimvol)
        lab(d.stimvol{i})=i;
    end

    run=[];
    for i=1:size(d.concatInfo.runTransition,1)
    run(d.concatInfo.runTransition(i,1):d.concatInfo.runTransition(i,2))=i;
    end
    
    d.data=squeeze(d.data);
    
    %pick out the eventstrings and average if requested
    idx = find(lab>0);
    if scanParams{scanNum}.eventLength==1
        m_ = d.data(:,idx+scanParams{scanNum}.hdLag);
        run=run(idx);
        lab=lab(idx);
    elseif scanParams{scanNum}.averageEvent %average across the TRs for each event
        for i=1:size(idx,2)
            m_(:,i)=mean(d.data(:,idx(i)+scanParams{scanNum}.hdLag:idx(i)+scanParams{scanNum}.eventLength+scanParams{scanNum}.hdLag-1),2);
        end
%         d.roi{1}.tSeries=m_;clear m_
        run=run(idx);
        lab=lab(idx);
    elseif scanParams{scanNum}.eventLength>1 %create instance from each TR in the stim duration
        for i=1:size(idx,2)
            m_(:,idx(i):idx(i)+scanParams{scanNum}.eventLength-1)=d.data(:,idx(i)+scanParams{scanNum}.hdLag:idx(i)+scanParams{scanNum}.eventLength+scanParams{scanNum}.hdLag-1);
            l_(idx(i):idx(i)+scanParams{scanNum}.eventLength-1)=repmat(lab(idx(i)),1,scanParams{scanNum}.eventLength);
            r_(idx(i):idx(i)+scanParams{scanNum}.eventLength-1)=repmat(run(idx(i)),1,scanParams{scanNum}.eventLength);
        end
        lab=l_;
        run=r_;
        clear l_ r_
    end
    
    % works out which roi coords are indexed by the spotlights
    disppercent(-inf, '(searchlightClassification) Creating 4D TSeries....');
    for i=1:length(d.roiCoords{1})
        mm(d.roiCoords{1}(1,i),d.roiCoords{1}(2,i),d.roiCoords{1}(3,i),:) = m_(i,:);
        disppercent(i/length(d.roiCoords{1}));
    end
    disppercent(inf);
    
    %initialise overlays per scan
    accDiag{scanNum}=nan(viewGet(thisView,'scanDims',scanNum));
    accDiag_Class{scanNum}=nan([viewGet(thisView,'scanDims',scanNum),length(d.stimvol)]);
    
    if params.sigTest
        pDiag{scanNum}=nan(viewGet(thisView,'scanDims',scanNum));
        pDiag_Class{scanNum}=nan([viewGet(thisView,'scanDims',scanNum),length(d.stimvol)]);
    end
    
    disppercent(-inf,'(searchlightClassification) Classifying based on spotlight....');
    
    for i_sphere=1:size(d.roiCoords{1},2)
    
        %find voxels in spotlight
        idx=repmat(d.roiCoords{1}(:,i_sphere),1,size(params.offsets,2))+params.offsets;
        [~,j] = find((idx(1,:)<min(d.roiCoords{1}(1,:)) | idx(1,:)>max(d.roiCoords{1}(1,:))) | (idx(2,:)<min(d.roiCoords{1}(2,:)) | idx(2,:)>max(d.roiCoords{1}(2,:))) | (idx(3,:)<min(d.roiCoords{1}(3,:)) | idx(3,:)>max(d.roiCoords{1}(3,:))));
        idx=idx(:,setdiff(1:size(idx,2),j));

        %create patterns based on spotlight
        for i_vol=1:size(idx,2)
            xxx(i_vol,:)=mm(idx(1,i_vol),idx(2,i_vol),idx(3,i_vol),:);
        end
        [I,~]=find(xxx>0);
        xxx=xxx(unique(I),:);

        
        %run classification on searchlight patterns
        class_lab=nan(1,length(run));
        corr=nan(1,length(run));
        
        for i=1:size(d.concatInfo.runTransition,1)
            class_lab(run==i)=classify(xxx(:,run==i)',xxx(:,run~=i)',lab(run~=i),'diagLinear')';
            corr(run==i)=class_lab(run==i)==lab(run==i);    
        end
        
        if params.sigTest
            pDiag{scanNum}(d.roiCoords{1}(1,i_sphere),d.roiCoords{1}(2,i_sphere),d.roiCoords{1}(3,i_sphere)) = binomTest(sum(corr),length(corr),1/length(d.stimvol),'Greater');
            for i=1:length(d.stimvol)
                pDiag_Class{scanNum}(d.roiCoords{1}(1,i_sphere),d.roiCoords{1}(2,i_sphere),d.roiCoords{1}(3,i_sphere),i) =  binomTest(sum(corr(lab==i)),length(corr(lab==i)),1/length(d.stimvol),'Greater');
            end
        end
        
        %calculate global and class specific accuracies
        accDiag{scanNum}(d.roiCoords{1}(1,i_sphere),d.roiCoords{1}(2,i_sphere),d.roiCoords{1}(3,i_sphere))=mean(corr);
        for i=1:length(d.stimvol)
            accDiag_Class{scanNum}(d.roiCoords{1}(1,i_sphere),d.roiCoords{1}(2,i_sphere),d.roiCoords{1}(3,i_sphere),i)=mean(corr(lab==i));
        end
        

        xxx=[];

        disppercent(i_sphere/length(d.roiCoords{1}));
    end
      
%     [acc{scanNum} class_acc{scanNum} acc_p{scanNum}]=classify_searchlight(thisView,d,params.radius,params.roiMask,params.scanParams{scanNum});
    
    if params.sigTest
        [pDiag{scanNum} fdr{scanNum} fwe{scanNum}] = transformStatistic(pDiag{scanNum},precision,params);
        [pDiag_Class{scanNum} fdr_Class{scanNum} fwe_Class{scanNum}] = transformStatistic(pDiag_Class{scanNum},precision,params);
    end

end



%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
classAnal.name = params.saveName;
classAnal.type = 'searchClass';
classAnal.groupName = params.groupName;
classAnal.function = 'searchlightClasfficiation';
classAnal.reconcileFunction = 'defaultReconcileParams';
classAnal.mergeFunction = 'defaultMergeParams';
classAnal.guiFunction = 'searchlightClassGUI';
classAnal.params = params;
classAnal.date = dateString;


%--------------------------------------------------------- Output overlay structures
nScans = viewGet(thisView,'nScans');
% create generic parameters 
defaultOverlay.groupName = params.groupName;
defaultOverlay.function = 'searchlightClasfficiation';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.date = dateString;
defaultOverlay.params = cell(1,nScans);
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = 'timecoursePlot';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [0 1];
defaultOverlay.clip = [0 1];
defaultOverlay.alphaOverlay='';
defaultOverlay.alphaOverlayExponent=1;
defaultOverlay.data = cell(1,nScans);
defaultOverlay.name = '';
for iScan = params.scanNum
   defaultOverlay.data{iScan} = NaN(scanDims); %to make values outside the box transparent
end

%------------------------------------save overlay
overlays = defaultOverlay;
overlays.name = ['mean accuracy (radius = ',num2str(params.radius),')'];

for iScan = params.scanNum
   overlays.data{iScan}=accDiag{iScan};
   overlays.params{iScan} = params.scanParams{iScan};
end

thisOverlay = defaultOverlay;
for i=1:params.numberEvents
        overlays(end+1)=thisOverlay;
        overlays(end).name=[params.EVnames{i},'_acc (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=accDiag_Class{iScan}(:,:,:,i);
            overlays(end).params=params.scanParams{iScan};
        end
end

if params.sigTest
    thisOverlay=defaultOverlay;
    overlays(end+1)=thisOverlay;
    overlays(end).colormap = statsColorMap(256);
    overlays(end).name=['probability (radius = ',num2str(params.radius),')'];
    for iScan = params.scanNum
        overlays(end).data{iScan}=pDiag{iScan};
        overlays(end).params=params.scanParams{iScan};
    end
    for i=1:params.numberEvents
        overlays(end+1)=thisOverlay;
        overlays(end).colormap = statsColorMap(256);
        overlays(end).name=[params.EVnames{i},'_prob (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=pDiag_Class{iScan}(:,:,:,i);
            overlays(end).params=params.scanParams{iScan};
        end
    end
    if params.fdrAdjustment
        overlays(end+1)=thisOverlay;
        overlays(end).colormap = statsColorMap(256);
        overlays(end).name=['fdr probability (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=fdr{iScan};
            overlays(end).params=params.scanParams{iScan};
        end
        for i=1:params.numberEvents
            overlays(end+1)=thisOverlay;
            overlays(end).colormap = statsColorMap(256);
            overlays(end).name=[params.EVnames{i},'_fdr_prob (radius = ',num2str(params.radius),')'];
            for iScan = params.scanNum
                overlays(end).data{iScan}=fdr_Class{iScan}(:,:,:,i);
                overlays(end).params=params.scanParams{iScan};
            end
        end
    end
    if params.fweAdjustment
        overlays(end+1)=thisOverlay;
        overlays(end).colormap = statsColorMap(256);
        overlays(end).name=['fwe probability (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=fwe{iScan};
            overlays(end).params=params.scanParams{iScan};
        end
        for i=1:params.numberEvents
            overlays(end+1)=thisOverlay;
            overlays(end).colormap = statsColorMap(256);
            overlays(end).name=[params.EVnames{i},'_fwe_prob (radius = ',num2str(params.radius),')'];
            for iScan = params.scanNum
                overlays(end).data{iScan}=fwe_Class{iScan}(:,:,:,i);
                overlays(end).params=params.scanParams{iScan};
            end
        end
    end
end


classAnal.overlays = overlays;
thisView = viewSet(thisView,'newAnalysis',classAnal);
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
toc

%-------------------------------------------------------- Save the analysis
saveAnalysis(thisView,classAnal.name);

oneTimeWarning('nonZeroHrfStart',0);
oneTimeWarning('tfceOutputsZeros',0);
set(viewGet(thisView,'figNum'),'Pointer','arrow');
refreshMLRDisplay(viewGet(thisView,'viewNum'));



return