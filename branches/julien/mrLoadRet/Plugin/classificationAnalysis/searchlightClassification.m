% searchlightClassification.m
%
%        $Id: searchlightClassification.m 1839 2010-11-14 17:45:36Z julien $
%      usage: view = eventRelated(view,params)
%         by: alex beckett
%       date: 10/20/06
%    purpose: event related data analysis
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

% 
% if ~isfield(params, 'applyFiltering')
%   params.applyFiltering = 0;
% end

% % inplace concatenation is handeled by a different function
% if isfield(params, 'inplaceConcat')
%     if params.inplaceConcat
%         [thisView d] = eventRelatedMultiple(thisView,params);
%         return
%     end
% end


% for scanNum = params.scanNum
%     d = loadScan(view,scanNum,[],0)
%     d = getStimvol(d,params.scanParams{scanNum});
%     m_d=classify_sphere_beta(radius,d)
% end

% create the parameters for the overlay
% dateString = datestr(now);
% acc.name = ['acc_',num2str(params.radius)];
% acc.groupName = params.groupName;
% acc.function = 'searchClass';
% acc.reconcileFunction = 'defaultReconcileParams';
% acc.data = cell(1,viewGet(thisView,'nScans'));
% acc.date = dateString;
% acc.params = cell(1,viewGet(thisView,'nScans'));
% acc.range = [0 1];
% acc.clip = [0 1];
% % colormap is made with a little bit less on the dark end
% acc.colormap = hot(312);
% acc.colormap = acc.colormap(end-255:end,:);
% acc.alpha = 1;
% acc.colormapType = 'setRangeToMax';
% acc.interrogator = 'timeSeriesPlot';
% acc.mergeFunction = 'defaultMergeParams';

tic
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
    scanDims = viewGet(thisView,'dims',scanNum);
    d = loadScan(thisView,scanNum,[],0);
    d = getStimvol(d,params.scanParams{scanNum});
    
%   do any called for preprocessing
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
      
    [acc{scanNum} class_acc{scanNum}]=classify_searchlight(thisView,d,params.radius,params.roiMask,params.scanParams{scanNum});
    

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
   overlays.data{iScan}=acc{iScan};
   overlays.params{iScan} = params.scanParams{iScan};
end
thisOverlay = defaultOverlay;
for i=1:params.numberEvents
        overlays(end+1)=thisOverlay;
        overlays(end).name=[params.EVnames{i},'_acc (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=class_acc{iScan}(:,:,:,i);
            overlays(end).params=params.scanParams{iScan};
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










% install analysis
erClass.name = params.saveName;
erClass.type = 'erClass';
erClass.groupName = params.groupName;
erClass.function = 'eventRelatedClassification';
erClass.reconcileFunction = 'defaultReconcileParams';
erClass.mergeFunction = 'defaultMergeParams';
erClass.guiFunction = 'eventRelatedClassGUI';
erClass.params = params;
erClass.overlays = acc;
erClass.curOverlay = 1;
erClass.date = dateString;
view = viewSet(thisView,'newAnalysis',erClass);
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end

% Save it
saveAnalysis(view,erClass.name);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    erAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(erAnal.d) == 1
    d = erAnal.d{1};
  else
    d = erAnal.d;
  end
end

