% roiClassification.m
%
%        $Id: roiClassification.m 1839 2010-11-14 17:45:36Z julien $
%      usage: view = roiClassification(view,params)
%         by: alex beckett
%       date: 10/20/06
%    purpose: roi based classification data analysis
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = roiClassification(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = roiClassification(v,[],'justGetParams=1');
%
function [view d] = roiClassification_old(view,params,varargin)

d = [];

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8])
  help roiClassification.m
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin,[],[],{'justGetParams','defaultParams','scanList'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end

roi_n = selectROIs(view);
if isempty(roi_n)
    mrWarnDlg('(roiClassification) No ROI selected!');
  return
end

% First get parameters
if ieNotDefined('params')
  % put up the gui
  if defaultParams
    params = roiClassGUI('groupName',viewGet(view,'groupName'),'useDefault=1','scanList',scanList);
  else
    params = roiClassGUI('groupName',viewGet(view,'groupName'),'scanList',scanList);
  end
end

% just return parameters
if justGetParams
  d = params;
  return
end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = mrParamsReconcile([],params);

% Abort if params empty
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

if ~isfield(params, 'applyFiltering')
  params.applyFiltering = 0;
end

% inplace concatenation is handeled by a different function
if isfield(params, 'inplaceConcat')
    if params.inplaceConcat
        [view d] = eventRelatedMultiple(view,params);
        return
    end
end


% for scanNum = params.scanNum
%     d = loadScan(view,scanNum,[],0)
%     d = getStimvol(d,params.scanParams{scanNum});
%     m_d=classify_sphere_beta(radius,d)
% end

% create the parameters for the overlay
dateString = datestr(now);
acc.name = 'acc';
acc.groupName = params.groupName;
acc.function = 'roiClassification';
acc.reconcileFunction = 'defaultReconcileParams';
acc.data = cell(1,viewGet(view,'nScans'));
acc.date = dateString;
acc.params = cell(1,viewGet(view,'nScans'));
acc.range = [0 1];
acc.clip = [0 1];
% colormap is made with a little bit less on the dark end
acc.colormap = hot(312);
acc.colormap = acc.colormap(end-255:end,:);
acc.alpha = 1;
acc.colormapType = 'setRangeToMax';
acc.interrogator = 'timecoursePlot';
acc.mergeFunction = 'defaultMergeParams';

tic
set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
    d = loadScan(view,scanNum,[],0);
    d = getStimvol(d,params.scanParams{scanNum});
     % do any called for preprocessing
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
      
    acc.data{scanNum}=classify_roi(view,d,params.scanParams{scanNum},roi_n,params.select_vox,params.numShuff);
    acc.params{scanNum} = params.scanParams{scanNum};
   
    %so far this is as far as we go, ends at keyboard propmt in the above
    %function. Is it worth making an overlay for this analysis? accuracy
    %for each ROI? output weights of the classifiers? or just display the
    %rsults and save a file with the outputs?
    

  % get the actual size of the data (not just the size of the last
  % slice/set of rows we were working on).
  d.dim(1:3) = size(acc.data{scanNum});

  % save the r2 overlay
%   r2.data{scanNum} = d.r2;
%   r2.params{scanNum} = params.scanParams{scanNum};
  
  % save other eventRelated parameters
%   erClass.d{scanNum}.ver = d.ver;
%   erClass.d{scanNum}.filename = d.filename;
%   erClass.d{scanNum}.filepath = d.filepath;
%   erClass.d{scanNum}.dim = d.dim;
%   erClass.d{scanNum}.ehdr = d.ehdr;
%   erClass.d{scanNum}.ehdrste = d.ehdrste;
%   erClass.d{scanNum}.nhdr = d.nhdr;
%   erClass.d{scanNum}.hdrlen = d.hdrlen;
%   erClass.d{scanNum}.tr = d.tr;
%   erClass.d{scanNum}.stimvol = d.stimvol;
%   erClass.d{scanNum}.stimNames = d.stimNames;
%   erClass.d{scanNum}.scm = d.scm;
%   erAnal.d{scanNum}.expname = d.expname;
%   erAnal.d{scanNum}.fullpath = d.fullpath;
end
toc

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
view = viewSet(view,'newAnalysis',erClass);
if ~isempty(viewGet(view,'fignum'))
  refreshMLRDisplay(viewGet(view,'viewNum'));
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

