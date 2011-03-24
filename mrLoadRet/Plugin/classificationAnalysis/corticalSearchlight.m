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
function [view d] = corticalSearchlight(view,params,varargin)
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
    params = surfaceClassGUI('thisView',view,'params',params,'defaultParams',defaultParams,'scanList',scanList);
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
%   return
end


% set the group
view = viewSet(view,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

% for scanNum = params.scanNum
%     d = loadScan(view,scanNum,[],0)
%     d = getStimvol(d,params.scanParams{scanNum});
%     m_d=classify_sphere_beta(radius,d)
% end

% create the parameters for the overlay
dateString = datestr(now);
acc.name = ['acc_',num2str(params.radius)];
acc.groupName = params.groupName;
acc.function = 'eventRelated';
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
acc.interrogator = 'timeSeriesPlot';
acc.mergeFunction = 'defaultMergeParams';

tic
set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
    d = loadScan(view,scanNum,[],0);
    d = getStimvol(d,params.scanParams{scanNum});
     % do any called for preprocessing
      d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
      
    acc.data{scanNum}=classify_corticalSearchlight(view,d,params.radius,params.scanParams{scanNum});
    acc.params{scanNum} = params.scanParams{scanNum};
    
    
%   % decide how many slices to do at a time, this is done
%   % simply to save memory -- currently our system is limited
%   % to 2G of memory and for large concatenations, you need
%   % to break up the analysis into smaller portions of the data
%   numSlices = viewGet(view,'nSlices',scanNum);
%   numVolumes = viewGet(view,'nFrames',scanNum);
%   dims = viewGet(view,'dims',scanNum);
%   % choose how many slices based on trying to keep a certain
%   % amount of data in the memory
%   [numSlicesAtATime rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,dims);
%   currentSlice = 1;
%   ehdr = [];ehdrste = [];thisr2 = [];
% 
%   for i = 1:ceil(numSlices/numSlicesAtATime)
%     % calculate which slices we are working on
%     thisSlices = [currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)];
%     % set the row we are working on
%     currentRow = 1;
%     % clear variables that will hold the output for the slices we
%     % are working on
%     sliceEhdr = [];sliceEhdrste = [];sliceR2 = [];
%     for j = 1:ceil(dims(1)/numRowsAtATime)
%       % load the scan
%       thisRows = [currentRow min(dims(1),currentRow+numRowsAtATime-1)];
%       d = loadScan(view,scanNum,[],thisSlices,precision,thisRows);
%       % get the stim volumes, if empty then abort
%       d = getStimvol(d,params.scanParams{scanNum});
%       if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
%       % do any called for preprocessing
%       d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
%       % make a stimulation convolution matrix
%       d = makescm(d,ceil(params.scanParams{scanNum}.hdrlen/d.tr),params.applyFiltering);
%       % compute the estimated hemodynamic responses
%       d = getGlmStatistics(d);
%       % update the current row we are working on
%       currentRow = currentRow+numRowsAtATime;
%       % if we are calculating full slice, then just pass that on
%       if numRowsAtATime == dims(1)
% 	sliceEhdr = d.ehdr;
% 	sliceEhdrste = d.ehdrste;
% 	sliceR2 = d.r2;
%       % working on a subset of rows, cat together with what
%       % has been computed for other rows
%       else
% 	sliceEhdr = cat(1,sliceEhdr,d.ehdr);
% 	sliceEhdrste = cat(1,sliceEhdrste,d.ehdrste);
% 	sliceR2 = cat(1,sliceR2,d.r2);
%       end
%     end
%     % update the current slice we are working on
%     currentSlice = currentSlice+numSlicesAtATime;
%     % cat with what has already been computed for other slices
%     ehdr = cat(3,ehdr,sliceEhdr);
%     ehdrste = cat(3,ehdrste,sliceEhdrste);
%     thisr2 = cat(3,thisr2,sliceR2);
%   end

  % now put all the data from all the slices into the structure
%   d.ehdr = single(ehdr);
%   d.ehdrste = single(ehdrste);
%   d.r2 = single(thisr2);

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

