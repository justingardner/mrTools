% eventRelatedGlm.m
%
%      usage: view = eventRelatedGlm(view,params)
%         by: farshad moradi
%       date: 06/14/07
%    purpose: same as eventRelated, but uses canonical hrf instead of
%             deconvolution
%
function [view d] = eventRelatedGlm(view,params)

d = [];

% check arguments
if ~any(nargin == [1 2])
  help eventRelated
  return
end

mrGlobals;

% First geteventRelatedGlmReconcileParams parameters
if ieNotDefined('params')
  % put up the gui
  params = eventRelatedGlmGUI('groupName',viewGet(view,'groupName'));
end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = mrParamsReconcile([],params);

% Abort if params empty
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

if ~isempty(params.contrast)
    contrast = str2num(params.contrast);
    % check if the contrast is defined correctly
    if isempty(contrast),
        mrErrorDlg('invalid contrast. must be a vector');
    end
else
    contrast = [];
end

% create the parameters for the overlay
dateString = datestr(now);
if isempty(contrast)
    r2.name = 'r2';
else
    r2.name = 'r2c';
end
r2.groupName = params.groupName;
r2.function = 'eventRelatedGlm';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(view,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.interrogator = 'glmContrastPlot';
r2.mergeFunction = 'defaultMergeParams';

if ~isempty(contrast)
    amps = cell(1,size(contrast,1));
    for i = 1:size(contrast,1)
        amp = r2;
        amp.name = ['c ', num2str(contrast(i,:))];
        amp.range = [-0.25 0.25];
        amp.clip = [-0.25 0.25];
        amp.colormap = jet(256);
        amps{i} = amp;
    end
end

r2.colormapType = 'setRangeToMax';

tic
set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
  % decide how many slices to do at a time, this is done
  % simply to save memory -- currently our system is limited
  % to 2G of memory and for large concatenations, you need
  % to break up the analysis into smaller portions of the data
  numSlices = viewGet(view,'nSlices',scanNum);
  numVolumes = viewGet(view,'nFrames',scanNum);
  dims = viewGet(view,'dims',scanNum);
  % choose how many slices based on trying to keep a certain
  % amount of data in the memory
  numSlicesAtATime = getNumSlicesAtATime(numVolumes,dims);
  currentSlice = 1;
  ehdr = [];ehdrste = [];thisr2 = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % load the scan
    d = loadScan(view,scanNum,[],[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);
    params.hrfParams.tmax = params.scanParams{scanNum}.hdrlen+d.tr/2;
    
    d.supersampling = params.trSupersampling;
    % use the duration of stimuli/events in the design matrix
    d.impulse = 0; 

    % get the stim volumes, if empty then abort
    d = getStimvol(d,params.scanParams{scanNum});
    if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    
    % do any called for preprocessing
    hrf = feval(params.hrfModel, d.tr/d.supersampling, params.hrfParams);
    
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
    % make a stimulation convolution matrix
    
    d = makeglm(d,hrf);
    
    if isempty(contrast)
        % compute the estimated hemodynamic responses
        d = getr2(d);
    else
        % check that contrast has the same number of columns as the
        % design matrix
        if size(d.scm,2)~=size(contrast,2)*size(hrf,2)
            ehdr = [];
            ehdrste = [];
            thisr2 = [];
            mrErrorDlg( sprintf('contrast incompatible with scan %d', scanNum) );
        end
        % expand the contrast to a matrix with one row per hrf column
        nhr = size(hrf,2);
        c = zeros( size(contrast)*nhr );
        for h=1:nhr
            c(h:nhr:end, h:nhr:end) = contrast;
        end
        % compute the estimated hemodynamic responses for the given
        % contrast, discounting any effect of other orthogonal contrasts
        d = getGlmContrast(d, c);
    end
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    ehdrste = cat(3,ehdrste,d.ehdrste);
    thisr2 = cat(3,thisr2,d.r2);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = ehdr;
  d.ehdrste = ehdrste;
  d.r2 = thisr2;

  d.dim(3) = size(d.r2,3);

  % save the r2 overlay
  r2.data{scanNum} = d.r2;
  r2.params{scanNum} = params.scanParams{scanNum};
  
  if ~isempty(contrast)
      stimNames = cell(1,size(contrast,1));
      % for all contrasts report the amplitude of the 1st component of the hrf
      % first, calculate the amplitude of modulation of the 1st component
      temp = d.simulatedhrf(:,1);
      sc = max(temp)-min(temp);
      for j = 1:size(contrast,1)
          temp = squeeze(d.ehdr(:,:,:,j,1))*sc;
          amps{j}.data{scanNum} = temp;
          mx = max(temp(:));
          mn = min(temp(:));
          amps{j}.clip = [ min([mn,amps{j}.clip(1)]), max([amps{j}.clip(2),mx])];
          amps{j}.params{scanNum} = params.scanParams{scanNum};
          
          stimNames{j} = ['c ', num2str(contrast(j,:))];
      end
  else
      stimNames = d.stimNames;
  end

  
  % save other eventRelated parameters
  erAnal.d{scanNum}.hrf = d.simulatedhrf;
  erAnal.d{scanNum}.actualhrf = hrf;
  erAnal.d{scanNum}.trsupersampling = d.supersampling;
  erAnal.d{scanNum}.ver = d.ver;
  erAnal.d{scanNum}.filename = d.filename;
  erAnal.d{scanNum}.filepath = d.filepath;
  erAnal.d{scanNum}.dim = d.dim;
  erAnal.d{scanNum}.ehdr = d.ehdr;
  erAnal.d{scanNum}.ehdrste = d.ehdrste;
  erAnal.d{scanNum}.nhdr = d.nhdr;
  erAnal.d{scanNum}.hdrlen = d.hdrlen;
  erAnal.d{scanNum}.tr = d.tr;
  erAnal.d{scanNum}.stimNames = stimNames;
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;

  stimvol = d.stimvol;
  for i=1:length(stimvol)
      stimvol{i} = unique(ceil(stimvol{i}/d.supersampling));
  end
  erAnal.d{scanNum}.stimvol = stimvol;

end
toc

% install analysis
erAnal.name = params.saveName;
if isempty(contrast)
    erAnal.type = 'glmAnal';
else
    erAnal.type = 'glmcAnal';
end
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelatedGlm';
erAnal.reconcileFunction = 'defaultReconcileParams';
erAnal.mergeFunction = 'defaultMergeParams';
erAnal.guiFunction = 'eventRelatedGlmGUI';
erAnal.params = params;
% erAnal.overlays = r2;
% erAnal.curOverlay = 1;
erAnal.date = dateString;
view = viewSet(view,'newAnalysis',erAnal);
view = viewSet(view,'newoverlay',r2);

if ~isempty(contrast)
    for i = 1:size(contrast,1)
        view = viewSet(view,'newoverlay',amps{i});
    end
end

% Save it
saveAnalysis(view,erAnal.name);

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

