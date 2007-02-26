% eventRelatedMultiple.m
%
%      usage: view = eventRelatedMultiple(view,params)
%         by: farshad moradi (based on code by justin gardner)
%       date: 02/05/07
%    purpose: eventrelated for multiple files
%
function view = eventRelatedMultiple(view,params)

% check arguments
if ~any(nargin == [1 2])
  help(mfilename)
  return
end

% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
%  params = eventRelatedGUI('groupName',viewGet(view,'groupName'));
  params.scanNum = selectScans(view);
  params = eventRelatedReconcileParams(viewGet(view,'groupName'),params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params = eventRelatedReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('eventRelated cancelled');
  return
end

if ~isfield(params, 'includeScans')
    params.includeScans = params.scanNum;
end;

% step 1: estimate hrf using all included scans
ehdr = 0;
[scms, volumes, nhdr] = getFullDesign(view, params);
if nhdr<=0
    mrMsgBox('incompatible number of conditions');
    return
end

for i=params.includeScans
    % load the scan
    d = loadScan(view,i);
    d.scms = scms;
    d.volumes = volumes{i};
    % compute the hemodynamic responses
    d = calcleastsqestimate(d);
    ehdr = ehdr+d.ehdr;
end
% step 2: calculate r2
unexplainedVariance = 0;
totalVariance = 0;

for i=params.includeScans
    d = loadScan(view,i);
    d.ehdr = ehdr;
    d.scm = scms{i};
    d.volumes = volumes{i};
    d = calcVariances(d);
    unexplainedVariance = unexplainedVariance+d.unexplainedVariance;
    totalVariance = totalVariance+d.totalVariance;
end

d = rmfield(d,'data');
d.nhdr = nhdr;
d.hdrlen = ceil(params.hdrlen/d.tr);
d.r2 = 1-unexplainedVariance./totalVariance;

% create the r2 overlay
dateString = datestr(now);
r2.name = 'co';
r2.function = 'corAnal';
r2.groupName = params.groupName;
r2.reconcileFunction = 'eventRelatedReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
for i=params.includeScans(1)
    r2.data{i}=d.r2;
end;
r2.date = dateString;
r2.params = params;
r2.range = [0 1];
r2.clip = [0 1];
r2.colormap = hot(256);
r2.alpha = 1;
r2.interrogator = 'eventRelatedPlot';

% install analysis
erAnal.name = 'erAnal';  % This can be reset by editAnalysisGUI
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'eventRelatedReconcileParams';
erAnal.guiFunction = 'eventRelatedGUI';
erAnal.params = params;
erAnal.overlays = r2;
erAnal.curOverlay = 1;
erAnal.date = dateString;
erAnal.ehdr = ehdr;
erAnal.nhdr = nhdr;
erAnal.d = d;

view = viewSet(view,'newAnalysis',erAnal);

% Save it
if ~isfield(params, 'suppress_save')
    saveAnalysis(view,erAnal.name);
end


%%%%%%%%%%%%%%%%%%%
% loadScan
%%%%%%%%%%%%%%%%%%%
function d = loadScan(view,scanNum)

% load parameters
d.tr = viewGet(view,'framePeriod',scanNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum);
d.dim = viewGet(view,'scandims',scanNum);
d.nFrames = viewGet(view,'nFrames',scanNum);
d.dim(4) = d.nFrames;
d.scanNum = scanNum;
% Load data
mrDisp(sprintf('Loading scan %i from %s\n',scanNum,viewGet(view,'groupName')));
d.data = loadTSeries(view,scanNum,'all');
% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum);
d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);
% load dicom header
d.dicom = viewGet(view,'dicom',scanNum);
d.baseline = ones(d.dim);

%%%%%%%%%%%%%%%%%%%
% getFullDesign
%%%%%%%%%%%%%%%%%%%
function [scms, volumes, nhdr] = getFullDesign(view, params)
scms=[];
volumes = [];
nhdr = 0;
for i=params.includeScans
    stimfile = viewGet(view,'stimfile',i);
    if length(stimfile) == 1
      stimfile = stimfile{1};
      switch stimfile.filetype,
          case 'traces',
            % get acquisition times. 2 means volume acq (we don't know slice acq)
            stimfile.acq = [0 2*(diff(d.traces(1,:))==1)];
          case 'eventtimes',
          otherwise,
            mrMsgBox(sprintf('Unknown of invalid stimfile type %s',stimfile.filetype));
      end
    else
        disp(i);
        error('invalid number of stimulus files');
    end
    d.filetype = stimfile.filetype;
    stimfile.tr = viewGet(view,'framePeriod',i);
    d = getStimVol(stimfile);
    d.dim = viewGet(view,'scandims',i);
    d.nFrames = viewGet(view,'nFrames',i);
    d.dim(4) = d.nFrames;
    d = makescm(d,ceil(params.hdrlen/d.tr));
    scms{i}=d.scm;
    volumes{i}=d.volumes;
    if nhdr==0
        nhdr=d.nhdr;
    end
    if nhdr~=d.nhdr
        nhdr = -1;
        return
    end
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getstimvol.m
%
%      usage: d = getStimVol(d)
%         by: justin gardner
%       date: 12/21/05
%    purpose: gets the stimulus vols from the traces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getStimVol(d,keepstim)

% check arguments
if ~any(nargin == [1 2])
  help getstimvol
  return
end
if exist('keepstim','var')~=1,keepstim = []; end;
switch d.filetype,
    case 'traces',
        d = getStimVolFromTraces(d);
    case 'eventtimes',
        d = getStimVolFromEventTimes(d);
end;
% get the first volume
d = getfirstvol(d);

% check if we need to only keep some vols
if ~isempty(keepstim)
  if (min(keepstim) < 1) || (max(keepstim) > length(d.stimvol))
    mrDisp(sprintf('UHOH: Keepstim out of range. Ignoring'));
  else
    % go through and only keep the stim values asked for
    stimvol = [];
    for i = 1:length(keepstim)
      stimvol{i} = d.stimvol{keepstim(i)};
    end
    d.stimvol = stimvol;
  end
end


function d = getStimVolFromTraces(d)
% get the stim times
stimraw = d.traces(d.stimtrace,:);
stimraw(stimraw < 0) = 0;
stimtimes = find([0 (diff(stimraw)~=0)]);

% get the image number
acqnum = cumsum(d.acq>1);

% set the beginning acqnum to 1, so that
% any event that happens before the first
% acquistion pulse is assumed to happen
% during the first acquisition pulse.
acqnum(1:first(find(acqnum == 1))) = 1;

% sort into stimuli
nhdr = max(stimraw);
for i = 1:nhdr
  d.stimtimes{i} = stimtimes(stimraw(stimtimes) == i);
  d.pulselens(i) = i;
  d.stimvol{i} = acqnum(d.stimtimes{i});
end


function d = getStimVolFromEventTimes(d)
% sort into stimuli
nhdr = length(d.stimtimes_s);
for i = 1:nhdr
  d.stimtimes{i} = d.stimtimes_s{i};
  d.pulselens(i) = i;
  d.stimvol{i} = round(d.stimtimes_s{i} / d.tr);
end

%%%%%%%%%%%%%%%%%%%
% calcVariances
%
%      usage: d = calcVariances(d)
%         by: farshad moradi
%       date: 02/05/07
%    purpose: calculates unexplained and total variance
%
%%%%%%%%%%%%%%%%%%%
function d = calcVariances(d)

% no roi, so just do whole volume
slices = 1:d.dim(3);
yvals = 1:d.dim(2);
yvaln = length(yvals);
% initialize values
uv = [];
tv = [];
% preallocate memory
d.unexplainedVariance = zeros(d.dim(1),d.dim(2),d.dim(3));
d.totalVariance = zeros(d.dim(1),d.dim(2),d.dim(3));
% display string
disppercent(-inf,'Calculating goodness of fit');
% cycle through images calculating the estimated hdr and r^s of the 
% estimate.
%
onesmatrix = ones(length(d.volumes),1);
for j = yvals
    disppercent(max((j-min(yvals))/yvaln,0.1));
    for k = slices
        ehdr = squeeze(d.ehdr(:,j,k,:))';
        % get the time series we are working on
        timeseries = squeeze(d.data(:,j,k,d.volumes))';
        % subtract off column means
        colmeans = mean(timeseries,1);
        timeseries = timeseries - onesmatrix*colmeans;
        % calculate variance accounted for by the estimated hdr
        uv{j,k} = sum((timeseries-d.scm*ehdr).^2);
        tv{j,k} = sum(timeseries.^2);
    end
end
disppercent(inf);
% reshape matrix. 
for j = yvals
  for k = slices
    % now reshape into a matrix
    d.unexplainedVariance(:,j,k) = uv{j,k};
    d.totalVariance(:,j,k) = tv{j,k};
  end
end

%%%%%%%%%%%%%%%%%%%
% calcleastsqestimate
%
%      usage: d = calcleastsqestimate(d)
%         by: farshad moradi
%       date: 02/05/07
%    purpose: estimates hrf
%
%%%%%%%%%%%%%%%%%%%
function d = calcleastsqestimate(d)

% init some variables
ehdr=[];

% precalculate the normal equation (this dramatically speeds up things)
% Use the full design
n=zeros(1, 1+length(d.scms));
scm = [];
for i=1:length(d.scms),
  n(i+1)=size(d.scms{i}, 1);
  scm = [scm; d.scms{i}];
end;
n = cumsum(n);
precalcmatrix = pinv(scm);
precalcmatrix = precalcmatrix(:, n(d.scanNum)+1:n(d.scanNum+1));
% no roi, so just do whole volume
slices = 1:d.dim(3);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);
  
% preallocate memory
d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),size(precalcmatrix,1));

% display string
disppercent(-inf,'Calculating hdr');
% cycle through images calculating the estimated hdr and r^s of the 
% estimate.
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  disppercent(max((j-min(yvals))/yvaln,0.1));
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % get hdr for the each voxel
    ehdr{j,k} = precalcmatrix*timeseries;
  end
end
disppercent(inf);
% reshape matrix. this also seems the fastest way to do things. we
% could have made a matrix in the above code and then reshaped here
% but the reallocs needed to continually add space to the matrix
% seems to be slower than the loops needed here to reconstruct
% the matrix from the {} arrays.
disppercent(-inf,'Reshaping matrices');
for i = xvals
    disppercent((i-min(xvals))/xvaln);
    for j = yvals
        for k = slices
            % now reshape into a matrix
            d.ehdr(i,j,k,:) = ehdr{j,k}(:,i);
        end
    end
end

% display time took
disppercent(inf);

