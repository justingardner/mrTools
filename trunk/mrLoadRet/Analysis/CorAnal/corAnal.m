function [view params] = corAnal(view,params,varargin)
%
% view = corAnal(view,[params])
% 
% Loops throughs scans and slices, loads corresponding tSeries, computes
% correlation analysis from tSeries, and saves the resulting co, amp, and
% ph to the corAnal.mat file along with the analysis parameters.
%
% Checks to see of co, amp, and ph are already loaded. If so, it uses the
% existing maps and recomputes as specified.
%
% params: Optional initial parameters. Default: user is prompted via
%    corAnalGUI. If a corAnal already exists and is loaded then the
%    corAnalGUI is initialized with the existing parameters. Params must be
%    a structure with all of the following fields.
% groupName: group of scans that will be averaged.
%    Default: current group of view.
% recompute: vector of 1 and 0 specifying which scans to compute.
%    Default: all of the scans.
% ncycles: vector specifying number of cycles per scan.
%    Default: 1
% detrend: cell array of strings as in percentTSeries.
%    Default: 'None'
% spatialnorm: cell array of strings as in percentTSeries.
%    Default: 'None'
% tseriesfile: cell array of strings specifying tseries filenames. Or
%    'any' to allow any match with the tseries files. This is useful so
%    that you can run the analysis on multiple data sets using the same
%    params.
%
%
% Examples:
%
% params = corAnalGUI('groupName','Raw');
% view = newView;
% view = corAnal(view,params);
%
% view = corAnal(view);
%
% To just get parameters
% [view params] = corAnal(view,[],'justGetParams=1','defaultParams=1');
%
% djh, 5/2005, updated to mrLoadRet-4.0
% $Id$	

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help corAnal
  return
end

if ~isview(view)
    help corAnal
    mrErrorDlg('(corAnal) Invalid view.')
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

% If corAnal is loaded, then use it. Otherwise, viewGet returns [];
corAnal = viewGet(view,'corAnal');
co = viewGet(view,'co');
amp = viewGet(view,'amp');
ph = viewGet(view,'ph');
if ~isempty(corAnal)
    oldparams = corAnal.params;
else
    oldparams = [];
end

% if we are just getting default parameters then
% get them by calling the reconcile function
if defaultParams
  params = corAnalReconcileParams(viewGet(view,'groupName'),[]);
end

% Get analysis parameters from corAnalGUI, using co.params if it exists
if ieNotDefined('params')
  if ~isempty(oldparams)
    % Initialize analysis parameters with previous values from loaded
    % coherence map and reconcile with current status of group.
    params = corAnalGUI('groupName',viewGet(view,'groupName'),'params',oldparams);        
  else
    % Initialize analysis parameters with default values
    params = corAnalGUI('groupName',viewGet(view,'groupName'));
  end
else
  % Reconcile params with current status of group and ensure that params
  % has the required fields.
  params = corAnalReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
    mrMsgBox('corAnal cancelled',1);
    return
end

% if just getting params then return
if justGetParams,return,end

% Change group, get nScans
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end
nScans = viewGet(view,'nScans',groupNum);

% If co/amp/ph do not exist (above viewGet calls return [])...
% Intialize structure and initialize data to cell array of length nScans.
if isempty(co)
    co.name = 'co';
    co.function = 'corAnal';
    co.reconcileFunction = 'corAnalReconcileParams';
    co.mergeFunction = 'corAnalMergeParams';
    co.data = cell(1,nScans);
end
if isempty(amp)
    amp.name = 'amp';
    amp.function = 'corAnal';
    amp.reconcileFunction = 'corAnalReconcileParams';
    amp.mergeFunction = 'corAnalMergeParams';
    amp.data = cell(1,nScans);
end
if isempty(ph)
    ph.name = 'ph';
    ph.function = 'corAnal';
    ph.reconcileFunction = 'corAnalReconcileParams';
    ph.mergeFunction = 'corAnalMergeParams';
    ph.data = cell(1,nScans);
end

% Compute it
[co,amp,ph] = computeCorAnal(view,params,co,amp,ph);

% Set params field, merging with previous params
if ~isempty(oldparams)
    params = corAnalMergeParams(groupName,oldparams,params);
end
co.params = params;
amp.params = params;
ph.params = params;
    
% Fill range fields, fixed values for co and ph
co.range = [0 1];
ph.range = [0 2*pi];

% Fill range field for amp
ampMin = realmax;
ampMax = 0;
for scan=1:nScans
    if ~isempty(amp.data{scan})
        ampMin = min([ampMin min(amp.data{scan}(:))]);
        ampMax = max([ampMax max(amp.data{scan}(:))]);
    end
end
if (ampMin <= ampMax)
    amp.range = [ampMin ampMax];
else
    % if amp data is empty, need to make sure min < max
    amp.range = [0 1]
end

% Fill colormap fields, keeping previous if already exists
if ~isfield(co,'colormap')
    co.colormap = jet(256);
end
if ~isfield(amp,'colormap')
    amp.colormap = hot(256);
end
if ~isfield(ph,'colormap')
    ph.colormap = hsv(256);
end

% Set date field
dateString = datestr(now);
co.date = dateString;
amp.date = dateString;
ph.date = dateString;

% Set interrogrator function
co.interrogator = 'corAnalPlot';
amp.interrogator = 'corAnalPlot';
ph.interrogator = 'corAnalPlot';

% Set groupName
co.groupName = params.groupName;
amp.groupName = params.groupName;
ph.groupName = params.groupName;

% Install corAnal in the view
corAnal.name = 'corAnal';  % This can be reset by editAnalysisGUI
corAnal.type = 'corAnal';
corAnal.groupName = params.groupName;
corAnal.function = 'corAnal';
corAnal.reconcileFunction = 'corAnalReconcileParams';
corAnal.mergeFunction = 'corAnalMergeParams';
corAnal.guiFunction = 'corAnalGUI';
corAnal.overlayInterpFunction = 'corAnalInterp';
corAnal.params = params;
corAnal.date = dateString;
view = viewSet(view,'newanalysis',corAnal);
view = viewSet(view,'newoverlay',co);
view = viewSet(view,'newoverlay',amp);
view = viewSet(view,'newoverlay',ph);
if ~isempty(viewGet(view,'fignum'))
  refreshMLRDisplay(viewGet(view,'viewNum'));
end

% Save it
saveAnalysis(view,corAnal.name);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [co,amp,ph] = computeCorAnal(view,params,co,amp,ph)
% Required fields in params: 'recompute','ncycles','detrend','spatialnorm'

% Get scanList fromp params.recompute field
scanList = find(params.recompute(:));

disp('Computing corAnal...');
waitHandle = mrWaitBar(0,'Computing corAnal matrices from the tSeries.  Please wait...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
    scanNum = scanList(scanIndex);
    disp(['Processing scan ', int2str(scanNum),'...']);
    
    % sliceDims: [ydim xdim] for single slice
    % volDims; [ydim xdim nslices] for single scan
    sliceDims = viewGet(view,'sliceDims',scanNum);
    volDims = viewGet(view,'dims',scanNum);
    
    % Initialize data with NaNs
    co.data{scanNum} = NaN*ones(volDims);
    amp.data{scanNum} = NaN*ones(volDims);
    ph.data{scanNum} = NaN*ones(volDims);
    
    nslices = viewGet(view,'nslices',scanNum);
    for sliceNum = 1:nslices
        [coSeries,ampSeries,phSeries] = computeCorAnalSeries(view,scanNum,sliceNum,params);
        switch view.viewType
        case {'Volume'}
            co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
            amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
            ph.data{scanNum}(:,:,sliceNum) = reshape(phSeries,sliceDims);
        case {'Surface'}
            co.data{scanNum} = coSeries;
            amp.data{scanNum} = ampSeries;
            ph.data{scanNum} = phSeries;
        case {'Flat'}
            co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
            amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
            ph.data{scanNum}(:,:,sliceNum) = reshape(phSeries,sliceDims);
        end
    end
    % Update waitbar
    mrWaitBar(scanIndex/length(scanList),waitHandle);
end
warning('on','MATLAB:divideByZero');
mrCloseDlg(waitHandle);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [co, amp, ph] = computeCorAnalSeries(view, scan, slice, params)
% function [co, amp, ph] = ...
%     computeCorAnalSeries(view, scan, slice, params)
%
% Computes correlation analysis for the given scan and slice.
%
% Required fields in params: 'ncycles','detrend','spatialnorm'

% Analysis parameters for this scan
junkframes = viewGet(view,'junkframes',scan);
nframes = viewGet(view,'nframes',scan);
ncycles = params.ncycles(scan);
detrend = params.detrend{scan};
spatialnorm = params.spatialnorm{scan};

% Load tSeries
tSeries = loadTSeries(view, scan, slice);

% Reshape the tSeries
% ATTN: added reshapeTSeries function, since loadTSeries not longer reshapes when it loads -eli
tSeries = reshapeTSeries(tSeries);

% check that junkframes and nframes settings are ok
if size(tSeries,1) < (junkframes+nframes)
  mrErrorDlg(sprintf('(corAnal) Number of junkframes (%i) plus nframes (%i) should not be larger than number of volumes in scan %i',junkframes,nframes,size(tSeries,1)));
end

% Remove junkFrames
tSeries = tSeries(junkframes+1:junkframes+nframes,:);

% Set highpassPeriod
highpassPeriod = round(nframes/ncycles);

% Remove dc, convert to percent, detrend, and spatial normalization
warnState = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:divideByZero');
ptSeries = percentTSeries(tSeries,...
    'detrend', detrend,...
    'highpassPeriod', highpassPeriod,...
    'spatialNormalization', spatialnorm,...
    'subtractMean', 'Yes',...
    'temporalNormalization', 'No');
warning(warnState.state,warnState.identifier);

% Compute Fourier transform
ft = fft(ptSeries);
ft = ft(1:1+fix(size(ft, 1)/2), :);
ampFT = 2*abs(ft)/nframes;

% Compute co and amp (avoiding divide by zero)
amp = ampFT(ncycles+1,:);
co = zeros(size(amp));
sumAmp = sqrt(sum(ampFT.^2));
nonzeroIndices = find(sumAmp >0);
co(nonzeroIndices) = ampFT(ncycles+1,nonzeroIndices) ./ sumAmp(nonzeroIndices);

% Calculate phase:
% 1) add pi/2 so that it is in sine phase.
% 2) minus sign because sin(x-phi) is shifted to the right by phi.
% 3) Add 2pi to any negative values so phases increase from 0 to 2pi.
ph = -(pi/2) - angle(ft(ncycles+1,:));
ph(ph<0) = ph(ph<0)+pi*2;
