% pRFModel.m
%
%        $Id:$ 
%      usage: pRFModel()
%         by: justin gardner
%       date: 11/21/15
%    purpose: 
%
function retval = pRFModel()

% TODO:
% Overlay is not displaying?
% stim image flipping
% hemodynamic response pull out from this code (check)
% dispfit
% betaEachScan and beta parameters
% displaying of stimulus from dialog
% pRFFit work again with click on voxel
% pull out of filtering (check)
% Check against original 

% register the two models - one which just fits
% a guasssian receptive field 
pRFRegisterModel('gaussian',@initFit,@getStim,@initScan,@initVoxel,@getModel,@endScan,@endFit);

%%%%%%%%%%%%%%%%%
%    initFit    %
%%%%%%%%%%%%%%%%%
function [fitParams overlays] = initFit(fitParams,defaultOverlay)

% function does initial processing. Should set up
% overlays that will get filled by the analysis and
% set some info in fitParams (see below)

% init the overlays
overlays = initOverlays(defaultOverlay);

% set fields that specify what the params are
% number of parameters
fitParams.nParams = 3;
% parameter names
fitParams.paramNames = {'x','y','rfWidth'};
% which param goes with each overlay (needs to be same length as
% the number of overlays add above in initOverlays)
fitParams.overlayParamNum = [1 2 3];

%%%%%%%%%%%%%%%%%%
%    initScan    %
%%%%%%%%%%%%%%%%%%
function [fitParams prefitParams] = initScan(fitParams)

% this function is used to set a nxk list of prefitParams
% where n is arbitrary and k is the number of params. This
% params list will then be used to compute model responses
% at each of these paramter values and the time series will
% be correlated with the output. The best match will then
% be used as the initial parameters for the nonlinear fit.

% you can also do any other setup needed for each scan.

% compute prefit
[fitParams prefitParams] = getPrefit(fitParams);

%%%%%%%%%%%%%%%%%%%
%    initVoxel    %
%%%%%%%%%%%%%%%%%%%
function [fitParams tSeries initParams] = initVoxel(fitParams,tSeries,initParams)

% optional: function does initial processing for each voxle.
% Any preprocessing of tSeries can happen here.

% also, if not using prefit, then set the initParams here


%%%%%%%%%%%%%%%%%%
%    getModel    %
%%%%%%%%%%%%%%%%%%
function [rfModel modelResponse] = getModel(params,fitParams)

% gets the model time series. rfModel is the model (i.e. receptive
% field of the voxel)
%
% modelResponse is the expected model response time series
% (usually of length n = num time points). This should *NOT* 
% be convolved with a hemodynamic response. You do *NOT* need
% to handle junk frames or filtering.

% convert parameter array into a parameter strucutre
p = getFitParams(params,fitParams);

% compute an RF
rfModel = getRFModel(p,fitParams);

% init model response
modelResponse = [];

% create the model for each concat
for i = 1:fitParams.concatInfo.n

  % get model response
  nFrames = fitParams.concatInfo.runTransition(i,2)-fitParams.concatInfo.runTransition(i,1)+1;
  thisModelResponse = convolveModelWithStimulus(rfModel,fitParams.stim{i},nFrames);

  % make into a column array
  modelResponse = [modelResponse;thisModelResponse(:)];
end

%%%%%%%%%%%%%%%%%
%    endScan    %
%%%%%%%%%%%%%%%%%
function [overlayParams d] = endScan(fitParams,overlayParams,d)

% optional: for each scan before, do any final
% processing of params before they go into the 
% overlay. Also, can change d structure if necessary

% change from cartesian to polar angle
[overlayParams(:,1) overlayParams(:,2)] = cart2pol(overlayParams(:,1),overlayParams(:,2));

%%%%%%%%%%%%%%%%
%    endFit    %
%%%%%%%%%%%%%%%%
function a = endFit(fitParams,a)

% optional: do any final processing of analysis before it gets saved

%%%%%%%%%%%%%%%%%%%%%%
%    initOverlays    %
%%%%%%%%%%%%%%%%%%%%%%
function overlays = initOverlays(defaultOverlay)

% create the parameters for the polarAngle overlay
overlays.polarAngle = defaultOverlay;
overlays.polarAngle.name = 'polarAngle';
overlays.polarAngle.range = [-pi pi];
overlays.polarAngle.clip = [-pi pi];
overlays.polarAngle.colormapType = 'normal';
overlays.polarAngle.colormap = hsv(256);

% create the parameters for the eccentricity overlay
overlays.eccentricity = defaultOverlay;
overlays.eccentricity.name = 'eccentricity';
overlays.eccentricity.range = [0 15];
overlays.eccentricity.clip = [0 inf];
overlays.eccentricity.colormapType = 'normal';
overlays.eccentricity.colormap = copper(256);

% create the paramteres for the rfHalfWidth overlay
overlays.rfHalfWidth = defaultOverlay;
overlays.rfHalfWidth.name = 'rfHalfWidth';
overlays.rfHalfWidth.range = [0 15];
overlays.rfHalfWidth.clip = [0 inf];
overlays.rfHalfWidth.colormapType = 'normal';
overlays.rfHalfWidth.colormap = pink(256);

%%%%%%%%%%%%%%%%%
%    getStim    %
%%%%%%%%%%%%%%%%%
function [fitParams stim] = getStim(fitParams,v,scanNum)

% get stimfile
stimfile = viewGet(v,'stimfile',scanNum);

% get volume to trigger ratio
volTrigRatio = viewGet(v,'auxParam','volTrigRatio',scanNum);

% display that we are computing
disp(sprintf('(pRFFit) Computing stim image'));

% create stim image
stim = pRFGetStimImageFromStimfile(stimfile,'volTrigRatio',volTrigRatio,'xFlip',fitParams.xFlipStimulus,'yFlip',fitParams.yFlipStimulus,'timeShift',fitParams.timeShiftStimulus,'verbose',fitParams.verbose,'saveStimImage',fitParams.saveStimImage,'recomputeStimImage',fitParams.recomputeStimImageAndPrefit);
  
% check for averages
stim = checkStimForAverages(v,scanNum,viewGet(v,'curGroup'),stim,fitParams.concatInfo,fitParams.stimImageDiffTolerance);
if isempty(stim),return,end
  
% make into cell array
stim = cellArray(stim);

% get stimulus x,y and t
fitParams.stimX = stim{1}.x;
fitParams.stimY = stim{1}.y;
fitParams.stimT = stim{1}.t;

% set stimulus extents
fitParams.stimExtents(1) = min(fitParams.stimX(:));
fitParams.stimExtents(3) = max(fitParams.stimX(:));
fitParams.stimExtents(2) = min(fitParams.stimY(:));
fitParams.stimExtents(4) = max(fitParams.stimY(:));
fitParams.stimWidth = fitParams.stimExtents(3)-fitParams.stimExtents(1);
fitParams.stimHeight = fitParams.stimExtents(4)-fitParams.stimExtents(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimForAverages    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim ignoreMismatchStimfiles] = checkStimForAverages(v,scanNum,groupNum,stim,concatInfo,stimImageDiffTolerance)

ignoreMismatchStimfiles = false;  

% this function will check for some bad casses (like concat of concats etc)
% it will also check that all the component scans of an average have the
% same stim image and warn if they do not. It will then replace the stim cell
% array for the average with a single stim file, so that processing
% can continue as normal for pRFFit

% if not a cell, then ok, return
if ~iscell(stim),return,end

% first check for bad shiftList or refverseLIst
p = viewGet(v,'params',scanNum,groupNum);
if isfield(p,'params') && isfield(p.params,'shiftList') && any(p.params.shiftList~=0)
  disp(sprintf('(pRFFit) Component scan %s:%i has a shiftList that is non-zero (%s). pRFFit does not handle non-zero shifts in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
  keyboard
end
if isfield(p,'params') && isfield(p.params,'reverseList') && any(p.params.reverseList~=0)
  disp(sprintf('(pRFFit) Component scan %s:%i has a reverseList that is non-zero (%s). pRFFit does not handle time-reversed time series in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
  keyboard
end

% if is a cell, check to see if this is a concat or not
if ~isempty(concatInfo) && (concatInfo.isConcat)
  % this is a concat, so check each one of the elements
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  for i = 1:length(stim)
    % get concatInfo for original scan
    concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
    if ~isempty(concatInfo)
      disp(sprintf('(pRFFit:checkStimForAverages) Detected concatenation of concatenations. pRFFit not implemented yet to handle this'));
      stim = [];
      keyboard
      return;
    end
    % check this next scan
    [stim{i} ignoreMismatchStimfiles] = checkStimForAverages(v,originalScanNum(i),originalGroupNum(i),stim{i},concatInfo,stimImageDiffTolerance);
    % if user has accepted all then set stimImageDiffTOlerance to infinity
    if isinf(ignoreMismatchStimfiles),stimImageDiffTolerance = inf;end
    if isempty(stim{i}),stim = [];return,end
  end
else
  % this for orignals
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  % if it is an original than check each element
  if ~isempty(originalScanNum)
    % check that this is not an average of a concat
    for i = 1:length(stim)
      % get concatInfo for original scan
      concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
      if ~isempty(concatInfo)
	disp(sprintf('(pRFFit:checkStimForAverages) Detected average of a concatenations. pRFFit not implemented yet to handle this'));
	keyboard
	stim = [];
	return;
      end
      % see if it is an average of an average
      originalOfOriginalScanNum = viewGet(v,'originalScanNum',originalScanNum(i),originalGroupNum(i));
      if length(originalOfOriginalScanNum) > 1
	disp(sprintf('(pRFFit:checkStimForAverages) Detected average of an average. pRFFit not implemented yet to handle this'));
	keyboard
	stim = [];
	return;
      end
    end
    % ok, not an average of a concatenation/average so check all the stim files 
    % and warn if there are any inconsistencies
    for i = 1:length(stim)
      if ~isequalwithequalnans(stim{1}.im,stim{i}.im)    
	dispHeader
	disp(sprintf('(pRFFit:checkStimForAverages) !!! Average for %s:%i component scan %i does not match stimulus for other scans. If you wish to continue then this will use the stimfile associated with the first scan in the average !!!',viewGet(v,'groupName',groupNum),scanNum,originalScanNum(i)));
	% display which volumes are different
	diffVols = [];
	for iVol = 1:min(size(stim{1}.im,3),size(stim{i}.im,3))
	  if ~isequalwithequalnans(stim{1}.im(:,:,iVol),stim{i}.im(:,:,iVol))
	    diffVols(end+1) = iVol;
	  end
	end
	disp(sprintf('(pRFFit) Stimulus files are different at %i of %i vols (%0.1f%%): %s',length(diffVols),size(stim{1}.im,3),100*length(diffVols)/size(stim{1}.im,3),num2str(diffVols)));
	if 100*(length(diffVols)/size(stim{1}.im,3)) < stimImageDiffTolerance
	  disp(sprintf('(pRFFit) This could be for minor timing inconsistencies, so igorning. Set stimImageDiffTolerance lower if you want to stop the code when this happens'));
	else
	  % ask user if they want to continue (only if there is a difference of more than 10 vols	  
	  ignoreMismatchStimfiles = askuser('Do you wish to continue',1);
	  if ~ignoreMismatchStimfiles
	    stim = [];
	    return;
	  end
	end
	dispHeader
      end
    end
    % if we passed the above, this is an average of identical
    % scans, so just keep the first stim image since they are all the same
    stim = stim{1};
  end
end

%%%%%%%%%%%%%%%%%%%
%    getPrefit    %
%%%%%%%%%%%%%%%%%%%
function [fitParams prefitParams] = getPrefit(fitParams)

% set the values over which to first prefit
% the best of these parameters will then be used 
% to init the non-linear optimization. 
% Here values are expressed as percentile of screen size
if fitParams.quickPrefit
  if fitParams.verbose,disp(sprintf('(pRFFit) Doing quick prefit'));end
  % make sure here that x and y points go through 0 symmetrically
  [prefitX prefitY prefitRFHalfWidth] = ndgrid(-0.375:0.125:0.375,-0.375:0.125:0.375,[0.025 0.05 0.15 0.4]);
else
  [prefitX prefitY prefitRFHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
end

% convert the x/y and width parameters into sizes
% on the actual screen
prefitX = prefitX(:)*fitParams.stimWidth;
prefitY = prefitY(:)*fitParams.stimHeight;
prefitRFHalfWidth = prefitRFHalfWidth*max(fitParams.stimWidth,fitParams.stimHeight);

% set prefit params
prefitParams = [prefitX(:) prefitY(:) prefitRFHalfWidth(:)];

%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)

p.rfType = fitParams.rfType;

switch (fitParams.rfType)
  case 'gaussian'
    p.x = params(1);
    p.y = params(2);
    p.std = params(3);
otherwise 
    disp(sprintf('(pRFFit) Unknown rfType %s',rfType));
end

%%%%%%%%%%%%%%%%%%%%
%%   getRFModel   %%
%%%%%%%%%%%%%%%%%%%%
function rfModel = getRFModel(params,fitParams)

rfModel = [];

% now gernerate the rfModel
if any(strcmp(fitParams.rfType,{'gaussian','gaussian-hdr'}))
  rfModel = makeRFGaussian(params,fitParams);
else
  disp(sprintf('(pRFFit:getRFModel) Unknown rfType: %s',fitParams.rfType));
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeRFGaussian   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFGaussian(params,fitParams)

% compute rf
rfModel = exp(-(((fitParams.stimX-params.x).^2)/(2*(params.std^2))+((fitParams.stimY-params.y).^2)/(2*(params.std^2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = convolveModelWithStimulus(rfModel,stim,nFrames)

% get number of frames
nStimFrames = size(stim.im,3);

% preallocate memory
modelResponse = zeros(1,nFrames);

for frameNum = 1:nStimFrames
  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
  modelResponse(frameNum) = sum(sum(rfModel.*stim.im(:,:,frameNum)));
end

