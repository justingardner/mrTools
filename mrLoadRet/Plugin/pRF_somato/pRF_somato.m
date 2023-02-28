% pRF_somato.m
%
%
%      usage: pRF_somato(v,params,varargin)
%         by: only slightly modified from pRF.m by justin gardner
%       date:
%    purpose: compute pRF analysis on MLR data
%
%             if you just want a default parameter structure you
%             can do:
%
%             v = newView;
%             [v params] = pRF_somato(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = pRF_somato(v,[],'justGetParams=1');
%
function [v d] = pRF_somato(v,params,varargin)

% check arguments
if nargin < 1
  help pRF_somato
  return
end

d = [];
% a version number in case we make major changes
pRFVersion = 1;

% params defaults to empty
if nargin < 2,params =[];end

% other arguments
justGetParams=[];defaultParams=[];scanList=[];
groupNum=[];
getArgs(varargin,{'justGetParams=0','defaultParams=0','scanList=[]','groupNum=[]', 'crossVal=[]'});

% first get parameters
if isempty(params)
  % get group
  if isempty(groupNum),groupNum = viewGet(v,'curGroup');end
  % put up the gui
  params = pRF_somatoGUI('v',v,'groupNum',groupNum,'defaultParams',defaultParams,'scanList',scanList);
end

% just return parameters
if justGetParams,d = params;return,end

% Reconcile params with current status of group and ensure that it has
% the required fields.
params = defaultReconcileParams([],params);

% Abort if params empty
if isempty(params),return,end

% check the params
%params = checkPRFparams(params);

% set the group
v = viewSet(v,'curGroup',params.groupName);

% create the parameters for the r2 overlay

% mod = 'somato';
% overlayNames = getMetaData(v,params,mod,'overlayNames');
% theOverlays = getMetaData(v,params,mod,'theOverlays');


dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'pRF_somato';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(v,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(v,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
%colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'normal';
r2.interrogator = 'myOverlayStats';
r2.mergeFunction = 'pRFMergeParams';

% at this point we need to decide on which parameters we want to estimate
% from data

% for pRF_somato e.g.
%       prefDigit (1, 2, 3)
%       rfHalfWidth...?
%       etc.
% 
% create the parameters for the prefDigit overlay
prefDigit = r2;
prefDigit.name = 'prefDigit';

prefDigit.range = [1 4];
prefDigit.clip = [1 4];
prefDigit.colormapType = 'normal';
%prefDigit.colormap = rainbow_colors(4); %This colour map and 1-3 range look much better when delineating 3 fingers. Change for more fingers!
prefDigit.colormap = digits(4);

prefPD = r2;
prefPD.name = 'prefPD';
prefPD.range = [1 4];
prefPD.clip = [1 4];
prefPD.colormapType = 'normal';
%prefPD.colormap = rainbow_colors(4);
prefPD.colormap = digits(256);


% create the paramteres for the rfHalfWidth overlay
% deal with the sigma.

if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
    rfHalfWidthX = r2;
    rfHalfWidthX.name = 'rfHalfWidthX';
    rfHalfWidthX.range = [0 5];
    rfHalfWidthX.clip = [0 5];
    rfHalfWidthX.colormapType = 'normal';
    rfHalfWidthX.colormap = pink(256);
    
    rfHalfWidthY = r2;
    rfHalfWidthY.name = 'rfHalfWidthY';
    rfHalfWidthY.range = [0 5];
    rfHalfWidthY.clip = [0 5];
    rfHalfWidthY.colormapType = 'normal';
    rfHalfWidthY.colormap = pink(256);
else
    rfHalfWidth = r2;
    rfHalfWidth.name = 'rfHalfWidth';
    rfHalfWidth.range = [0 5];
    rfHalfWidth.clip = [0 5];
    rfHalfWidth.colormapType = 'normal';
    rfHalfWidth.colormap = pink(256);
    
end




% % consider creating other parameters for the somatosensory
% % Maybe get the haemodynamic delay?

% hrfDelay = r2;
% hrfDelay.name = 'hrfDelay';
% hrfDelay.range = [0 5];
% hrfDelay.clip = [0 inf];
% hrfDelay.colormapType = 'normal';
% hrfDelay.colormap = hot(256);

% % make space to keep rawParams in d
%rawParametersFromFit = cell(1,viewGet(v,'nScans'));


% get number of workers
nProcessors = mlrNumWorkers;

% code snippet for clearing precomputed prefit
%global gpRFFitStimImage;gpRFFitStimImage = [];

dispHeader
disp(sprintf('(pRF_somato) Running on scans %s:%s (restrict %s)',params.groupName,num2str(params.scanNum,'%i '),params.restrict ));

for scanNum = params.scanNum
  % see how long it took
  tic;

  % get voxels that we are restricted to
  [x y z] = getVoxelRestriction(v,params,scanNum);
  if isempty(x)
    disp(sprintf('(pRF_somato) No voxels to analyze with current restriction'));
    return
  end

  % get total number of voxels
  n = length(x);

  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum);

  % init overlays
  r2.data{scanNum} = nan(scanDims);
  prefDigit.data{scanNum} = nan(scanDims);
  prefPD.data{scanNum} = nan(scanDims);
  if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
      rfHalfWidthX.data{scanNum} = nan(scanDims);
      rfHalfWidthY.data{scanNum} = nan(scanDims);
  else
      rfHalfWidth.data{scanNum} = nan(scanDims);
  end
  
%   for iOverlay = 1:numel(overlayNames)
%       %
%       theOverlays{iOverlay}.data{scanNum} = nan(scanDims);
%   end


  % default all variables that will be returned
  % by pRFFIt, so that we can call it the
  % second time and save some time
  concatInfo = [];
  stim = [];

  % save pRF parameters
  pRFAnal.d{scanNum}.ver = pRFVersion;
  pRFAnal.d{scanNum}.linearCoords = [];
  pRFAnal.d{scanNum}.params = [];

  % get some information from pRFFit that will be used again in
  % the fits, including concatInfo, stim, prefit, etc.
  fit = pRF_somatoFit(v,scanNum,[],[],[],'fitTypeParams',params.pRFFit,'returnPrefit',true);
  if isempty(fit),return,end
  
  % here we now how many fit params there will be, so make space.
  
  
  stim = fit.stim;
  pRFAnal.d{scanNum}.stim = cellArray(stim);
  pRFAnal.d{scanNum}.stimX = fit.stimX;
  pRFAnal.d{scanNum}.stimY = fit.stimY;
  pRFAnal.d{scanNum}.stimT = fit.stimT;
  concatInfo = fit.concatInfo;
  pRFAnal.d{scanNum}.concatInfo = fit.concatInfo;
  prefit = fit.prefit;
  paramsInfo = fit.paramsInfo;
  pRFAnal.d{scanNum}.paramsInfo = paramsInfo;
  % grab all these fields and stick them onto a structure called paramsInfo
  % preallocate some space
  % fudge this for now
    tf = strcmpi('gaussian-tips', params.pRFFit.rfType);
    if tf == 1
        rawParams = nan(fit.nParams+1,n);
    else 
        rawParams = nan(fit.nParams,n);
    end
  %thisWeights = nan(numel(fit.stimX),n);
  rawParams = nan(fit.nParams,n);
  thisRawParamsCoords = nan(3,n);
  r = nan(n,fit.concatInfo.n);
  thisr2 = nan(1,n);
  
  if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
      thisRfHalfWidthX = nan(1,n);
      thisRfHalfWidthY = nan(1,n);
      
  else
      thisRfHalfWidth = nan(1,n);
  end
  
  thisResid = nan(numel(concatInfo.whichScan),n); % this may break if using Nelder-Mead
  thistSeries = nan(numel(concatInfo.whichScan),n);
  thismodelResponse = nan(numel(concatInfo.whichScan),n); 
  
  %thisData = nan(numel(overlayNames), n);

  % get some info about the scan to pass in (which prevents
  % pRFFit from calling viewGet - which is problematic for distributed computing
  framePeriod = viewGet(v,'framePeriod');
  junkFrames = viewGet(v,'junkFrames',scanNum);

  % compute pRF for each voxel in the restriction
  if params.pRFFit.prefitOnly,algorithm='prefit-only';else algorithm=params.pRFFit.algorithm;end

  % disp info about fitting
  dispHeader;
  disp(sprintf('(pRF_somato) Scan %s:%i (restrict %s) running on %i processor(s)',params.groupName,scanNum,params.restrict,nProcessors));
  disp(sprintf('(pRF_somato) Computing %s fits using %s for %i voxels',params.pRFFit.rfType,algorithm,n));
  dispHeader;

  % this is a bit arbitrary but is the number of voxels to read in at a time.
  % should probably be either calculated based on memory demands or a
  % user settings. The bigger the number the less overhead and will run faster
  % but consume more memory. The overhead is not terribly significant though
  % as tested on my machine - maybe a few percent faster with full n, but
  % on many machines without enough memory that will crash it so keeping
  % this preliminary value in for now.
  blockSize = n;
  tic;
  % break into blocks of voxels to go easy on memory
  % if blockSize = n then this just does on block at a time.
  for blockStart = 1:blockSize:n

    % display information about what we are doing
    % get blockEnd
    blockEnd = min(blockStart + blockSize-1,n);
    blockSize = blockEnd-blockStart+1;

    % load ROI
    loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
    loadROI.coords(1,1:blockSize) = x(blockStart:blockEnd);
    loadROI.coords(2,1:blockSize) = y(blockStart:blockEnd);
    loadROI.coords(3,1:blockSize) = z(blockStart:blockEnd);
    % load all time series for block, we do this to pass into pRFFit. Generally
    % the purpose here is that if we run on distributed computing, we
    % can't load each voxel's time series one at a time. If this is
    % too large for memory then you can comment this out and not
    % pass it into pRFFit and pRFFit will load the tSeries itself
    loadROI = loadROITSeries(v,loadROI,scanNum,params.groupName);
    % reorder x,y,z coordinates since they can get scrambled in loadROITSeries
    
    blockEnd = size(loadROI.scanCoords,2); % HACK TO STOP NANS
    blockSize = blockEnd;
    n = blockEnd;
    
    x(blockStart:blockEnd) = loadROI.scanCoords(1,1:blockSize);
    y(blockStart:blockEnd) = loadROI.scanCoords(2,1:blockSize);
    z(blockStart:blockEnd) = loadROI.scanCoords(3,1:blockSize);
    % keep the linear coords
    pRFAnal.d{scanNum}.linearCoords = [pRFAnal.d{scanNum}.linearCoords sub2ind(scanDims,x(blockStart:blockEnd),y(blockStart:blockEnd),z(blockStart:blockEnd))];

    if blockStart ~= 1
      % display time update
      dispHeader(sprintf('(pRF_somato) %0.1f%% done in %s (Estimated time remaining: %s)',100*blockStart/n,mlrDispElapsedTime(toc),mlrDispElapsedTime((toc*n/blockStart) - toc)));
    end
    
    
    %import hrfprf code from pRF.m
    if params.pRFFit.HRFpRF == 1
        disp('Give me your HRFs. Remember, these should be outputted from prfhrfRefit and then ideally from deconvRealDataWiener')
        myfilename_hrf = uigetfile;
        thehrfs = load(myfilename_hrf);
        
        myVar = thehrfs.hrf_struct.yf;
   
    end
    
    
    if params.pRFFit.HRFpRF == 1
        
        %sliceFix = 128.*128.*12;
        %thehrfs.hrf_struct.volumeIndices = thehrfs.hrf_struct.volumeIndices + sliceFix;
        
        for ii = blockStart:blockEnd
            myVoxel = find(thehrfs.hrf_struct.volumeIndices == sub2ind(scanDims,x(ii),y(ii),z(ii)));
            if isempty(myVoxel)
                fprintf('\ncaught an empty, x %d y %d z %d, idx %f\n', x(ii), y(ii), z(ii), myVoxel);
                
                fit = [];
            elseif myVoxel > length(thehrfs.hrf_struct.yf)
                disp('caught one')
                fit = [];
            else
                
                fit = pRF_somatoFit(v,scanNum,x(ii),y(ii),z(ii),'stim',stim,'concatInfo',concatInfo,...
                    'prefit',prefit,'fitTypeParams',params.pRFFit,'dispIndex',ii,'dispN',n,...
                    'tSeries',loadROI.tSeries(ii-blockStart+1,:)','framePeriod',framePeriod,'junkFrames',junkFrames,...
                    'paramsInfo',paramsInfo, 'hrfprf', myVar(:,myVoxel));
            end
            if ~isempty(fit)
                thisr2(ii) = fit.r2;
                thisPrefDigit(ii) = fit.prefDigit;
                thisPrefPD(ii) = fit.prefPD;
                thisRfHalfWidth(ii) = fit.std;
                thisResid(:,ii) = fit.residual;
                thistSeries(:,ii) = fit.tSeries;
                thismodelResponse(:,ii) = fit.modelResponse;

%                 
%                 tempVar = zeros(length(overlayNames),1);
%                 for iOverlay = 1:numel(overlayNames)
%                     
%                     test = strcmpi(fieldnames(fit), overlayNames(iOverlay) );
%                     %pos = find(test==1);
%                     bla = struct2cell(fit);
%                     val = cell2mat(bla(test==1));
%                     % this is temporary, gets overwritten each time
%                     tempVar(iOverlay,1) = val;
%                 end
%                 % now put the values for this voxel into some sort of order :)
%                 thisData(:,ii) = tempVar;
                
                % keep parameters
                rawParams(:,ii) = fit.params(:);
                r(ii,:) = fit.r;
                %thisr2(ii) = fit.r2;
                thisRawParamsCoords(:,ii) = [x(ii) y(ii) z(ii)];
                %myrawHrfs(:,ii) = fit.myhrf.hrf; %save out prfs hrfs
            end
        end
        
    else
        
        for ii = blockStart:blockEnd
            
            fit = pRF_somatoFit(v,scanNum,x(ii),y(ii),z(ii),'stim',stim,'concatInfo',concatInfo, ...
                'prefit',prefit,'fitTypeParams',params.pRFFit,'dispIndex',ii,'dispN',n,...
                'tSeries',loadROI.tSeries(ii-blockStart+1,:)','framePeriod',framePeriod,...
                'junkFrames',junkFrames,'paramsInfo',paramsInfo);
            %fit = pRF_somatoFit(v,scanNum,x(ii),y(ii),z(ii),'stim',stim,'concatInfo',concatInfo,'prefit',prefit,'fitTypeParams',params.pRFFit,'dispIndex',ii,'dispN',n,'tSeries',loadROI.tSeries(ii-blockStart+1,:)','framePeriod',framePeriod,'junkFrames',junkFrames,'paramsInfo',paramsInfo, 'crossVal', myVar(:,ii));
            
            
            
            if ~isempty(fit)
                %                 tempVar = zeros(length(overlayNames),1);
                %                 for iOverlay = 1:numel(overlayNames)
                %
                %                     test = strcmpi(fieldnames(fit), overlayNames(iOverlay) );
                %                     %pos = find(test==1);
                %                     bla = struct2cell(fit);
                %                     val = cell2mat(bla(test==1));
                %                     % this is temporary, gets overwritten each time
                %                     tempVar(iOverlay,1) = val;
                %                 end
                %                 % now put the values for this voxel into some sort of order :)
                %                 thisData(:,ii) = tempVar;
                
                thisr2(ii) = fit.r2;
                thisPrefDigit(ii) = fit.prefDigit;
                thisPrefPD(ii) = fit.prefPD;
                if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
                    thisRfHalfWidthX(ii) = fit.rfHalfWidthX;
                    thisRfHalfWidthY(ii) = fit.rfHalfWidthY;
                else
                    thisRfHalfWidth(ii) = fit.rfHalfWidth;
                end
                
                
                % keep parameters
                rawParams(:,ii) = fit.params(:);
                r(ii,:) = fit.r;
                %thisr2(ii) = fit.r2;
                thisRawParamsCoords(:,ii) = [x(ii) y(ii) z(ii)];
                if ~strcmpi(algorithm,'nelder-mead')
                    thisResid(:,ii) = fit.residual;
                end
                thistSeries(:,ii) = fit.tSeries;
                thismodelResponse(:,ii) = fit.modelResponse;
                %myrawHrfs(:,ii) = fit.myhrf.hrf; %save out prfs hrfs
            end
            
        end
        
    end
      %% debugging, show rf model each time
%       if strcmpi('gaussian', fit.rfType)
%           tt = exp(-(((pRFAnal.d{2}.stimX-fit.x).^2)/(2*(fit.std^2))+((pRFAnal.d{2}.stimY-fit.y).^2)/(2*(fit.std^2))));
%           if ii == 1
%               figure
%               plot(pRFAnal.d{2}.stimX, tt)
%           elseif ii > 1
%               hold on
%               plot(pRFAnal.d{2}.stimX, tt)
%           end
%           
%       elseif strcmpi('gaussian-tips', fit.rfType)
%           X = pRFAnal.d{2}.stimX;
%           pone = [fit.amp fit.meanOne fit.std 0];
%           Z = gauss(pone,X);
%           if ii == 1
%               figure
%               plot(X,Z)
%           elseif ii > 1
%               hold on
%               plot(X,Z)
%           end
%           
%           
%       end
      %% 
    % set overlays and info for d
    for ii = 1:n
              r2.data{scanNum}(x(ii),y(ii),z(ii)) = thisr2(ii);
              prefDigit.data{scanNum}(x(ii),y(ii),z(ii)) = thisPrefDigit(ii);
              prefPD.data{scanNum}(x(ii),y(ii),z(ii)) = thisPrefPD(ii);
              if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
                  rfHalfWidthX.data{scanNum}(x(ii),y(ii),z(ii)) = thisRfHalfWidthX(ii);
                  rfHalfWidthY.data{scanNum}(x(ii),y(ii),z(ii)) = thisRfHalfWidthY(ii);
              else
                  rfHalfWidth.data{scanNum}(x(ii),y(ii),z(ii)) = thisRfHalfWidth(ii);
              end
             
%         for iOverlay = 1:length(overlayNames)
%             theOverlays{iOverlay}.data{scanNum}(x(ii),y(ii),z(ii)) = thisData(iOverlay,ii);
%         end
    end
  end
  % display time update
  dispHeader;
  disp(sprintf('(pRF_somato) Fitting %i voxels took %s.',n,mlrDispElapsedTime(toc)));
  dispHeader;
  
  
  pRFAnal.d{scanNum}.time = toc; %speed testing
  
  pRFAnal.d{scanNum}.params = rawParams;
  pRFAnal.d{scanNum}.r = r;
  pRFAnal.d{scanNum}.r2 = thisr2;
  pRFAnal.d{scanNum}.maxr2 = max(thisr2); % saves out maximum voxel peak (for curiosity)
  pRFAnal.d{scanNum}.rawCoords = thisRawParamsCoords; % this is where we save it, so we can access it via the d structure
  %pRFAnal.d{scanNum}.weights = thisWeights;
  if ~strcmpi(algorithm,'nelder-mead')
      pRFAnal.d{scanNum}.myresid = thisResid;
  end
  pRFAnal.d{scanNum}.mytSeries = thistSeries;
  pRFAnal.d{scanNum}.mymodelResp = thismodelResponse;
  

  iScan = find(params.scanNum == scanNum);
  thisParams.scanNum = params.scanNum(iScan);
%   for iOverlay = 1:length(overlayNames)
%       theOverlays{iOverlay}.params{scanNum} = thisParams;
%   end
  
  r2.params{scanNum} = thisParams;
  prefDigit.params{scanNum} = thisParams;
  prefPD.params{scanNum} = thisParams;
  if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
      rfHalfWidthX.params{scanNum} = thisParams;
      rfHalfWidthY.params{scanNum} = thisParams;
      
  else
      rfHalfWidth.params{scanNum} = thisParams;
      
  end
  

  % display how long it took
  disp(sprintf('(pRF_somato) Fitting for %s:%i took in total: %s',params.groupName,scanNum,mlrDispElapsedTime(toc)));
end

% install analysis
pRFAnal.name = params.saveName;
pRFAnal.type = 'pRFAnal';
pRFAnal.groupName = params.groupName;
pRFAnal.function = 'pRF_somato';
pRFAnal.reconcileFunction = 'defaultReconcileParams';
pRFAnal.mergeFunction = 'pRFMergeParams';
pRFAnal.guiFunction = 'pRF_somatoGUI';
pRFAnal.params = params;

if strcmpi(params.pRFFit.rfType,'gaussian-hdr-double')
    pRFAnal.overlays = [r2 prefDigit prefPD rfHalfWidthX rfHalfWidthY];
else
    pRFAnal.overlays = [r2 prefDigit prefPD rfHalfWidth ];
end
%pRFAnal.overlays = [];
% for iOverlay = 1:numel(theOverlays)
%     eval(sprintf('%s = struct(theOverlays{iOverlay});',overlayNames{iOverlay}));    
%     eval(sprintf('pRFAnal.overlays = [pRFAnal.overlays %s];',overlayNames{iOverlay}));
% end

pRFAnal.curOverlay = 1;
pRFAnal.date = date;
v = viewSet(v,'newAnalysis',pRFAnal);

% if we are going to merge, temporarily set overwritePolicy
if isfield(params,'mergeAnalysis') && params.mergeAnalysis
  saveMethod = mrGetPref('overwritePolicy');
  mrSetPref('overwritePolicy','Merge');
end
% Save it
saveAnalysis(v,pRFAnal.name);
% now set policy back
if isfield(params,'mergeAnalysis') && params.mergeAnalysis
  mrSetPref('overwritePolicy',saveMethod);
end

if ~isempty(viewGet(v,'fignum'))
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

%set(viewGet(v,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for ii = 1:length(d)
    pRFAnal.d{ii}.r2 = r2.data{ii};
  end
  % make d strucutre
  if length(pRFAnal.d) == 1
    d = pRFAnal.d{1};
  else
    d = pRFAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getVoxelRestriction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y z] = getVoxelRestriction(v,params,scanNum)

x = [];y = [];z = [];

if strncmp(params.restrict,'Base: ',6)
  % get the base name
  baseName = params.restrict(7:end);
  baseNums = [];
  if strcmp(baseName,'ALL')
    for iBase = 1:viewGet(v,'numBase')
      % if the base is a surface or flat then add to the list
      if any(viewGet(v,'baseType',iBase) == [1 2])
	baseNums(end+1) = iBase;
      end
    end
  else
    baseNums = viewGet(v,'baseNum',baseName);
  end
  % cycle through all bases that we are going to run on
  scanCoords = [];
  for iBase = 1:length(baseNums)
    % get the baseNum
    baseNum = baseNums(iBase);
    if isempty(baseNum)
      disp(sprintf('(pRF_somato) Could not find base to restrict to: %s',params.restrict));
      continue
    end
    % get the base
    base = viewGet(v,'base',baseNum);
    if isempty(base)
      disp(sprintf('(pRF_somato) Could not find base to restrict to: %s',params.restrict));
      return;
    end
    % if flat or surface
    if any(base.type == [1 2])
      % get base coordinates from the coordMap
      for corticalDepth = 0:0.1:1
	if base.type == 1
	  % flat map
	  baseCoords = (base.coordMap.innerCoords + corticalDepth * (base.coordMap.outerCoords-base.coordMap.innerCoords));
	  baseCoords = reshape(baseCoords,prod(size(base.data)),3)';
	else
	  % surface
	  baseCoords = (base.coordMap.innerVtcs + corticalDepth * (base.coordMap.outerVtcs-base.coordMap.innerVtcs))';
	end
	% convert to 4xn array
	baseCoords(4,:) = 1;
	% and convert to scan coordinates
	base2scan = viewGet(v,'base2scan',scanNum,params.groupName,baseNum);
	scanCoords = [scanCoords round(base2scan*baseCoords)];
      end
    end
  end
  % check against scandims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  scanCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
  % remove duplicates and nans
  scanCoords = scanCoords(~isnan(scanCoords));
  scanCoords = unique(scanCoords);
  % convert back to x,y,z coordinates
  [x y z] = ind2sub(scanDims,scanCoords);
elseif strncmp(params.restrict,'ROI: ',5)
  % get the roi name
  roiName = params.restrict(6:end);
  scanCoords = getROICoordinates(v,roiName,scanNum,params.groupName,'straightXform=1');
  if isempty(scanCoords),return,end
  x = scanCoords(1,:);y = scanCoords(2,:);z = scanCoords(3,:);
elseif strncmp(params.restrict,'None',4)
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  [x y z]  = ndgrid(1:scanDims(1),1:scanDims(2),1:scanDims(3));
  x = x(:);y = y(:);z = z(:);
else
  return
end

%check if we have already computed Voxels
if isfield(params,'computedVoxels') && (length(params.computedVoxels)>=scanNum) && ~isempty(params.computedVoxels{scanNum})
  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  % convert x, y, z to linear coords
  linearCoords = sub2ind(scanDims,x,y,z);
  % get new ones
  newLinearCoords = setdiff(linearCoords,params.computedVoxels{scanNum});
  if length(newLinearCoords) ~= length(linearCoords)
    % show what we are doing
    disp(sprintf('(pRF) Dropping %i voxels that have been already computed',length(linearCoords)-length(newLinearCoords)));
    % convert back to x, y, z
    [x y z] = ind2sub(scanDims,newLinearCoords);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%    checkPRFparams    %
%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkPRFparams(params)


% check the pRFFit params
checkFields = {{'stimImageDiffTolerance',5}};
for iFit = 1:length(params.pRFFit)

  % set defaults
  for iField = 1:length(checkFields)
    if ~isfield(params.pRFFit(iFit),checkFields{iField}{1})
      params.pRFFit(iFit).(checkFields{iField}{1}) = checkFields{iField}{2};
    end
  end
end
