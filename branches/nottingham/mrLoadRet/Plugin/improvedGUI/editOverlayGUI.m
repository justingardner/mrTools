% editOverlayGUImrParams.m
%
%        $Id: editOverlayGUImrParams.m 2511 2012-05-11 12:19:53Z julien $	
%      usage: editOverlayGUImrParams(viewNum)
%         by: eli merriam
%       date: 09/04/07
%    purpose: 
%
function retval = editOverlayGUImrParams(viewNum)

% check arguments
  if ~any(nargin == [1])
    help editOverlayGuimrParams
    return
  end
  
  if isview(viewNum)
    thisView = viewNum;
    viewNum = thisView.viewNum;
  else
    thisView = viewGet(viewNum,'view');
  end

  % Get the original overlay
  analysisNum = viewGet(thisView,'currentAnalysis');
  if isempty(analysisNum),mrWarnDlg('(editOverlayGUI) No current analysis');return,end
  overlayNum = viewGet(thisView,'currentOverlay', analysisNum);
  if isempty(overlayNum),mrWarnDlg('(editOverlayGUI) No current overlay');return,end
  if length(overlayNum)>1
    mrWarnDlg('(editOverlayGUI) Not implemented for several overlays');
    %close aEdit Overlay Dialog if it's on (that could be the case if overlay is changed while the edit overlay params dialog is already on)
    global gParams
    if ~isempty(gParams) && strcmp(gParams.figlocstr{1},'mrParamsDialog_Change_overlay_colormap')
      closeHandler;
    end
    return;
  end
  overlay = viewGet(thisView, 'overlay', overlayNum, analysisNum);
%   overlayUsefulRange = viewGet(thisView,'overlayRange', overlayNum, analysisNum);
%   overlayColorRange = viewGet(thisView,'overlayColorRange', overlayNum, analysisNum);
  overlayRange = viewGet(thisView,'overlayRange', overlayNum, analysisNum);
%   overlayClipRange = viewGet(thisView,'overlayClip', overlayNum, analysisNum);
  overlayClip = viewGet(thisView,'overlayClip', overlayNum, analysisNum);
  overlayName = viewGet(thisView, 'overlayName', overlayNum, analysisNum);
  overlayType = viewGet(thisView, 'overlayType', overlayNum,analysisNum);
  alphaOverlay = viewGet(thisView,'alphaOverlay');
  if isempty(alphaOverlay)
    alphaOverlay = 'none';
  end
  alphaOverlayMenu = putOnTopOfList(alphaOverlay,setdiff(['none',viewGet(thisView,'overlayNames')],overlayName));
  alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
  interrogator = viewGet(thisView,'interrogator',overlayNum,analysisNum);
  
  numColors = size(viewGet(thisView,'overlaycmap',overlayNum),1);
  overlayColormapTypeMenu = putOnTopOfList(viewGet(thisView,'overlayctype',overlayNum),...
    {'normal', 'setRangeToMax', 'setRangeToMaxAroundZero', 'setRangeToMaxAcrossSlices', 'setRangeToMaxAcrossSlicesAndScans'});
  
  % colormaps
  colormaps = {'default','hot','hsv','pink','cool','bone','copper','flag','gray','grayCirc','twoCondCmap','twoCondCircCmap','hsvDoubleCmap','cmapExtendedHSV','overlapCmap','redGreenCmap','rygbCmap','bicolorCmap' 'coolCmap','hotColdCmap'};
  altColormaps = viewGet(thisView,'colormaps');
  if ~isempty(altColormaps)
    colormaps = {colormaps{:} altColormaps{:}};
  end
  
  % set up params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'overlayCmap', colormaps,'type=popupmenu','List of possible colormaps'};
  paramsInfo{end+1} = {'userDefinedCmap','','Allows you to call a user defined function to set the overla colormap'};
  paramsInfo{end+1} = {'numColors', numColors, 'first argument to the colormap function'};
  paramsInfo{end+1} = {'numGrays', 0, 'second argument to the colormap function'};
  paramsInfo{end+1} = {'flipColormap', 0, 'type=checkbox', 'check this box to reverse the direction of the colormap'};
  paramsInfo{end+1} = {'shiftColormap', 0, 'incdec=[-16 16]', 'shift the colormap -- this can be useful for retinotopy scans with circular colormaps'}; 
  paramsInfo{end+1} = {'overlayColormapType',overlayColormapTypeMenu , 'type=popupmenu',...
      '''normal'' scales the colormap to the value specified by colormap range; ''setRangeToMax'' scales the colormap to the min amnd max of the displayed overlay slice and ignores the color range (like R2 maps)'};
%   paramsInfo{end+1} = {'overlayColorRange', overlayColorRange, 'The lower and upper bound on the colormap when overlayColormapType=''normal'''};
%   paramsInfo{end+1} = {'overlayClipRange', overlayClipRange, 'callback',{@checkCmapParams,'cliprange'},'passCallbackOutput=1','passValue=1','passParams=1',...
%       'The lower and upper clip points beyond which the overlay is masked. These should be inside the useful range. If clip(1)>clip(2), then values inside the clip range are masked.'};
%   paramsInfo{end+1} = {'overlayUsefulRange', overlayUsefulRange, 'callback',{@checkCmapParams,'usefulrange'},'passCallbackOutput=1','passValue=1','passParams=1',...
%       'The lower and upper bound on the clip slider. These should be lower/higher than the clip values'};
%   paramsInfo{end+1} = {'setUsefulRange', 0, 'type=pushbutton','callback',@mrCmapSetUsefulRange,'callbackArg',viewNum,'buttonString=Set useful range to overlay min/max','passParams=1','passCallbackOutput=0',...
%       'Sets the useful range to the min/max values of this overlay accross scans'};
  paramsInfo{end+1} = {'overlayRange', overlayRange, 'callback',{@checkCmapParams,'range'},'passCallbackOutput=1','passValue=1','passParams=1',...
      'The lower and upper bound on the clip slider. These should be lower/higher than the clip values'};
  paramsInfo{end+1} = {'overlayClip', overlayClip, 'callback',{@checkCmapParams,'clip'},'passCallbackOutput=1','passValue=1','passParams=1',...
      'The lower and upper clip points beyond which the overlay is masked. These should be inside the useful range. If clip(1)>clip(2), then values inside the clip range are masked.'};
  paramsInfo{end+1} = {'overlayType', overlayType, 'The type of overlay (ph, amp, co ...'};
  paramsInfo{end+1} = {'interrogator', interrogator, 'Sets the overlay default interrogator function'};
  paramsInfo{end+1} = {'alphaOverlay', alphaOverlayMenu, 'You can specify the name of another overlay in the analysis to use as an alpha map. For instance, you might want to display one overlay with the alpha set to the r2 or the coherence value.'};
  paramsInfo{end+1} = {'alphaOverlayExponent', alphaOverlayExponent,'incdec=[-0.1 0.1]','If you are using an alphaOverlay, this sets an exponent on the alphaOverlay to pass through. For example, if you just want the value from the overlay to be the alpha value then set this to 1. If you want to have it so that lower values get accentuated (usually this is the case), set the exponent to less than 1, but greater than 0. The alpha values are passed through the function alpha = alpha.^alphaOverlayExponent'};
  paramsInfo{end+1} = {'setManyOverlays', 0, 'type=pushbutton','Set many overlays to have them same settings as the current overlay. This cannot be cancelled.','callback',@mrCmapSetManyOverlaysCallback,'callbackArg',viewNum,'buttonString=Set many overlays'};
%   paramsInfo{end+1} = {'overlayName', overlayName, 'The name for the overlay (e.g., co, am, or ph)'};

  % display dialog
  mrParamsDialog(paramsInfo,'Change overlay colormap','modal=0','callback',@mrCmapCallback,'callbackArg',viewNum,'cancelCallback',{@mrCmapParamsCancel,overlay});

  return;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapSetUsefulRange   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapSetUsefulRange(viewNum,params)

thisView = viewGet(viewNum,'view');

minOverlayData = floor(double(viewGet(thisView,'minOverlayData'))*1e6)/1e6;
maxOverlayData = ceil(double(viewGet(thisView,'maxOverlayData'))*1e6)/1e6;

if isempty(maxOverlayData) || isempty(minOverlayData)
  mrWarnDlg('(editOverlayGUImrParams) overlay seems to be empty');
% elseif all(abs(params.overlayUsefulRange-[minOverlayData maxOverlayData])<5e-7)
%   mrWarnDlg('(editOverlayGUImrParams) Useful range is already set to min/max');
elseif all(abs(params.overlayRange-[minOverlayData maxOverlayData])<5e-7)
  mrWarnDlg('(editOverlayGUImrParams) Useful range is already set to min/max');
else
%   params.overlayUsefulRange = [minOverlayData maxOverlayData];
  params.overlayRange = [minOverlayData maxOverlayData];
  %make sure clip values are within this range
%   params.overlayClipRange(1) = max(params.overlayClipRange(1),minOverlayData);
%   params.overlayClipRange(2) = min(params.overlayClipRange(2),maxOverlayData);
  params.overlayClip(1) = max(params.overlayClip(1),minOverlayData);
  params.overlayClip(2) = min(params.overlayClip(2),maxOverlayData);
  mrParamsSet(params,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapParamsCancel   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapParamsCancel(oldOverlay,viewNum)

thisView = viewGet(viewNum,'view');
analysisNum = viewGet(thisView,'currentAnalysis');
overlayNum = viewGet(thisView,'currentOverlay',analysisNum);
currentOverlay = viewGet(thisView, 'overlay', overlayNum, analysisNum);

%iff the overlay has changed, put the old overlay params back
if ~isequalwithequalnans(oldOverlay,currentOverlay)
  disppercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
  % set the new overlay
  thisView = viewSet(thisView,'newOverlay', oldOverlay);
  % and refresh
  refreshMLRDisplay(thisView.viewNum);
  disppercent(inf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  checkCmapParams   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = checkCmapParams(params,value,whichParam)

switch(whichParam)
  case 'cliprange'
    %we don't know what value we got so we need to find out first
%     indexInArray = find(params.overlayClipRange==value);
    indexInArray = find(params.overlayClip==value);
    switch(indexInArray)
      case 1
%         if params.overlayClipRange(indexInArray)<params.overlayUsefulRange(1)
        if params.overlayClip(indexInArray)<params.overlayRange(1)
          value = [];
        end
      case 2
%         if params.overlayClipRange(indexInArray)>params.overlayUsefulRange(2)
        if params.overlayClip(indexInArray)>params.overlayRange(2)
          value = [];
        end
    end
    if isempty(value)
      mrWarnDlg('(editOverlayGUImrParams) clip range must be within useful range');
    end
  case 'usefulrange'
    %we don't know what value we got so we need to find out first
%     indexInArray = find(params.overlayUsefulRange==value);
    indexInArray = find(params.overlayRange==value);
    switch(indexInArray)
      case 1
%         if params.overlayUsefulRange(indexInArray)>params.overlayClipRange(1)
        if params.overlayRange(indexInArray)>params.overlayClip(1)
          value = [];
        end
      case 2
%         if params.overlayUsefulRange(indexInArray)<params.overlayClipRange(2)
        if params.overlayRange(indexInArray)<params.overlayClip(2)
          value = [];
        end
    end
    if isempty(value)
      mrWarnDlg('(editOverlayGUImrParams) useful range must contain clip range');
    end  
    %check that min<max
%     if diff(params.overlayUsefulRange)<0
    if diff(params.overlayRange)<0
      mrWarnDlg('(editOverlayGUImrParams) useful range must be increasing');
      value=[];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mrCmapCallback(params,viewNum)

   thisView = viewGet(viewNum,'view');
   thisView = viewSet(thisView,'overlayCache','init');
  
  % get the current overlay
  analysisNum = viewGet(thisView,'currentAnalysis');
  overlayNum = viewGet(thisView,'currentOverlay',analysisNum);
  currentOverlay = viewGet(thisView, 'overlay', overlayNum, analysisNum);
  newOverlay = currentOverlay;
  
  %for some parameters, there is an existing viewSet
  newOverlay.interrogator = params.interrogator;
  if ~isequal(newOverlay.interrogator,currentOverlay.interrogator)
    viewSet(thisView,'interrogator',newOverlay.interrogator);
    refreshMLRDisplay(viewNum);
    return;
  end
  %overlay type
  newOverlay.overlayType = params.overlayType;
  if ~isequal(newOverlay.overlayType,currentOverlay.type)
    viewSet(thisView,'overlayType',newOverlay.overlayType);
    refreshMLRDisplay(viewNum);
    return;
  end
  % set the overlay useful range 
%   newOverlay.range = [params.overlayUsefulRange(1) params.overlayUsefulRange(2)];
  newOverlay.range = [params.overlayRange(1) params.overlayRange(2)];
  %if the range has changed, we only need to update the slider
  if any(abs(newOverlay.range-currentOverlay.range)>5e-7)
    viewSet(thisView,'overlayRange',newOverlay.range, overlayNum);
    refreshMLRDisplay(viewNum);
    return;
  end
%   % set the overlay color range 
%   newOverlay.colorRange = [params.overlayColorRange(1) params.overlayColorRange(2)];
%   if any(abs(newOverlay.colorRange-currentOverlay.colorRange)>5e-7)% && ~strcmp(newOverlay.colormapType,'normal')
%     viewSet(thisView,'overlayColorRange',newOverlay.colorRange, overlayNum);
%     refreshMLRDisplay(viewNum);
%     return;
%   end
  
  %parameters for which the overlay has to be recomputed as a whole (this should be changed by adding cases to viewSet)
  % set which color cmap to use
  if ~strcmp(params.overlayCmap, 'default')
    if sum(strcmp(params.overlayCmap, {'hsvDoubleCmap','cmapExtendedHSV','cmapHSV','overlapCmap','redGreenCmap','rygbCmap','bicolorCmap','coolCmap'}))
      newOverlay.colormap = eval(sprintf('%s(%i,%i)', params.overlayCmap, params.numGrays, params.numColors));
    else
      newOverlay.colormap = eval(sprintf('%s(%i)', params.overlayCmap, params.numColors));
    end
  end

  % see if we need to call a function
  if ~isempty(params.userDefinedCmap)
    % look for the m function
    if exist(sprintf('%s.m',params.userDefinedCmap),'file')
      colormap = eval(sprintf('%s(%i)',params.userDefinedCmap,params.numColors));
      if isequal(size(colormap),[params.numColors 3])
	newOverlay.colormap = colormap;
      else
	disp(sprintf('(editOverlay) Function %s must return a %ix%i array',params.userDefinedCmap,params.numColors,3));
      end
    end
  end
    
  % flip the cmap
  if params.flipColormap
    newOverlay.colormap = flipud(newOverlay.colormap);
  end

  % shift the cmap
  if params.shiftColormap
    newOverlay.colormap = circshift(newOverlay.colormap, params.shiftColormap);
  end

  % scale to max, or not
  newOverlay.colormapType = params.overlayColormapType;

  % set the overlay clip
%   newOverlay.clip = [params.overlayClipRange(1) params.overlayClipRange(2)];
  newOverlay.clip = [params.overlayClip(1) params.overlayClip(2)];
  

%   % set the name of the overlay
%   newOverlay.name = params.overlayName;
  if strcmp(params.alphaOverlay,'none')
    newOverlay.alphaOverlay = '';
  else
    newOverlay.alphaOverlay = params.alphaOverlay;
  end
  newOverlay.alphaOverlayExponent = params.alphaOverlayExponent;

  
  %iff the overlay has changed, 
  if ~isequalwithequalnans(newOverlay,currentOverlay)
    disppercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
    % set the new overlay
    thisView = viewSet(thisView,'newOverlay', newOverlay);
    % and refresh
    refreshMLRDisplay(thisView.viewNum);
    disppercent(inf);
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapSetManyOverlaysCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mrCmapSetManyOverlaysCallback(viewNum)
retval = [];

thisView = viewGet(viewNum,'view');

% get the current overlay
currentOverlayNum = viewGet(thisView,'curOverlay');
currentOverlayName = viewGet(thisView,'overlayname',currentOverlayNum);
currentOverlay = viewGet(thisView,'overlay',currentOverlayNum);
nOverlays = viewGet(thisView,'numOverlays');
allOverlays = viewGet(thisView,'overlays');

%get a list of overlays to change
overlayList = selectInList(thisView,'overlays','Select overlays to modify');
allOverlays = allOverlays(overlayList);
nOverlays = length(overlayList);

% a list of fields which should be copied
fieldsToCopy = {'alpha','alphaOverlayExponent','interrogator','clip','colormap','colormapType','range','colorRange'};
%remove all other fields from the currentoverlay
fieldsToRemove = setdiff(fieldnames(currentOverlay),fieldsToCopy);
currentOverlay = rmfield(currentOverlay,fieldsToRemove);

% now copy the remaining fields to all overlays
allOverlays = copyFields(currentOverlay,allOverlays,1:nOverlays);

%and set them in the view
thisView = viewSet(thisView,'newOverlay',allOverlays);

% set back the current overlay to the original one
viewSet(thisView,'curOverlay',viewGet(thisView,'overlayNum',currentOverlayName));
refreshMLRDisplay(thisView.viewNum);


%%%%%%%%%%%%%%%%%%%%
% function to close a previoulsy existing Dialog, taken form mrParamsDialog
%%%%%%%%%%%%%%%%%%%%
function closeHandler(varargin)

global gParams;
if isempty(gParams),return,end

if isfield(gParams,'fignum') 
  if isfield(gParams,'figlocstr')
  % save figure locations .mrDefaults
    for iFig = 1:length(gParams.fignum)
      mrSetFigLoc(fixBadChars(gParams.figlocstr{iFig}),get(gParams.fignum(iFig),'Position'));
    end
  end
  % close figure
  delete(gParams.fignum);
end
saveMrDefaults;

clear global gParams;
drawnow


  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% list of useful cmaps taken from MLR3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmap = grayCirc(numGrays)
  if ~iseven(numGrays);
    numGrays = numGrays + 1;
  end
  
  cmapA = gray(numGrays/2 + 1);
  cmapB = flipud(cmapA);
  
  cmap = cat(1,cmapA(1:numGrays/2,:), cmapB(1:numGrays/2,:));

function cmap = hsvDoubleCmap(numGrays,numColors,symmetric)
%
% cmap = hsvDoubleCmap(numGrays,numColors,symmetric)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   hsv colors - numGrays+1:numGrays+numColors
%
% Double wrapping the colors is useful for visualize retinotopy
% so that 180 deg from one hemifield maps onto a full hsv
% colormap (instead of just half of it).  If symmetric, flip
% second hsv cmap
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end
  if ~exist('symmetric','var')
    symmetric=1;
  end

  cmap = zeros(numGrays+numColors,3);
  if symmetric
    cmap(1:numGrays+numColors,:) = ...
        [gray(numGrays);
         hsv(floor(numColors/2));
         flipud(hsv(ceil(numColors/2)))];
  else
    cmap(1:numGrays+numColors,:) = ...
        [gray(numGrays);
         hsv(floor(numColors/2));
         hsv(ceil(numColors/2))];
  end


  
function cmap = overlapCmap(numGrays,numColors)
%
% cmap = overlapCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   3-color overlap map - numGrays+1:numGrays+numColors
%
% ras 2/04

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);zeros(numColors,3)];
  cmap(numGrays+1:numGrays+numColors/3,1) = .7;
  cmap(numGrays+numColors/3+1:numGrays+2*numColors/3,2) = .7;
  cmap(numGrays+2*numColors/3:numGrays+numColors,[3]) = .7;

  return

function cmap = redGreenCmap(numGrays,numColors)
%
% cmap = redGreenCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   redGreen ramp - numGrays+1:numGrays+numColors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays,:) = gray(numGrays);
  for i=numGrays+1:numGrays+numColors
    cmap(i,:) = ...
        [((i-numGrays)/numColors)^.5, (1-(i-numGrays)/numColors)^.5, 0];
  end


function cmap = rygbCmap(numGrays,numColors)
%
% cmap = rygbCmap(numGrays)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   red, yellow, green, blue - next 4 colors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end

  nRs=floor(1/4*numColors);
  nYs=floor(1/2*numColors)-nRs;
  nGs=floor(3/4*numColors)-nRs-nYs;
  nBs=numColors-nGs-nYs-nRs;

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays,:) = gray(numGrays);
  cmap(numGrays+1:numGrays+nRs,:) = ones(nRs,1)*[1 0 0];
  cmap(numGrays+nRs+1:numGrays+nRs+nYs,:) = ones(nYs,1)*[1 1 0];
  cmap(numGrays+nRs+nYs+1:numGrays+nRs+nYs+nGs,:) = ones(nGs,1)*[0 1 0];
  cmap(numGrays+nRs+nYs+nGs+1:numGrays+numColors,:) = ones(nBs,1)*[0 0 1];


function cmap = bicolorCmap(numGrays,numColors)
%
% cmap = bicolorCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   cool colors - values in which map < 0
%   black - values in which map=0
%   hot colors - values in which map > 0
%
% This is useful for plotting contrast maps in which both
% positive and negative effects are displayed (for related updates, see
% loadParameterMap, computeContrastMap).
%
% djh 1/98
% ras, 03/04, off of hotCmap
%   numGrays  = 64;
%   numColors = 64;

%   hi = max(max(max(map)));
%   lo = min(min(min(map)));

  
% ATTN:  for now just assume that range is [-1 1]
  hi = 1;
  lo = -1;
  
  rng = linspace(lo,hi,numColors);
  
  %   if lo > 0 % all positive
  %     colors = hot(numColors);
  %   elseif hi < 0 % all negative 
  %     colors = cool(numColors);
  %   else        % crosses zero
  colors = zeros(numColors,3);
  neg = length(find(rng < 0));
  colors(neg,:) = [0 0 0]; % black when crosses
  colors(1:neg-1,:) = flipud(cool(neg-1));
  colors(neg+1:end,:) = hot(numColors-neg);
  %   end
  
  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);colors];
  
  return

  
function cmap = twoCondCmap(numGrays)
  minx = 0;maxx = 0.5;
  x = minx:(maxx-minx)/(numGrays/2):maxx;
  y = mygauss([1 0 0.25 0],x);
  
  red = [zeros(length(y),1) y' y'];
  y = fliplr(y);
  blue = [y' y' zeros(length(y),1)];
  
  cmap = 1-[red;blue];
return

function cmap = hotColdCmap(n)


  h = hot(floor(n/2));
  c(:,1) = h(:,3);
  c(:,2) = h(:,2);
  c(:,3) = h(:,1);
  
  if iseven(h)
    cmap = [flipud(c);h];
  else
    cmap = [flipud(c);[0 0 0];h];
  end

return

function cmap = twoCondCircCmap(numGrays)
  if isodd(numGrays)
    numGrays = numGrays +1;
  end
  cmapA = twoCondCmap(numGrays/2);
  cmapB = flipud(cmapA);
  cmap = cat(1,cmapA, cmapB);
return
  
function cmap = cmapExtendedHSV(numGrays,numColors,range)
%
% cmap = cmapExtendedHSV([numGrays=128],[numColors=96],[range=query user])
% 
% Makes a special purpose colormap that is useful for representing rotating
% wedge data.  The map created here are hsv maps where all of the colors
% can be placed in a subsection of the full color map.  In this way, the
% full range of colors spans less than 2pi, like in the double color map.
% Rather than compressing by a complete factor of 2, like the double color
% map, the compression factor can be a bit smaller.
%
%   There are numGrays gray scale entries.  They occupy the first part of the cmap, 1:numGrays
%   There are numColors hsv colors.  They fill the map entries following the gray, 
%   from (numGrays+1):numGrays+numColors
%   When the range is 1, the hsv map is hsv(numColors) and we insert it
%   into the cmap.
%   When the range is 1.5, we compute tmp = hsv(numColors/1.5) and we
%   create [tmp,tmp(1:needed)] to fill up numColors entries.
%
% Examples:
%   cmap = cmapExtendedHSV(128,96,1.1);
%   cmap = cmapExtendedHSV(128,96);    -- 128 gray levels, 96 color levels,
%          query use for compression
%   cmap = cmapExtendedHSV;
%   cmap = cmapExtendedHSV(128,96,2);  -- Same as hsvDoubleCmap
%


  if ~exist('numGrays','var'),    numGrays=128; end
  if ~exist('numColors','var'),   numColors=96; end
  if ~exist('range','var'),       range = readRange;  end

  if (range < 1) | (range > 2),  error('Range must be betweem 1 and 2.'); end

  hsvColors = round(numColors/range);
  hsvColorsExtra = numColors - hsvColors;
  hsvMap = cmapHSV(hsvColors);

  cmap = zeros(numGrays+numColors,3);

  % If you want the map symmetric at the boundary, you should do this.
  % We could trap range == 2 and do it then ... which would be backwards
  % compatible?
  % cmap = [gray(numGrays); hsvMap; flipud(hsvMap(1:hsvColorsExtra,:))];
  cmap = [hsvMap; hsvMap(1:hsvColorsExtra,:)];
  shiftSize = round(hsvColorsExtra/2);
  hsvMap = circshift(cmap,shiftSize);

  cmap = [gray(numGrays); cmap];

  return;
  
function hsvMap = cmapHSV(hsvColors);
%
%   hsvMap = cmapHSV(hsvColors);
% 
%Author: AB/BW
%Purpose:
%   Create an hsv map such that the magenta/blue boundary is in the middle.
%   This makes manipulation of the colors easier for wedge maps.
%

% shiftSize to make magenta  the middle color is in the center.
  magenta = [.1 0 1];
  mp = hsv(hsvColors);

  % Obscure bit of code.  Have fun.
  [val,idx] = min(sum( abs(mp - repmat(magenta,hsvColors,1))'));
  shiftSize = -round((idx - hsvColors/2));

  hsvMap = circshift(mp,shiftSize);

  return;

function cmap = coolCmap(numGrays,numColors)
%
% cmap = coolCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   cool colors - numGrays+1:numGrays+numColors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);cool(numColors)];



%----------------------------------------
function range = readRange

  prompt={'Enter compression range for the hsv map (2 = double color map)'};
  def={'1.2'};
  dlgTitle='Color map compression factor';
  lineNo=1;
  range=inputdlg(prompt,dlgTitle,lineNo,def);
  range = str2num(range{1});

  return;


% Gaussian
%
% usage: gauss(p,X,Y);
%   p is an array of parameters:
%     p(1) = height of Gaussian
%     p(2) = center x
%     p(3) = center y
%     p(4) = SD in x dimension
%     p(5) = SD in y dimension
%     p(6) = offset
%     p(7) = rotation in radians
%   X and Y are the position on which to evaluate the gaussian
%     to evaluate at a matrix of points use,
%     e.g. [X,Y] = meshgrid(-1:.1:1,-1:.1:1);
%
%  the function can also be called as follows for 1D
%  usage: gauss(p,X);
%
%     p(1) = height of Gaussian
%     p(2) = center
%     p(3) = SD
%     p(4) = offset
% 
%   by: justin gardner
% date: 6/6/97
function G=mygauss(p,X,Y)

% 2D Gaussian 
if nargin == 3

  % rotate coordinates
  % note that the negative sign is because
  % we are rotating the coordinates and not the function
  X1 = cos(-p(7)).*X - sin(-p(7)).*Y;
  Y1 = sin(-p(7)).*X + cos(-p(7)).*Y;

  % calculate the Gaussian
  G = p(1) * exp(-((((X1-p(2)).^2)/(2*p(4)^2))+(((Y1-p(3)).^2)/(2*p(5)^2))))+p(6);
  
% 1D Gaussian
elseif nargin == 2

  % calculate the Gaussian
  G = p(1) * exp(-(((X-p(2)).^2)/(2*p(3)^2)))+p(4);
  
else 
   % usage error
   disp('USAGE: gauss(parameters, X, Y)');
end
