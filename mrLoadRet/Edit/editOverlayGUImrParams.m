% editOverlayGUImrParams.m
%
%      usage: editOverlayGUImrParams(viewNum)
%         by: eli merriam, modified by julien besle (in 2010 and 2016)
%       date: 09/04/07
%    purpose: 
%
function editOverlayGUImrParams(viewNum)

% check arguments
  if ~any(nargin == [1])
    help editOverlayGUImrParams
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
  if isempty(analysisNum),mrWarnDlg('(editOverlayGUImrParams) No current analysis');return,end
  overlayNum = viewGet(thisView,'currentOverlay', analysisNum);
  if isempty(overlayNum),mrWarnDlg('(editOverlayGUImrParams) No current overlay');return,end
  if length(overlayNum)>1
    mrWarnDlg('(editOverlayGUImrParams) Not implemented for several overlays');
    %close aEdit Overlay Dialog if it's on (that could be the case if overlay is changed while the edit overlay params dialog is already on)
    global gParams
    if ~isempty(gParams) && strcmp(gParams.figlocstr{1},'mrParamsDialog_Change_overlay_colormap')
      closeHandler;
    end
    return;
  end
  overlay = viewGet(thisView, 'overlay', overlayNum, analysisNum);
  cmapCallbackArg.viewNum = viewNum;
  cmapCallbackArg.originalColormap = overlay.colormap;
  if strcmp(mrGetPref('overlayRangeBehaviour'),'New')
      overlayColorRange = viewGet(thisView,'overlayColorRange', overlayNum, analysisNum);
  end
  overlayRange = viewGet(thisView,'overlayRange', overlayNum, analysisNum);
  overlayClipRange = viewGet(thisView,'overlayClip', overlayNum, analysisNum);
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
  colormaps = {'default','hot','hsv','pink','cool','bone','copper','flag','gray','jet'};
  altColormaps = viewGet(thisView,'colormaps');
  if ~isempty(altColormaps)
    colormaps = {colormaps{:} altColormaps{:}};
  end
  
  % set up params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'overlayCmap', colormaps,'type=popupmenu','List of possible colormaps'};
  paramsInfo{end+1} = {'userDefinedCmap','','Allows you to call a user defined function to set the overlap colormap. You can specify additional input arguments after the function name, separating with commas. The user-defined function must be on the path, accept an integer representing the number of colors as its last input argument, and output a params.numColors*3 RGB colormap'};
  paramsInfo{end+1} = {'numColors', numColors, 'first argument to the colormap function'};
  paramsInfo{end+1} = {'numGrays', 0, 'second argument to the colormap function'};
  paramsInfo{end+1} = {'flipColormap', 0, 'type=checkbox', 'check this box to reverse the direction of the colormap'};
  paramsInfo{end+1} = {'shiftColormap', 0, 'incdec=[-16 16]', 'shift the colormap -- this can be useful for retinotopy scans with circular colormaps'}; 
  paramsInfo{end+1} = {'overlayColormapType',overlayColormapTypeMenu , 'type=popupmenu',...
      '''normal'' scales the colormap to the value specified by colormap range; ''setRangeToMax'' scales the colormap to the min amnd max of the displayed overlay slice and ignores the color range (like R2 maps)'};
  switch(mrGetPref('overlayRangeBehaviour'))
    case 'Classic'
      paramsInfo{end+1} = {'overlayRange', overlayRange, ...
          'The lower and upper bounds on the clip slider, and of the colormap when overlayColormapType=''normal''.'};
      paramsInfo{end+1} = {'overlayClip', overlayClipRange, ...
          'The lower and upper bounds beyond which the overlay is masked, and of the colormap when overlayColormapType=''setRangeToMax''.. If clip(1)>clip(2), then values inside the clip range are masked.'};
    case 'New'
      paramsInfo{end+1} = {'overlayColorRange', overlayColorRange, ...
          'The lower and upper bounds of the colormap when overlayColormapType=''normal'''};
      paramsInfo{end+1} = {'overlayClipRange', overlayClipRange, ...
          'The lower and upper bounds beyond which the overlay is masked. These should be inside the sliders range. If clip(1)>clip(2), then values inside the clip range are masked.'};
      paramsInfo{end+1} = {'overlaySlidersRange', overlayRange, ...
          'The lower and upper bounds on the clip sliders. These should be lower/higher than the clip values'};
      paramsInfo{end+1} = {'setOverlaySlidersRange', 0, 'type=pushbutton','callback',@mrCmapsetOverlaySlidersRange,'callbackArg',cmapCallbackArg,'buttonString=Set sliders range to overlay min/max','passParams=1','passCallbackOutput=0',...
          'Sets the sliders range to the min/max values of this overlay accross scans'};
  end
  paramsInfo{end+1} = {'overlayType', overlayType, 'The type of overlay (ph, amp, co ...'};
  paramsInfo{end+1} = {'interrogator', interrogator, 'Sets the overlay default interrogator function'};
  paramsInfo{end+1} = {'alphaOverlay', alphaOverlayMenu, 'You can specify the name of another overlay in the analysis to use as an alpha map. For instance, you might want to display one overlay with the alpha set to the r2 or the coherence value.'};
  paramsInfo{end+1} = {'alphaOverlayExponent', alphaOverlayExponent,'incdec=[-0.1 0.1]','If you are using an alphaOverlay, this sets an exponent on the alphaOverlay to pass through. For example, if you just want the value from the overlay to be the alpha value then set this to 1. If you want to have it so that lower values get accentuated (usually this is the case), set the exponent to less than 1, but greater than 0. The alpha values are passed through the function alpha = alpha.^alphaOverlayExponent'};
  paramsInfo{end+1} = {'setManyOverlays', 0, 'type=pushbutton','Set many overlays to have them same settings as the current overlay. This cannot be cancelled.','callback',@mrCmapSetManyOverlaysCallback,'callbackArg',viewNum,'buttonString=Set many overlays'};
%   paramsInfo{end+1} = {'overlayName', overlayName, 'The name for the overlay (e.g., co, am, or ph)'};

  % display dialog
  mrParamsDialog(paramsInfo,'Edit overlay parameters','modal=0','callback',@mrCmapCallback,'callbackArg',cmapCallbackArg,'cancelCallback',{@mrCmapParamsCancel,overlay});

  return;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapsetOverlaySlidersRange   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapsetOverlaySlidersRange(cmapCallbackArg,params)

viewNum = cmapCallbackArg.viewNum;
thisView = viewGet(viewNum,'view');

minOverlayData = floor(double(viewGet(thisView,'minOverlayData'))*1e6)/1e6;
maxOverlayData = ceil(double(viewGet(thisView,'maxOverlayData'))*1e6)/1e6;

if isempty(maxOverlayData) || isempty(minOverlayData)
  mrWarnDlg('(editOverlayGUImrParams) overlay seems to be empty');
elseif all(abs(params.overlaySlidersRange-[minOverlayData maxOverlayData])<5e-7)
  mrWarnDlg('(editOverlayGUImrParams) Sliders range is already set to min/max');
else
  params.overlaySlidersRange = [minOverlayData maxOverlayData];
  %make sure clip values are within this range
  params.overlayClipRange(1) = max(params.overlayClipRange(1),minOverlayData);
  params.overlayClipRange(2) = min(params.overlayClipRange(2),maxOverlayData);
  mrParamsSet(params);
  mrCmapCallback(params,cmapCallbackArg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapParamsCancel   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapParamsCancel(oldOverlay,cmapCallbackArg)

viewNum = cmapCallbackArg.viewNum;
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
function valid = checkCmapParams(params)

valid = false;

switch(mrGetPref('overlayRangeBehaviour'))
  case 'Classic'
    if numel(params.overlayRange) ~= 2 || diff(params.overlayRange) <= 0
      mrWarnDlg('(editOverlayGUImrParams) Range values cannot be equal or decreasing.');
      return
    end

  case 'New'
    if numel(params.overlayColorRange) ~= 2 || diff(params.overlayColorRange) <= 0
      mrWarnDlg('(editOverlayGUImrParams) Color range values cannot be equal or decreasing.');
      return
    end
    if numel(params.overlaySlidersRange) ~= 2 || diff(params.overlaySlidersRange) < 0
      mrWarnDlg('(editOverlayGUImrParams) Sliders range cannot be decreasing.');
      return
    end
    if numel(params.overlayClipRange) ~= 2 || ...
        any(params.overlayClipRange < params.overlaySlidersRange(1)) || ...
        any(params.overlayClipRange > params.overlaySlidersRange(2))
      mrWarnDlg('(editOverlayGUImrParams) Clip range must be within sliders range.');
      return
    end
end

valid = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mrCmapCallback(params,cmapCallbackArg)

  viewNum = cmapCallbackArg.viewNum;

  % Validate complete range pairs here rather than validating each edit box
  % separately. This lets users edit both endpoints before applying them.
  if ~checkCmapParams(params)
    return
  end

  thisView = viewGet(viewNum,'view');
  thisView = viewSet(thisView,'overlayCache','init');
  
  % get the current overlay
  analysisNum = viewGet(thisView,'currentAnalysis');
  overlayNum = viewGet(thisView,'currentOverlay',analysisNum);
  currentOverlay = viewGet(thisView, 'overlay', overlayNum, analysisNum);
  newOverlay = currentOverlay;
  newOverlay.colormap = cmapCallbackArg.originalColormap;
  
  % Apply all dialog values to one overlay structure so one callback commits
  % the complete edit and refreshes only once.
  newOverlay.interrogator = params.interrogator;
  newOverlay.type = params.overlayType;

  % set the overlay (sliders) range and the overlay color range 
  switch(mrGetPref('overlayRangeBehaviour'))
    case 'Classic'
      % in the classic behaviour, we make sure the colour range follows the sliders range
      newOverlay.range = [params.overlayRange(1) params.overlayRange(2)];
      newOverlay.colorRange = [params.overlayRange(1) params.overlayRange(2)];
    case 'New'
      newOverlay.range = [params.overlaySlidersRange(1) params.overlaySlidersRange(2)];
      newOverlay.colorRange = [params.overlayColorRange(1) params.overlayColorRange(2)];
  end
  
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
    %parse the function name and its arguments
    cMapFunction=textscan(params.userDefinedCmap,'%s','delimiter',',');
    cMapFunction=cMapFunction{1};
    
    % look for the m function
    if exist(sprintf('%s.m',cMapFunction{1}),'file')
      cMapFunction{1} = str2func(cMapFunction{1}); %convert function string fo function handl
      for iArg =2:length(cMapFunction) %if there are additional arguments, convert numerical ones
        if ~isempty(str2num(cMapFunction{iArg}))
          cMapFunction{iArg} = str2num(cMapFunction{iArg});
        end
      end
      cMapFunction{end+1}=params.numColors; %add number of colors as last argument
      colormap = callbackEval(cMapFunction);
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
  switch(mrGetPref('overlayRangeBehaviour'))
    case 'Classic'
      newOverlay.clip = [params.overlayClip(1) params.overlayClip(2)];
    case 'New'
      newOverlay.clip = [params.overlayClipRange(1) params.overlayClipRange(2)];
  end
  

%   % set the name of the overlay
%   newOverlay.name = params.overlayName;
  if strcmp(params.alphaOverlay,'none')
    newOverlay.alphaOverlay = '';
  else
    newOverlay.alphaOverlay = params.alphaOverlay;
  end
  newOverlay.alphaOverlayExponent = params.alphaOverlayExponent;

  
  % If the overlay has changed, reinstall it. Refresh unconditionally below
  % because the overlay cache was cleared at the start of this callback.
  overlayChanged = ~isequalwithequalnans(newOverlay,currentOverlay);
  if overlayChanged
    disppercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
    % set the new overlay
    thisView = viewSet(thisView,'newOverlay', newOverlay);
  end

  % Redraw even when the edit is a no-op (for example, a zero colormap
  % shift), so that the cleared overlay cache is actually recomputed.
  refreshMLRDisplay(thisView.viewNum);

  if overlayChanged
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
