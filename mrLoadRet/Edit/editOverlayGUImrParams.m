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
  % Also, get all colormap functions located in the mrLoadRet colormap folder
  functionsDirectory = [fileparts(which('mrLoadRet')) '/colormapFunctions/'];
  colormapFunctionFiles =  dir([functionsDirectory '*.m']);
  for iFile=1:length(colormapFunctionFiles)
    colormapFunctions{iFile} = stripext(colormapFunctionFiles(iFile).name);
  end
  colormaps = union(colormaps,colormapFunctions,'stable');
  
  % set up params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'overlayCmap', colormaps,'type=popupmenu','List of possible colormaps'};
  paramsInfo{end+1} = {'userDefinedCmap','','Allows you to call a user defined function to set the overlap colormap. The user-defined function must output a params.numColors x 3 RGB colormap and can be either an anonymous function accepting the number of color as single input argument or a function on the Matlab path accepting the number of colors as its last argument. In the latter case, additional input arguments must be specified after the function name, separated with commas. '};
  paramsInfo{end+1} = {'numColors', numColors, 'first argument to the colormap function'};
  paramsInfo{end+1} = {'numGrays', 0, 'second argument to the colormap function'};
  paramsInfo{end+1} = {'flipColormap', 0, 'type=checkbox', 'check this box to reverse the direction of the colormap'};
  paramsInfo{end+1} = {'shiftColormap', 0, 'incdec=[-16 16]', 'shift the colormap -- this can be useful for retinotopy scans with circular colormaps'}; 
  paramsInfo{end+1} = {'overlayColormapType',overlayColormapTypeMenu , 'type=popupmenu',...
      '''normal'' scales the colormap to the value specified by colormap range; ''setRangeToMax'' scales the colormap to the min amnd max of the displayed overlay slice and ignores the color range (like R2 maps)'};
  switch(mrGetPref('overlayRangeBehaviour'))
    case 'Classic'
      paramsInfo{end+1} = {'overlayRange', overlayRange, 'callback',{@checkCmapParams,'range'},'passCallbackOutput=1','passValue=1','passParams=1',...
          'The lower and upper bounds on the clip slider, and of the colormap when overlayColormapType=''normal''.'};
      paramsInfo{end+1} = {'overlayClip', overlayClipRange, ...
          'The lower and upper bounds beyond which the overlay is masked, and of the colormap when overlayColormapType=''setRangeToMax''.. If clip(1)>clip(2), then values inside the clip range are masked.'};
    case 'New'
      paramsInfo{end+1} = {'overlayColorRange', overlayColorRange, 'callback',{@checkCmapParams,'colorrange'},'passCallbackOutput=1','passValue=1','passParams=1',...
          'The lower and upper bounds of the colormap when overlayColormapType=''normal'''};
      paramsInfo{end+1} = {'overlayClipRange', overlayClipRange, 'callback',{@checkCmapParams,'cliprange'},'passCallbackOutput=1','passValue=1','passParams=1',...
          'The lower and upper bounds beyond which the overlay is masked. These should be inside the sliders range. If clip(1)>clip(2), then values inside the clip range are masked.'};
      paramsInfo{end+1} = {'overlaySlidersRange', overlayRange, 'callback',{@checkCmapParams,'slidersrange'},'passCallbackOutput=1','passValue=1','passParams=1',...
          'The lower and upper bounds on the clip sliders. These should be lower/higher than the clip values'};
      paramsInfo{end+1} = {'setOverlaySlidersRange', 0, 'type=pushbutton','callback',@mrCmapsetOverlaySlidersRange,'callbackArg',viewNum,'buttonString=Set sliders range to overlay min/max','passParams=1','passCallbackOutput=0',...
          'Sets the sliders range to the min/max values of this overlay accross scans'};
  end
  paramsInfo{end+1} = {'overlayType', overlayType, 'The type of overlay (ph, amp, co ...'};
  paramsInfo{end+1} = {'interrogator', interrogator, 'Sets the overlay default interrogator function'};
  paramsInfo{end+1} = {'alphaOverlay', alphaOverlayMenu, 'You can specify the name of another overlay in the analysis to use as an alpha map. For instance, you might want to display one overlay with the alpha set to the r2 or the coherence value.'};
  paramsInfo{end+1} = {'alphaOverlayExponent', alphaOverlayExponent,'incdec=[-0.1 0.1]','If you are using an alphaOverlay, this sets an exponent on the alphaOverlay to pass through. For example, if you just want the value from the overlay to be the alpha value then set this to 1. If you want to have it so that lower values get accentuated (usually this is the case), set the exponent to less than 1, but greater than 0. The alpha values are passed through the function alpha = alpha.^alphaOverlayExponent'};
  paramsInfo{end+1} = {'setManyOverlays', 0, 'type=pushbutton','Set many overlays to have them same settings as the current overlay. This cannot be cancelled.','callback',@mrCmapSetManyOverlaysCallback,'callbackArg',viewNum,'buttonString=Set many overlays'};
%   paramsInfo{end+1} = {'overlayName', overlayName, 'The name for the overlay (e.g., co, am, or ph)'};

  % display dialog
  mrParamsDialog(paramsInfo,'Edit overlay parameters','modal=0','callback',@mrCmapCallback,'callbackArg',viewNum,'cancelCallback',{@mrCmapParamsCancel,overlay});

  return;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  mrCmapsetOverlaySlidersRange   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapsetOverlaySlidersRange(viewNum,params)

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
  mlrDispPercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
  % set the new overlay
  thisView = viewSet(thisView,'newOverlay', oldOverlay);
  % and refresh
  refreshMLRDisplay(thisView.viewNum);
  mlrDispPercent(inf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  checkCmapParams   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = checkCmapParams(params,value,whichParam)

switch(whichParam)
  case 'colorrange'
    if diff(params.overlayColorRange)<=0
      mrWarnDlg('(editOverlayGUImrParams) colorange values cannot be equal or decreasing.');
      value=[];
    end  
  case 'cliprange'
    %we don't know what value we got so we need to find out first
    indexInArray = find(params.overlayClipRange==value);
    switch(indexInArray)
      case 1
        if params.overlayClipRange(indexInArray)<params.overlaySlidersRange(1)
          value = [];
        end
      case 2
        if params.overlayClipRange(indexInArray)>params.overlaySlidersRange(2)
          value = [];
        end
    end
    if isempty(value)
      mrWarnDlg('(editOverlayGUImrParams) Clip range must be within sliders range');
    end
  case 'slidersrange'
    %we don't know what value we got so we need to find out first
    indexInArray = find(params.overlaySlidersRange==value);
    % for new overlay range behaviour, check that clip sliders range includes clip range
    switch(indexInArray)
      case 1
        if params.overlaySlidersRange(indexInArray)>min(params.overlayClipRange)
          value = [];
        end
      case 2
        if params.overlaySlidersRange(indexInArray)<max(params.overlayClipRange)
          value = [];
        end
    end
    if isempty(value)
      mrWarnDlg('(editOverlayGUImrParams) Sliders range must contain clip range');
    end  
    %check that min<max
    if diff(params.overlaySlidersRange)<0
      mrWarnDlg('(editOverlayGUImrParams) Sliders range must be increasing');
      value=[];
    end
  case 'range'
    if diff(params.overlayRange)<=0
      mrWarnDlg('(editOverlayGUImrParams) Range values cannot be equal or decreasing.');
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
  newOverlay.type = params.overlayType;
  if ~isequal(newOverlay.type,currentOverlay.type)
    viewSet(thisView,'overlayType',newOverlay.type);
    refreshMLRDisplay(viewNum);
    return;
  end
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
  % set 
  if any(abs(newOverlay.colorRange-currentOverlay.colorRange)>5e-7)% && ~strcmp(newOverlay.colormapType,'normal')
    if strcmp(mrGetPref('overlayRangeBehaviour'),'Classic')
      % in the 'Classic' behaviour, we set both the sliders and the colour range (to the same values)
      viewSet(thisView,'overlayRange',newOverlay.range, overlayNum);
    end
    % in the 'New behaviour, we only set the color range
    viewSet(thisView,'overlayColorRange',newOverlay.colorRange, overlayNum);
    refreshMLRDisplay(viewNum);
    return;
  end
  % in the 'New' behaviour, the sliders range could have changed independently of the color range
  if any(abs(newOverlay.range-currentOverlay.range)>5e-7)
      % but we only need to update the slider
      viewSet(thisView,'overlayRange',newOverlay.range, overlayNum);
%       refreshMLRDisplay(viewNum); %so no need to refresh the display
    return;
  end
  
  %parameters for which the overlay has to be recomputed as a whole (this should be changed by adding cases to viewSet)
  % set which color cmap to use
  if ~strcmp(params.overlayCmap, 'default')
    if sum(strcmp(params.overlayCmap, {'hsvDoubleCmap','cmapExtendedHSV','overlapCmap','redGreenCmap','rygbCmap','bicolorCmap','coolCmap'}))
      newOverlay.colormap = eval(sprintf('%s(%i,%i)', params.overlayCmap, params.numGrays, params.numColors));
    else
      newOverlay.colormap = eval(sprintf('%s(%i)', params.overlayCmap, params.numColors));
    end
  end

  % see if we need to call a function
  if ~isempty(params.userDefinedCmap)
    %parse the function name and its arguments
    if params.userDefinedCmap(1)=='@' % if this is an anonyous function
      try
        cMapfunction = str2func(params.userDefinedCmap);
        colormap = cMapfunction(params.numColors);
      catch
        fprintf('(editOverlay) Anonymous function %s returned an error\n',params.userDefinedCmap);
      end
    else
      cMapFunction=textscan(params.userDefinedCmap,'%s','delimiter',',');
      cMapFunction=cMapFunction{1};
      % look for the m function
      if exist(sprintf('%s.m',cMapFunction{1}),'file')
        cMapFunction{1} = str2func(cMapFunction{1}); %convert function string fo function handle
        for iArg =2:length(cMapFunction) %if there are additional arguments, convert numerical ones
          if ~isempty(str2num(cMapFunction{iArg}))
            cMapFunction{iArg} = str2num(cMapFunction{iArg});
          end
        end
        cMapFunction{end+1}=params.numColors; %add number of colors as last argument
        colormap = callbackEval(cMapFunction);
      end
    end
    if isequal(size(colormap),[params.numColors 3])
      newOverlay.colormap = colormap;
    else
      fsprintf('(editOverlay) Function %s must return a %ix%i array\n',params.userDefinedCmap,params.numColors,3);
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

  
  %if the overlay has changed, 
  if ~isequalwithequalnans(newOverlay,currentOverlay)
    mlrDispPercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
    % set the new overlay
    thisView = viewSet(thisView,'newOverlay', newOverlay);
    % and refresh
    refreshMLRDisplay(thisView.viewNum);
    mlrDispPercent(inf);
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

