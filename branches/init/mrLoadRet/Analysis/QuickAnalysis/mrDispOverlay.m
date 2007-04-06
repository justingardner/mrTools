% mrDispOverlay.m
%
%      usage: mrDispOverlay(overlay,scanNum,groupNum,<view>,<parameters>
%      usage: mrDispOverlay(overlay,scanNum,analysisStructure,<view>,<parameters>)
%         by: justin gardner
%       date: 04/04/07
%    purpose: displays an overlay in MrLoadRet
%             if called without a view will run
%             mrLoadRet to display a view. If you pass
%             scan and group, it will load the overlay into
%             a bogus analysis and then display it.
% 
%             If you pass it an analysis sturcture then
%             it will add the overlay to the analysis strucutre
%             and display that. 
%          
%             Some parameters you can set
%             'overlayName=name'
%             'cmap',colormap
%             'colormapType=setRangeToMax'
%             'range',[0 1]
%             'clip',[0 1]
%             'interrogator=interrogatorName'
%             'analName=analysisName'
%
%             to save instead of display set
%             'saveName=filename'
%
function analOut = mrDispOverlay(overlay,scanNum,groupNum,v,varargin)

% check arguments
if nargin < 3
  help mrDispOverlay
  return
end

% get the arguments
eval(evalargs(varargin));

% start up a mrLoadRet if we are not passed in an open one
if ieNotDefined('v')
  if ieNotDefined('saveName') && (nargout < 1)
    v = mrLoadRet;
  else
    % if we are saving, then don't bring up mrLoadRet
    v = newView('Volume');
  end
end

% if groupNum is actually a structure it means we were passed
% in an analysis strututre
if isstruct(groupNum)
  anal = groupNum;
  groupNum = viewGet(v,'groupNum',anal.groupName);
  v = viewSet(v,'curGroup',groupNum);
  v = viewSet(v,'newAnalysis',anal);
end

% now set to the current view and group
v = viewSet(v,'curGroup',groupNum);
mlrGuiSet(v.viewNum,'scanText',scanNum(1));
groupName = viewGet(v,'groupName',groupNum);

% get some default parameters
if ieNotDefined('overlayName')
  overlayName = 'mrDispOverlay';
end
if ieNotDefined('cmap')
  % colormap is made with a little bit less on the dark end
  cmap = hot(312);
  cmap = cmap(end-255:end,:);
end
if ieNotDefined('colormapType')
  colormapType = 'setRangeToMax';
end
if ieNotDefined('interrogator')
  if ieNotDefined('anal')
    interrogator = 'mrDefaultInterrogator';
  else
    if isfield(anal,'overlays')
      interrogatorFound = 0;
      for i = 1:length(anal.overlays)
	if strcmp(anal.overlays(i).name,overlayName)
	  interrogator = anal.overlays(i).interrogator;
	end
      end
      if ~interrogatorFound
	interrogator = 'mrDefaultInterrogator';
      end
    end
  end
end
if ieNotDefined('range')
  range = [0 1];
end
if ieNotDefined('clip')
  clip = [0 1];
end
if ieNotDefined('analName')
  analName = 'mrDispOverlayAnal';
end

% create the parameters for the overlay
dateString = datestr(now);
o.(overlayName).name = overlayName;
o.(overlayName).function = '';
o.(overlayName).groupName = groupName;
o.(overlayName).reconcileFunction = 'quickReconcile';
o.(overlayName).data = cell(1,viewGet(v,'nScans'));
for i = 1:length(scanNum)
  % if overlay is a cell array then it means to set a different
  % overlay for each scan
  if iscell(overlay)
    o.(overlayName).data{scanNum(i)} = overlay{i};
  % otherwise set them all with the same scan
  else
    o.(overlayName).data{scanNum(i)} = overlay;
  end
end
o.(overlayName).date = dateString;
o.(overlayName).params = [];
o.(overlayName).range = range;
o.(overlayName).clip = clip;
o.(overlayName).colormap = cmap;
o.(overlayName).alpha = 1;
o.(overlayName).colormapType = colormapType;
o.(overlayName).interrogator = interrogator;

% see if we need to install a bogus analysis
if isempty(viewGet(v,'curAnalysis'))
  a.(analName).name = analName;  % This can be reset by editAnalysisGUI
  a.(analName).type = analName;
  a.(analName).groupName = groupName;
  a.(analName).function = '';
  a.(analName).reconcileFunction = 'quickReconcile';
  a.(analName).guiFunction = '';
  a.(analName).params = [];
  % if user didn't pass in an overlay, they are probably
  % using this to install something into the d strucutre
  % instead.
  if ~isempty(overlay)
    a.(analName).overlays = o.(overlayName);
  else
    a.(analName).overlays = [];
  end    
  a.(analName).curOverlay = 1;
  a.(analName).date = dateString;
  v = viewSet(v,'newAnalysis',a.(analName));
% or just install this overlay into the current analysis
else
  if ~isempty(overlay)
    v = viewSet(v,'newOverlay',o.(overlayName));
  end
end

% if we have an output argument then just return the new analysis
% structure, if we have a savename then save it, otherwise display it
if nargout == 1
  analOut = viewGet(v,'curAnalysis');
elseif ~ieNotDefined('saveName')
 v = viewSet(v,'analysisName',saveName);
 saveAnalysis(v,saveName);
else
  refreshMLRDisplay(v.viewNum);
end