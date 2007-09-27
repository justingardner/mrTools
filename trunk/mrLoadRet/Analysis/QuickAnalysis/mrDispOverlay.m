% mrDispOverlay.m
%
%      usage: [view analysis] = mrDispOverlay(overlay,scanNum,groupNum/analysisStruct,<view>)
%      usage: mrDispOverlay(overlay,scanNum,analysisStructure,<view>)
%         by: justin gardner
%       date: 04/04/07
%    purpose: displays an overlay in MrLoadRet
%
%             For example, to display a map for scan 1, group 2:
%             Where map has dimensions size(map) = [x y s]
% 
%             mrDispOverlay(map,1,2)
%
%             If you want to install maps for several scans, make
%             map into a cell array of maps of the correct
%             dimensions and do
%
%             mrDispOverlay(map,[1 2],2);
%           
%             If you want to install the overlay into an existing
%             analysis, rather than a group (say erAnal), pass
%             that analysis instead of the group
% 
%             mrDispOverlay(map,1,erAnal);
%
%             If you want to install the map into an existing
%             view, you can pass that (or empty if you don't want
%             to display at all)
%
%             mrDispOverlay(map,1,erAnal,[]);
%
%             You can also set optional parameters
%
%             Some parameters you can set
%             'overlayName=name'
%             'cmap',colormap
%             'colormapType=setRangeToMax'
%             'range',[0 1]
%             'clip',[0 1]
%             'interrogator=interrogatorName'
%             'analName=analysisName'
%             'd',d 
%
%             to save instead of display set
%             'saveName=filename'
%
%             e.g.
%             mrDispOverlay(map,1,erAnal,[],'overlayName=myMap');
%
%
function [v, analysis] = mrDispOverlay(overlay,scanNum,groupNum,v,varargin)

% check arguments
if nargin < 3
  help mrDispOverlay
  return
end

% get the arguments
eval(evalargs(varargin));

mrLoadRetViewing = 0;
viewShouldBeDeleted = 0;
% start up a mrLoadRet if we are not passed in a view
% or if we are returning an analysis struct (in this
if ieNotDefined('v')
  if ieNotDefined('saveName') && (nargout < 2)
    v = mrLoadRet;
    mrLoadRetViewing = 1;
  else
    % if we are saving, then don't bring up mrLoadRet
    v = newView('Volume');
    viewShouldBeDeleted = 1;
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

% make the overlay into a cell array, if it is not already.
overlay = cellArray(overlay);

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
o.(overlayName).reconcileFunction = 'defaultReconcileParams';
o.(overlayName).data = cell(1,viewGet(v,'nScans'));
o.(overlayName).params.scanNum = [];
for i = 1:length(scanNum)
  if length(overlay) >= i
    o.(overlayName).data{scanNum(i)} = overlay{i};
    o.(overlayName).params.scanNum(end+1) = i;
  % if we are passed in one overlay but multiple scans
  % then we install the same overlay for all the scans
  elseif length(overlay) == 1
    disp(sprintf('(mrDispOverlay) Only 1 overlay passed in. Installing that overlay for scan %i',scanNum(i)));
    o.(overlayName).data{scanNum(i)} = overlay{1};
    o.(overlayName).params.scanNum(end+1) = scanNum(i);
  % otherwise something is wrong
  else
    disp(sprintf('(mrDispOverlay) Only %i overlays found for %i scans',length(overlay),length(scanNum)));
  return
  end
end
o.(overlayName).date = dateString;
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
  a.(analName).guiFunction = '';
  a.(analName).params.scanNum = scanNum;
  a.(analName).reconcileFunction = 'defaultReconcileParams';
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
  if ~ieNotDefined('d')
    a.(analName).d{scanNum} = d;
  end
  v = viewSet(v,'newAnalysis',a.(analName));
% or just install this overlay into the current analysis
else
  % install overlay
  if ~isempty(overlay)
    v = viewSet(v,'newOverlay',o.(overlayName));
  end
  % install d structure
  if ~ieNotDefined('d')
    v = viewSet(v,'newd',d,scanNum);
  end
end

% either save it or refresh the display
if ~ieNotDefined('saveName')
  v = viewSet(v,'analysisName',saveName);
  saveAnalysis(v,saveName);
end

% update mr load ret if it is running
if mrLoadRetViewing
  v = viewSet(v,'curGroup',groupNum);
  mlrGuiSet(v,'group',groupNum);
  mlrGuiSet(v.viewNum,'scan',scanNum(1));
  refreshMLRDisplay(v.viewNum);
end

% prepare the output argument
if nargout >= 2
  analysis = viewGet(v,'analysis',viewGet(v,'curAnalysis'));
end

% delete view
if viewShouldBeDeleted
  deleteView(v);
end
