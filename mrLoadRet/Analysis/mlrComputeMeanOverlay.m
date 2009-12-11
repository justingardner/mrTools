% mlrComputeMeanOverlay.m
%
%        $Id:$ 
%      usage: mlrComputeMeanOverlay(sessionPath,groupNum,scanNum,analysisName,overlayName,baseAnatomy)
%         by: justin gardner
%       date: 12/11/09
%    purpose: To compute a mean overlay across subjects. The mean overlay is taken for a particular base
%             anatomy. So, if you want to take a mean across subjects, the base anatomy should be an 
%             atlas brain like one imported from Caret (see mlrImportCaret).
%
%             e.g. 
%             mlrComputeMeanOverlay({'S00320090717','S00920090717'},'Concatenation',1,'erAnal','r2','leftAtlasVeryInflated');
%
function retval = mlrComputeMeanOverlay(sessionPath, groupNum, scanNum, analysisName, overlayName, baseAnatomies)

% check arguments
if ~any(nargin == [6])
  help mlrComputeMeanOverlay
  return
end

% parse input arguments
sessions = parseArguments(sessionPath,groupNum,scanNum,analysisName,overlayName,baseAnatomies);
if isempty(sessions),return,end

% go through each session and load the overlays
overlays = loadOverlays(sessions);
if isempty(overlays),return,end

% compute the average overlays
meanOverlays = computeMeanOverlays(sessions,overlays);

% save the overlays back
saveMeanOverlays(sessions,meanOverlays);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   saveMeanOverlays   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveMeanOverlays(sessions,meanOverlays)

for iSession = 1:sessions.nSessions
  % open a view on to the session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % set the group
  v = viewSet(v,'curGroup',sessions.groupNum{iSession});
  % set the scan number
  v = viewSet(v,'curScan',sessions.scanNum{iSession});
  % get scan dims
  scanDims = viewGet(v,'scanDims');
  % now project the overlay back into the scan coordinates
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b = viewGet(v,'baseCoordMap');
    % get the transformation from base2scan coordinates
    base2scan = viewGet(v,'base2scan');
    % get the coordinates of each vertex
    baseCoords = squeeze(b.coords);
    baseCoords(:,4) = 1;
    % convert into scan coords
    scanCoords = round(base2scan*baseCoords');
    scanLinearCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
    % create an empty overlay
    overlay{iBase} = nan(scanDims);
    % get the overlay points that lay within the image
    overlayIm = meanOverlays{iBase}.overlayIm;
    % get points outside of scan and remove
    goodPoints = find(~isnan(scanLinearCoords));
    scanLinearCoords = scanLinearCoords(goodPoints);
    overlayIm = overlayIm(goodPoints);
    % and put the overlay into the scan coordinate iamge
    overlay{iBase}(scanLinearCoords) = overlayIm;
  end
  % and install the base into a custom analysis
  mrDispOverlay(overlay,sessions.scanNum{iSession},viewGet(v,'groupNum',sessions.groupNum{iSession}),[],'overlayName=meanOverlay','saveName=meanAnalysis');
  % close and delete
  deleteView(v);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeMeanOverlay    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meanOverlays = computeMeanOverlays(sessions,overlays)

% NOTE: This averaging won't work for amp/phase maps like the corAnal which have
% to be converted back to a complex number and averaged.

for iBase = 1:sessions.nBases
  for iSession = 1:sessions.nSessions
    % sum up overlays across all sessions
    if iSession == 1
      meanOverlays{iBase}.overlayIm = overlays.o(iSession,iBase).overlayIm;
    else
      meanOverlays{iBase}.overlayIm = meanOverlays{iBase}.overlayIm + overlays.o(iSession,iBase).overlayIm;
    end
  end
  % divide by n
  meanOverlays{iBase}.overlayIm = meanOverlays{iBase}.overlayIm/sessions.nSessions;
  % convert to RGB
  %clim = [min(meanOverlays{iBase}.overlayIm) max(meanOverlays{iBase}.overlayIm)];
  clim = overlays.o(iSession,iBase).range;
  meanOverlays{iBase}.overlayRGB = rescale2rgb(meanOverlays{iBase}.overlayIm,overlays.o(iSession,iBase).cmap,clim);
  % and display
  dispOverlay(overlays.b(iSession,iBase),meanOverlays{iBase}.overlayRGB,sprintf('Mean on %s',sessions.baseAnatomies{iBase}));
end


%%%%%%%%%%%%%%%%%%%%%%
%    loadOverlays    %
%%%%%%%%%%%%%%%%%%%%%%
function overlays = loadOverlays(sessions)

overlays = [];

for iSession = 1:sessions.nSessions
  mrQuit([]);
  % first make sure the directory exists
  if ~myisdir(sessions.sessionPath{iSession})
    disp(sprintf('(mlrComputeMeanOverlay:loadOverlays) Could not find directory %s',sessions.sessionPath{iSession}));
    return
  end
  % next make sure it has a session
  if ~isfile(fullfile(sessions.sessionPath{iSession},'mrSession.mat'))
    disp(sprintf('(mlrComputeMeanOverlay:loadOverlays) Cound not find mrSession.mat in %s',sessions.sessionPath{iSession}));
    return
  end
  % create a view for each session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % set the group
  v = viewSet(v,'curGroup',sessions.groupNum{iSession});
  % set the scan number
  v = viewSet(v,'curScan',sessions.scanNum{iSession});
  % load the analysis
  v = loadAnalysis(v,sessions.analysisName{iSession});
  if viewGet(v,'numAnalyses') < 1,return,end
  % load the overlay, just to make sure it is there
  v = viewSet(v,'curOverlay',sessions.overlayName{iSession});
  o{iSession} = viewGet(v,'overlay');
  if isempty(o{iSession}),return,end
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b(iSession,iBase) = viewGet(v,'baseCoordMap');
    if isempty(b(iSession,iBase)), disp(sprintf('(mlrComputeMeanOverlay:loadOverlays) Could not load base'));return,end
    % compute the overlay
    [overlayImage base roi overlay(iSession,iBase)] = refreshMLRDisplay(viewGet(v,'viewNum'));
    if isempty(overlay(iSession,iBase)),return,end
    % display in a figure
    dispOverlay(b(iSession,iBase),overlay(iSession,iBase).RGB,makeSessionName(v,sessions,iSession));
  end
  % delete the view
  deleteView(v);
end

% check that all the bases match
for iBase = 1:sessions.nBases
  for iSession = 2:sessions.nSessions
    if size(b(1,iBase).outerVtcs,1) ~= size(b(iSession,iBase).outerVtcs,1) 
      disp(sprintf('(mlrComputeMeanOverlay:loadOverlays) Base %s does not match number of vertices for %s',sessions.baseAnatomies{iBase},makeSessionName([],sessions,iSession)));
      return
    end
  end
end
% set return fields
overlays.b = b;
overlays.o = overlay;

%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeSessionName    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function s = makeSessionName(v,sessions,i)

if isempty(v)
  s = sprintf('%s: %s:%i %s:%s',sessions.sessionPath{i},sessions.groupNum{i},sessions.scanNum{i},sessions.analysisName{i},sessions.overlayName{i});
else
  s = sprintf('%s: %s:%i %s:%s',sessions.sessionPath{i},viewGet(v,'groupName',sessions.groupNum{i}),sessions.scanNum{i},sessions.analysisName{i},viewGet(v,'overlayName',sessions.overlayName{i}));
end
  
%%%%%%%%%%%%%%%%%%%%%
%    dispOverlay    %
%%%%%%%%%%%%%%%%%%%%%
function dispOverlay(b,overlay,titleStr)

if exist('smartfig') == 2
  smartfig(fixBadChars(sprintf('mlrComputeMeanOverlay:%s',titleStr)),'reuse');clf
else
  figure;
end
patch('vertices', b.innerVtcs, 'faces', b.tris,'FaceVertexCData', squeeze(overlay),'facecolor','interp','edgecolor','none');
axis equal;
axis off;
camva(6);
rotate3d on;
title(titleStr);
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%
%    parseArguments    %
%%%%%%%%%%%%%%%%%%%%%%%%
function sessions = parseArguments(sessionPath,groupNum,scanNum,analysisName,overlayName,baseAnatomies)

sessions = [];

% check length of arguments, if we only have one, then copy for all scans
nSessions = length(sessionPath);
if nSessions <= 1
  disp(sprintf('(mlrComputeMeanOverlay) Only %i sessions passed in',nSessions));
  return
end

% validate argument length of rest of arguments
groupNum = validateArgumentLength('groupNum',groupNum,nSessions);
if isempty(groupNum),return,end
scanNum = validateArgumentLength('scanNum',scanNum,nSessions);
if isempty(scanNum),return,end
analysisName = validateArgumentLength('analysisName',analysisName,nSessions);
if isempty(analysisName),return,end
overlayName = validateArgumentLength('overlayName',overlayName,nSessions);
if isempty(overlayName),return,end

% pack into output and return
sessions.nSessions = nSessions;
sessions.sessionPath = sessionPath;
sessions.groupNum = convertToCellArray(groupNum);
sessions.scanNum = convertToCellArray(scanNum);
sessions.analysisName = convertToCellArray(analysisName);
sessions.overlayName = convertToCellArray(overlayName);
sessions.baseAnatomies = cellArray(baseAnatomies);
sessions.nBases = length(sessions.baseAnatomies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateArgumentLength    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arg = validateArgumentLength(argname,arg,nSessions)


% copy arg to all, if it is a single string that is passed in
if isstr(arg)
  thisArg = arg;
  arg = {};
  for i = 1:nSessions
    arg{i} = thisArg;
  end
elseif length(arg) == 1
  % arg is a number
  if isnumeric(arg)
    arg(2:nSessions) = arg(1);
  end
end

% check length
if length(arg) ~= nSessions
  disp(sprintf('(mlrComputeMeanOverlay) %s length (%i) does not match number of sessions (%i)',argname,length(arg),nSessions));
  arg = [];
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    convertToCellArray    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = convertToCellArray(a)

if ~iscell(a)
  for i = 1:length(a)
    b{i} = a(i);
  end
  a = b;
end

%%%%%%%%%%%%%%%%%
%    myisdir    %
%%%%%%%%%%%%%%%%%
function tf = myisdir(dirname)

tf = (length(dir(dirname)) ~= 0);
  