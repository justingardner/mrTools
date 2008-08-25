% mlrDisplayEPI.m
%
%        $Id$
%      usage: mlrDisplayEPI(<v>,<groupNum>)
%         by: justin gardner
%       date: 08/16/07
%    purpose: display epis as a movie for inspection
%
function retval = mlrDisplayEPI(v,groupNum)

% check arguments
if ~any(nargin == [0 1 2])
  help mlrDisplayEPI
  return
end

global gMLRDisplayEPI;

if ieNotDefined('v'),
  v = newView;
  gMLRDisplayEPI.deleteViewAtEnd = 1;
else
  gMLRDisplayEPI.deleteViewAtEnd = 0;
end

if ~ieNotDefined('groupNum')
  v = viewSet(v,'curGroup',groupNum);
end

% keep v
gMLRDisplayEPI.v = v;

% setup a cache
gMLRDisplayEPI.c = mrCache('init',1000);

% set animating flag
gMLRDisplayEPI.animating = 0;
gMLRDisplayEPI.stopAnimating = 0;

% get interpMethod
gMLRDisplayEPI.interpMethod = mrGetPref('interpMethod');

% get max frames and slices
maxFrames = -inf;
maxSlices = -inf;
for scanNum = 1:viewGet(v,'nScans')
  maxFrames = max(maxFrames,viewGet(v,'nFrames',scanNum));
  maxSlices = max(maxSlices,viewGet(v,'nSlices',scanNum));
end

% get descriptions and scan2scan
needsWarping = 0;
for scanNum = 1:viewGet(v,'nScans')
  descriptions{scanNum} = viewGet(v,'description',scanNum);
  scan2scan{scanNum} = viewGet(v,'scan2scan',1,[],scanNum);
  if ~isequal(scan2scan{scanNum},eye(4))
    needsWarping = 1;
  end
end

% set up params dialog
paramsInfo = {};
paramsInfo{end+1} = {'scanNum',viewGet(v,'currentScan'),sprintf('minmax=[1 %i]',viewGet(v,'nScans')),sprintf('incdec=[-1 1]'),'round=1','Choose scan to view'};
paramsInfo{end+1} = {'scanDescription',descriptions,'editable=0','group=scanNum','type=string','Description for scan'};
paramsInfo{end+1} = {'scanNumMovie',0,'type=pushbutton','buttonString=Animate over scans','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','scanNum','Press to animate over scans'};
paramsInfo{end+1} = {'sliceNum',viewGet(v,'curSlice'),sprintf('minmax=[1 %i]',maxSlices),sprintf('incdec=[-1 1]'),'round=1','Choose slice number to view'};
paramsInfo{end+1} = {'sliceNumMovie',0,'type=pushbutton','buttonString=Animate over slices','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','sliceNum','Press to animate over slices'};
paramsInfo{end+1} = {'frameNum',1,sprintf('minmax=[1 %i]',maxFrames),sprintf('incdec=[-1 1]'),'round=1','Choose frame to view'};
paramsInfo{end+1} = {'frameNumMovie',0,'type=pushbutton','buttonString=Animate over frames','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','frameNum','Press to animate over frames'};
if needsWarping
  paramsInfo{end+1} = {'warp',0,'type=checkbox','Apply warping implied by scan2scan transform to each image. This allows you to preview what will happen when you apply warping in averageTSeries or concatTSeries'};
  paramsInfo{end+1} = {'warpBaseScan',1,sprintf('minmax=[1 %i]',viewGet(v,'nScans')),'incdec=[-1 1]','Choose which scan you want to warp the images to','contingent=warp'};
end

% display dialog
[gMLRDisplayEPI.f params] = mrParamsDialog(paramsInfo,sprintf('EPI images for %s',viewGet(v,'groupName')),[],@mlrDisplayEPICallback,[],@mlrDisplayEPIClose);

% and draw first frame
mlrDisplayEPIDispImage(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mlrDisplayEPIAnimate  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mlrDisplayEPIAnimate(type,params)

retval = [];

global gMLRDisplayEPI;
v = gMLRDisplayEPI.v;

% check to see whether we are already running an animation
if gMLRDisplayEPI.animating
  gMLRDisplayEPI.stopAnimating = 1;
  return
else
  gMLRDisplayEPI.animating = 1;
  gMLRDisplayEPI.stopAnimating = 0;
end

switch type
  case {'scanNum'}
   n = viewGet(v,'nScans');
  case {'sliceNum'}
   n = viewGet(v,'nSlices',params.scanNum);
  case {'frameNum'}
   n = viewGet(v,'nFrames',params.scanNum);
end

i = 1;
while (gMLRDisplayEPI.stopAnimating ~=1)
  % set the variable
  params.(type) = i;
  % set the dialog appropriately
  mrParamsSet(params);
  % load image
  mlrDisplayEPIDispImage(params);
  pause(.02);
  i = (mod(i,n))+1;
end

gMLRDisplayEPI.animating = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPICallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrDisplayEPICallback(params)

global gMLRDisplayEPI;
if gMLRDisplayEPI.animating
  gMLRDisplayEPI.stopAnimating = 1;
  return
end

% draw the image
mlrDisplayEPIDispImage(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPIGetImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epiImage = mlrDisplayEPIDispImage(params)

global gMLRDisplayEPI;
v = gMLRDisplayEPI.v;

if ~isfield(params,'warp'),params.warp = 0;end
  
% check input arguments
nFrames = viewGet(v,'nFrames',params.scanNum);
nSlices = viewGet(v,'nSlices',params.scanNum);
params.frameNum = min(params.frameNum,nFrames);
params.sliceNum = min(params.sliceNum,nSlices);

% check cache for image
if params.warp
  cacheStr = sprintf('%i_%i_%i_%i',params.scanNum,params.sliceNum,params.frameNum,params.warpBaseScan);
else
  cacheStr = sprintf('%i_%i_%i',params.scanNum,params.sliceNum,params.frameNum);
end
[epiImage gMLRDisplayEPI.c] = mrCache('find',gMLRDisplayEPI.c,cacheStr);

if isempty(epiImage)
  % apply warping if necessary
  if params.warp
    % get scan2scan
    scan2scan = viewGet(v,'scan2scan',params.warpBaseScan,[],params.scanNum);
    
    if ~isequal(scan2scan,eye(4))
      % load the volume
      epiVolume = loadTSeries(v,params.scanNum,[],params.frameNum);

      % swapXY is needed because warpAffine3 is in yx not xy
      swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];

      % check to see if we also want to try to do a motion comp correction
      %if params.warpWithMotionComp
      %baseVolume = loadTSeries(v,params.warpBaseScan,[],params.frameNum);
      %scan2scan
      %  scan2scan = estMotionIter3(baseVolume,epiVolume,3,scan2scan)
      %end

      % compute transform
      M = swapXY * scan2scan * swapXY;

      % display transformation
      for rownum = 1:4
	disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
      end
      disppercent(-inf,sprintf('Warping scan %i to match scan %i with transformation using %s',params.scanNum,params.warpBaseScan,gMLRDisplayEPI.interpMethod));
      epiVolume = warpAffine3(epiVolume,M,NaN,0,gMLRDisplayEPI.interpMethod);

      disppercent(inf);
      epiImage = epiVolume(:,:,params.sliceNum);
    else
      % scan2scan was identity, so no warping is necessary
      epiImage = loadTSeries(v,params.scanNum,params.sliceNum,params.frameNum);
    end
  else
    % no warping needed, just load from disk
    epiImage = loadTSeries(v,params.scanNum,params.sliceNum,params.frameNum);
  end

  
  % Choose a sensible clipping value
  histThresh = length(epiImage(:))/1000;
  [cnt, val] = hist(epiImage(:),100);
  goodVals = find(cnt>histThresh);
  if isempty(goodVals)
    clipMin = 0;clipMax = 0;
  else
    clipMin = val(min(goodVals));
    clipMax = val(max(goodVals));
  end

  % and convert the image
  epiImage(epiImage<clipMin) = clipMin;
  epiImage(epiImage>clipMax) = clipMax;
  if (clipMax-clipMin) > 0
    epiImage = 255*(epiImage-clipMin)./(clipMax-clipMin);
  end
  % and save in cache
  gMLRDisplayEPI.c = mrCache('add',gMLRDisplayEPI.c,cacheStr,epiImage);
end

selectGraphWin;

% always make it so that the h dimension is longer than v dimension
if size(epiImage,1) > size(epiImage,2)
  epiImage = epiImage';
end

% display image
image(epiImage);
colormap(gray(256));
axis equal; axis tight; axis off
title(sprintf('Scan: %i slice %i/%i frame %i/%i',params.scanNum,params.sliceNum,nSlices,params.frameNum,nFrames));
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPIClose   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrDisplayEPIClose

global gMLRDisplayEPI;

% delete the view if necessary
if gMLRDisplayEPI.deleteViewAtEnd
  deleteView(gMLRDisplayEPI.v);
end

% clear the global
gMLRDisplayEPI.stopAnimating = 1;
gMLRDisplayEPI.v = [];
gMLRDisplayEPI.c = [];

% close the graph window
selectGraphWin;
closeGraphWin;

