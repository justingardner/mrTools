% mrEpiMovie.m
%
%        $Id$
%      usage: mrEpiMovie(v)
%         by: justin gardner
%       date: 08/16/07
%    purpose: display epis as a movie for inspection
%
function retval = mrEpiMovie(v)

% check arguments
if ~any(nargin == [1])
  help mrEpiMovie
  return
end

global gMrEpiMovie;

% keep v
gMrEpiMovie.v = v;

% setup a cache
gMrEpiMovie.c = mrCache('init',1000);

% set animating flag
gMrEpiMovie.animating = 0;
gMrEpiMovie.stopAnimating = 0;

% get max frames and slices
maxFrames = -inf;
maxSlices = -inf;
for scanNum = 1:viewGet(v,'nScans')
  maxFrames = max(maxFrames,viewGet(v,'nFrames',scanNum));
  maxSlices = max(maxSlices,viewGet(v,'nSlices',scanNum));
end

% set up params dialog
paramsInfo = {};
paramsInfo{end+1} = {'scanNum',1,sprintf('minmax=[1 %i]',viewGet(v,'nScans')),sprintf('incdec=[-1 1]'),'round=1','Choose scan to view'};
paramsInfo{end+1} = {'scanNumMovie',0,'type=pushbutton','buttonString=Animate over scans','callback',@mrEpiMovieAnimate,'passParams=1','callbackArg','scanNum','Press to animate over scans'};
paramsInfo{end+1} = {'sliceNum',1,sprintf('minmax=[1 %i]',maxSlices),sprintf('incdec=[-1 1]'),'round=1','Choose slice number to view'};
paramsInfo{end+1} = {'sliceNumMovie',0,'type=pushbutton','buttonString=Animate over slices','callback',@mrEpiMovieAnimate,'passParams=1','callbackArg','sliceNum','Press to animate over slices'};
paramsInfo{end+1} = {'frameNum',1,sprintf('minmax=[1 %i]',maxFrames),sprintf('incdec=[-1 1]'),'round=1','Choose frame to view'};
paramsInfo{end+1} = {'frameNumMovie',0,'type=pushbutton','buttonString=Animate over frames','callback',@mrEpiMovieAnimate,'passParams=1','callbackArg','frameNum','Press to animate over frames'};

% display dialog
[gMrEpiMovie.f params] = mrParamsDialog(paramsInfo,'epiMovie',[],@mrEpiMovieCallback,[],@mrEpiMovieClose);

% and draw first frame
mrEpiMovieDispImage(params.scanNum,params.sliceNum,params.frameNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrEpiFrameNumMovie   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mrEpiMovieAnimate(type,params)

retval = [];

global gMrEpiMovie;
v = gMrEpiMovie.v;

% check to see whether we are already running an animation
if gMrEpiMovie.animating
  gMrEpiMovie.stopAnimating = 1;
  return
else
  gMrEpiMovie.animating = 1;
  gMrEpiMovie.stopAnimating = 0;
end

switch type
  case {'scanNum'}
   n = viewGet(v,'nScans');
  case {'sliceNum'}
   n = viewGet(v,'nSlices',params.scanNum);
  case {'frameNum'}
   n = viewGet(v,'nFrames',params.scanNum);
end

% set what to view
scanNum = params.scanNum;
sliceNum = params.sliceNum;
frameNum = params.frameNum;

for i = 1:n
  % check to see if we have to return
  if gMrEpiMovie.stopAnimating
    gMrEpiMovie.animating = 0;
    return
  end
  % set the variable
  eval(sprintf('%s = i;',type));
  % load image
  mrEpiMovieDispImage(scanNum,sliceNum,frameNum);
end

% go back to original frame
mrEpiMovieCallback(params);

gMrEpiMovie.animating = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrEpiMovieCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrEpiMovieCallback(params)

global gMrEpiMovie;
if gMrEpiMovie.animating
  gMrEpiMovie.stopAnimating = 1;
  return
end

% draw the image
mrEpiMovieDispImage(params.scanNum,params.sliceNum,params.frameNum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrEpiMovieGetImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epiImage = mrEpiMovieDispImage(scanNum,sliceNum,frameNum)

global gMrEpiMovie;
v = gMrEpiMovie.v;

% check input arguments
frameNum = min(frameNum,viewGet(v,'nFrames',scanNum));
sliceNum = min(sliceNum,viewGet(v,'nSlices',scanNum));

% check cache for image
cacheStr = sprintf('%i_%i_%i',scanNum,sliceNum,frameNum);
[epiImage gMrEpiMovie.c] = mrCache('find',gMrEpiMovie.c,cacheStr);

if isempty(epiImage)
  % read from disk
  epiImage = loadTSeries(v,scanNum,sliceNum,frameNum);
  % and save in cache
  gMrEpiMovie.c = mrCache('add',gMrEpiMovie.c,cacheStr,epiImage);
end

selectGraphWin;

% display image
imagesc(epiImage);
colormap('gray');
axis equal; axis off
title(sprintf('Scan: %i slice %i\nframe %i',scanNum,sliceNum,frameNum));
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrEpiMovieClose   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function mrEpiMovieClose

% clear the global
global gMrEpiMovie;
gMrEpiMovie.stopAnimating = 1;
gMrEpiMovie.v = [];
gMrEpiMovie.c = [];

% close the graph window
selectGraphWin;
closeGraphWin;
