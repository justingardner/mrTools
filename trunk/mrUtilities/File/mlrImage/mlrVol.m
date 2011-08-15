% mlrVol.m
%
%        $Id:$ 
%      usage: mlrVol(filename)
%         by: justin gardner
%       date: 08/15/11
%    purpose: Volume viewer to display volume multi-dimensional data (usually 3D volumes or 4D epi sequences)
%
function retval = mlrVol(filename,varargin)

% check arguments
if nargin < 1
  help mlrVol
  return
end

if isstr(filename)
  % init the system
  sysNum = initSystem(varargin);

  % load the volume
  loadVolume(filename,sysNum);

  % display the volume
  refreshDisplay(sysNum);
else
  % command argument
  switch (filename)
   case 1
    mouseDownHandler(varargin{1});
   case 2
    mouseMoveHandler(varargin{1});
   case 3
    textHandler(varargin{1});
   case 4
    incdecHandler(varargin{1},varargin{2},varargin{3});
   case 5
    buttonHandler(varargin{1},varargin{2});
   case 6
    closeHandler(varargin{1});
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%%   closeHandler   %%
%%%%%%%%%%%%%%%%%%%%%%
function closeHandler(sysNum)

global gSystem;

% close the figures
for i = 1:length(gSystem{sysNum}.fig)
  figure(gSystem{sysNum}.fig(i));closereq;
%  smartfig('close',gSystem{sysNum}.fig(i));
%  smartfig('close',gSystem{sysNum}.fig(i));
end

% remove the variable
gSystem = {gSystem{1:sysNum-1} gSystem{sysNum+1:end}};

%%%%%%%%%%%%%%%%%%%%%%%
%%   buttonHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function buttonHandler(textNum,sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% if we are already running the animation, then turn it off
if gSystem{sysNum}.animating == textNum
  gSystem{sysNum}.animating = 0;
  return
end

% get number of dimensions
nDim = gSystem{sysNum}.vols(1).h.nDim;

% start the animation
gSystem{sysNum}.animating = textNum;
coord = vol.coord;

% loop that runs the animation
while(gSystem{sysNum}.animating)
  % increment coord
  coord(textNum) = coord(textNum)+1;
  if coord(textNum) > vol.h.dim(textNum)
    coord(textNum) = 1;
  end
  % set the volume coord
  setVolCoord(sysNum,gSystem{sysNum}.n,coord);
  % set the text boxes
  for iCoord = 1:nDim
    set(gSystem{sysNum}.hCoordTextbox(iCoord),'String',coord(iCoord));
  end
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(coord(1),coord(2),coord(3),coord(4),coord(5)));
  % and refresh
  refreshDisplay(sysNum);
  % make it ability to be seen
  pause(0.02);
end
%%%%%%%%%%%%%%%%%%%%%%%
%%   incdecHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function incdecHandler(textNum,incdec,sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

gSystem{sysNum}.animating = false;

% inc or dec the text box
val = mystr2num(get(gSystem{sysNum}.hCoordTextbox(textNum),'String'));
val = val+incdec;

% if we have made a valid change then set it
if (val>=1) && (val<=vol.h.dim(textNum))
  set(gSystem{sysNum}.hCoordTextbox(textNum),'String',val);
  % now refresh
  textHandler(sysNum);
end

%%%%%%%%%%%%%%%%%%%%%
%%   textHandler   %%
%%%%%%%%%%%%%%%%%%%%%
function textHandler(sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% turn off animation
gSystem{sysNum}.animating = 0;

% get number of dimensions
nDim = gSystem{sysNum}.vols(1).h.nDim;

% get coordinate
for iCoord = 1:nDim
  coord(iCoord) = mystr2num(get(gSystem{sysNum}.hCoordTextbox(iCoord),'String'));
end
coord = round(coord);

if any(coord<1) || any(coord(:)>vol.h.dim(1:nDim))
  for iCoord = 1:nDim
    set(gSystem{sysNum}.hCoordTextbox(1),'String',vol.coord(iCoord));
  end
  return
end

% nDims hard coded to 5 here
coord(end+1:5) = 1;

% set the volume coord
setVolCoord(sysNum,gSystem{sysNum}.n,coord);

% nDims hard coded to 5 here
vol = gSystem{sysNum}.vols(1);
set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(vol.coord(1),vol.coord(2),vol.coord(3),vol.coord(4),vol.coord(5)));

% and redisplay
refreshDisplay(sysNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mouseMoveHadnler    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseMoveHandler(sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

coord = getMouseCoord(sysNum);
if isempty(coord)
  set(gcf,'Pointer','arrow');
  set(gSystem{sysNum}.hCoordTextbox(1),'String',vol.coord(1));
  set(gSystem{sysNum}.hCoordTextbox(2),'String',vol.coord(2));
  set(gSystem{sysNum}.hCoordTextbox(3),'String',vol.coord(3));
  % nDims hard coded to 5 here
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(vol.coord(1),vol.coord(2),vol.coord(3),vol.coord(4),vol.coord(5)));
else
  set(gcf,'Pointer','fullcrosshair');
  set(gSystem{sysNum}.hCoordTextbox(1),'String',coord(1));
  set(gSystem{sysNum}.hCoordTextbox(2),'String',coord(2));
  set(gSystem{sysNum}.hCoordTextbox(3),'String',coord(3));
  % nDims hard coded to 5 here
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(coord(1),coord(2),coord(3),coord(4),coord(5)));
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mouseDownHandler    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% stop any ongoing animation
gSystem{sysNum}.animating = 0;

% get the mouse pos
coord = getMouseCoord(sysNum);
if isempty(coord),return,end

% and set the volume coordinate
setVolCoord(sysNum,gSystem{sysNum}.n,coord);

% refresh the display
refreshDisplay(sysNum);

%%%%%%%%%%%%%%%%%%%%%%%
%    getMouseCoord    %
%%%%%%%%%%%%%%%%%%%%%%%
function coord = getMouseCoord(sysNum)

coord = [];viewNum = [];

global gSystem;
vol = gSystem{sysNum}.vols(1);

% figure out which axis we are on
pointerLoc = get(gcf,'CurrentPoint');
pos = get(gcf,'Position');
pos = pointerLoc./pos(3:4);
subplotNum = ceil(pos(1)*3);
if (subplotNum>0) && (subplotNum<=3)
  a = subplot(1,3,subplotNum);
else
  return
end
viewNum = find(gSystem{sysNum}.a == a);

% get pointer loc
pointerLoc = get(a,'CurrentPoint');
pointerLoc = round(pointerLoc(1,2:-1:1));
% check image boundaries
if (pointerLoc(1) < 1) | (pointerLoc(1) > vol.h.dim(vol.xDim(viewNum))),return,end
if (pointerLoc(2) < 1) | (pointerLoc(2) > vol.h.dim(vol.yDim(viewNum))),return,end

% now get the point on the image
coord(vol.xDim(viewNum)) = pointerLoc(1);
coord(vol.yDim(viewNum)) = pointerLoc(2);
coord(vol.viewDim(viewNum)) = vol.coord(vol.viewDim(viewNum));

% nDims hard coded to 5 here
coord(end+1:5) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%    refreshDisplay    %
%%%%%%%%%%%%%%%%%%%%%%%%
function refreshDisplay(sysNum)

global gSystem;

f = gcf;a = gca;
for iVol = 1:gSystem{sysNum}.n
  % display the coordinate
  dispVolume(gSystem{sysNum}.vols(iVol))
  % set that we have displayed this coordinate
  gSystem{sysNum}.vols(iVol).curCoord = gSystem{sysNum}.vols(iVol).coord;
end

drawnow;
figure(f);%axes(a);
%%%%%%%%%%%%%%%%%%%%
%    dispVolume    %
%%%%%%%%%%%%%%%%%%%%
function dispVolume(vol)

global gSystem;

for iView = 1:3
  % see if we need to redisplay
  if ~isequal(vol.curCoord(4:end),vol.coord(4:end)) || ~isequal(vol.curCoord(vol.viewDim(iView)),vol.coord(vol.viewDim(iView)))
    % get the slice
    % nDims hard coded to 5 here
    dispSlice = squeeze(vol.data(vol.viewIndexes{iView,1},vol.viewIndexes{iView,2},vol.viewIndexes{iView,3},vol.coord(4),vol.coord(5)));

    % clip
    %  dispSlice = clipImage(dispSlice);
  
    % get the correct axis to draw into
    figure(vol.fig(iView));
    a = subplot(1,3,vol.subplotNum(iView));

    % and display the image
    cla(a);imagesc(dispSlice,'Parent',a);
    axis equal
    axis off
  
    % set title
    titleStr = sprintf('%s (%i %i %i)',vol.h.filename,vol.coord(1),vol.coord(2),vol.coord(3));
    h = title(titleStr,'Interpreter','none');
  end
end

%%%%%%%%%%%%%%%%%%%
%    clipImage    %
%%%%%%%%%%%%%%%%%%%
function img = clipImage(img)

% Choose a sensible clipping value
histThresh = length(img(:))/1000;
[cnt, val] = hist(img(:),100);
goodVals = find(cnt>histThresh);
if isempty(goodVals)
  clipMin = 0;clipMax = 0;
else
  clipMin = val(min(goodVals));
  clipMax = val(max(goodVals));
end

% and convert the image
img(img<clipMin) = clipMin;
img(img>clipMax) = clipMax;
if (clipMax-clipMin) > 0
  img = 255*(img-clipMin)./(img-clipMin);
end

%%%%%%%%%%%%%%%%%%%%
%    loadVolume    %
%%%%%%%%%%%%%%%%%%%%
function loadVolume(filename,sysNum)

global gSystem;

[d h] = mlrLoadImage(filename);
if isempty(d),return,end

% updated number of volumes
gSystem{sysNum}.n = gSystem{sysNum}.n+1;

% update the fields in vols
gSystem{sysNum}.vols(gSystem{sysNum}.n).data = d;
gSystem{sysNum}.vols(gSystem{sysNum}.n).h = h;

% set which figure numbers this volume will display into
gSystem{sysNum}.vols(gSystem{sysNum}.n).fig(1:3) = gSystem{sysNum}.fig(1:3);
%gSystem{sysNum}.vols(gSystem{sysNum}.n).subplotNum(1:3) = gSystem{sysNum}.n;
gSystem{sysNum}.vols(gSystem{sysNum}.n).subplotNum(1:3) = 1:3;

% FIX: these should be set with reference to permutation matrix
gSystem{sysNum}.vols(gSystem{sysNum}.n).viewDim(1:3) = 1:3;

% for convenience set which dimension is x and y for each of the views
% given the viewDim info set above
for i = 1:3
  otherDims = setdiff(1:3,gSystem{sysNum}.vols(gSystem{sysNum}.n).viewDim(i));
  gSystem{sysNum}.vols(gSystem{sysNum}.n).xDim(i) = otherDims(1);
  gSystem{sysNum}.vols(gSystem{sysNum}.n).yDim(i) = otherDims(2);
end

% update the coordinate (start displaying in middle of volume)
coord = round(h.dim/2);
% nDims hard coded to 5 here
coord(end+1:5) = 1;
setVolCoord(sysNum,gSystem{sysNum}.n,coord);

% now this sets the indexes from the volume for which the
% image will be displayed
for iView = 1:3
  for jAxis = 1:3
    if jAxis ~= gSystem{sysNum}.vols(gSystem{sysNum}.n).viewDim(iView)
      gSystem{sysNum}.vols(gSystem{sysNum}.n).viewIndexes{iView,jAxis} = 1:gSystem{sysNum}.vols(gSystem{sysNum}.n).h.dim(jAxis);
    else
      gSystem{sysNum}.vols(gSystem{sysNum}.n).viewIndexes{iView,jAxis} = gSystem{sysNum}.vols(gSystem{sysNum}.n).coord(jAxis);
    end
  end
end

% set the current coordinates for the first time
gSystem{sysNum}.vols(gSystem{sysNum}.n).curCoord = [nan nan nan];

% print out the header of the image
dispHeaderInfo(gSystem{sysNum}.vols(gSystem{sysNum}.n));

% set the text boxes
if gSystem{sysNum}.n == 1
  names = {'x','y','z','slice','receiver'};
  for i = 1:h.nDim
    % make the inc/dec textboxes
    gSystem{sysNum}.hCoordTextbox(i) = makeTextboxIncDec(sysNum,1,coord(i),i,1,i);
    % get name for button
    if i < length(names)
      name = names{i};
    else 
      name = sprintf('dim %i',i);
    end
    % make the button
    gSystem{sysNum}.hButton(i) = makeButton(sysNum,1,name,i,2,i);
  end
  % make a text box for the value of the voxel
  gSystem{sysNum}.hValTextbox(1) = makeTextbox(sysNum,1,'',1,6);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispHeaderInfo   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function dispHeaderInfo(vol)

dispHeader(vol.h.filename);
disp(sprintf('type: %s (%s)',vol.h.type,vol.h.ext));
disp(sprintf('dim: [%s]',mynum2str(vol.h.dim(:)','sigfigs=0')));
disp(sprintf('pixdim: [%s]',mynum2str(vol.h.pixdim(:)')));
disp(sprintf('qform_code: %i',vol.h.qform_code));
disp(sprintf('qform:'));
disp(sprintf('%s',mynum2str(vol.h.qform44,'compact=0','sigfigs=-1')));
disp(sprintf('sform_code: %i',vol.h.sform_code));
disp(sprintf('sform:'));
disp(sprintf('%s',mynum2str(vol.h.sform44,'compact=0','sigfigs=-1')));
dispHeader;

%%%%%%%%%%%%%%%%%%%%%
%    setVolCoord    %
%%%%%%%%%%%%%%%%%%%%%
function setVolCoord(sysNum,iVol,coord)

global gSystem;

% set the current x,y,z coordinate (start in middle)
gSystem{sysNum}.vols(iVol).coord = coord;

% now this sets the indexes from the volume for which the
% image will be displayed
for iView = 1:3
  for jAxis = 1:3
    if jAxis ~= gSystem{sysNum}.vols(iVol).viewDim(iView)
      gSystem{sysNum}.vols(iVol).viewIndexes{iView,jAxis} = 1:gSystem{sysNum}.vols(iVol).h.dim(jAxis);
    else
      gSystem{sysNum}.vols(iVol).viewIndexes{iView,jAxis} = gSystem{sysNum}.vols(iVol).coord(jAxis);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%    initSystem    %
%%%%%%%%%%%%%%%%%%%%
function sysNum = initSystem(args)

global gSystem;
if isempty(gSystem)
  sysNum = 1;
else
  sysNum = length(gSystem)+1;
end

% number of loaded volumes
gSystem{sysNum}.n = 0;

% empty vols
vols = [];

% defaults for button sizes
gSystem{sysNum}.buttonWidth = 100;
gSystem{sysNum}.buttonWidthMargin = 20;
gSystem{sysNum}.buttonHeightMargin = 2;
gSystem{sysNum}.buttonHeight = 25;
gSystem{sysNum}.buttonLeftMargin = 10;
gSystem{sysNum}.buttonBottomMargin = 5;

% open the fig
gSystem{sysNum}.fig(1) = smartfig('mlrVol');
gSystem{sysNum}.fig(2) = gSystem{sysNum}.fig(1);
gSystem{sysNum}.fig(3) = gSystem{sysNum}.fig(1);
clf;colormap(gray);

% set the mouse functions
set(gSystem{sysNum}.fig(1),'WindowButtonDownFcn',sprintf('mlrVol(1,%i)',sysNum));
set(gSystem{sysNum}.fig(1),'WindowButtonMotionFcn',sprintf('mlrVol(2,%i)',sysNum));
set(gSystem{sysNum}.fig(1),'CloseRequestFcn',sprintf('mlrVol(6,%i)',sysNum));

% set up the axis
gSystem{sysNum}.a(1) = subplot(1,3,1);cla;axis off;
gSystem{sysNum}.a(2) = subplot(1,3,2);cla;axis off;
gSystem{sysNum}.a(3) = subplot(1,3,3);cla;axis off;

% no animation is running
gSystem{sysNum}.animating = false;

% update display
drawnow

%%%%%%%%%%%%%%%%%%%%%
%%   makeTextbox   %%
%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(sysNum,viewNum,displayString,rownum,colnum);

h = uicontrol('Style','edit','String',displayString,'Position',getUIControlPos(sysNum,viewNum,rownum,colnum,1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(3,%i)',sysNum));

%%%%%%%%%%%%%%%%%%%%%%%
%%   makePushbuton   %%
%%%%%%%%%%%%%%%%%%%%%%%
function h = makeButton(sysNum,viewNum,displayString,textboxNum,rownum,colnum)

h = uicontrol('Style','pushbutton','String',displayString,'Position',getUIControlPos(sysNum,viewNum,rownum,colnum,1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(5,%i,%i)',textboxNum,sysNum));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeTextboxIncDec   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = makeTextboxIncDec(sysNum,viewNum,displayString,textboxNum,rownum,colnum)

% make textbox
h = uicontrol('Style','edit','String',displayString,'Position',getUIControlPos(sysNum,viewNum,rownum,colnum+.125,0.75),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(3,%i)',sysNum));

% put up incdec buttons
uicontrol('Style','pushbutton','String','<','Position',getUIControlPos(sysNum,viewNum,rownum,colnum,.1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(4,%i,-1,%i)',textboxNum,sysNum));
uicontrol('Style','pushbutton','String','>','Position',getUIControlPos(sysNum,viewNum,rownum,colnum+.9,.1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(4,%i,1,%i)',textboxNum,sysNum));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(sysNum,viewNum,rownum,colnum,uisize)

global gSystem;

% get figure position
figpos = get(gSystem{sysNum}.fig(viewNum),'Position');

% set this buttons width
thisButtonWidth = gSystem{sysNum}.buttonWidth*uisize;

% set the position for the button
pos(1) = (gSystem{sysNum}.buttonWidth+gSystem{sysNum}.buttonWidthMargin)*(floor(colnum)-1)+gSystem{sysNum}.buttonLeftMargin+(colnum-floor(colnum))*gSystem{sysNum}.buttonWidth;
pos(2) = gSystem{sysNum}.buttonBottomMargin + (gSystem{sysNum}.buttonHeight+gSystem{sysNum}.buttonHeightMargin)*(rownum-1)+gSystem{sysNum}.buttonHeight;
pos(3) = thisButtonWidth;
pos(4) = gSystem{sysNum}.buttonHeight;
