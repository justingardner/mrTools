% mlrVol.m
%
%        $Id:$ 
%      usage: mlrVol(filename,<filename2>,'imageOrientation=0')
%         by: justin gardner
%       date: 08/15/11
%    purpose: Volume viewer to display volume multi-dimensional data (usually 3D volumes or 4D epi sequences)
%             If qform information exists (from a nifti header) it will orient the axis in LPI and label axis
%             correctly (left/right, posterior/anterior, inferior/superior). If you want to view the
%             image in its native orientation (i.e. order it is written into the matrix, then set the
%             input argument:
%
%             'imageOrienation=1'
%
%             With two filename arguments, will show the alignment between the two images, using
%             the sforms if they are both set / or the qforms otherwise. If neither one is set
%             will display with an identity alignment. A dialog box allows you to fiddle manually
%             with the alignment
%
%             Other options:
%
%             'verbose=1' Set verbose level, 0 for quiet. 1 for normal. 2 for detailed info.
function retval = mlrVol(filename,varargin)

% check arguments
if nargin < 1
  help mlrVol
  return
end

global gSystem;

if isstr(filename)
  % init the system
  [sysNum otherFilenames] = initSystem(varargin);

  % load the volume
  if ~loadVolume(filename,sysNum),closeHandler(sysNum);return,end

  % load the other volumes
  for i = 1:length(otherFilenames)
    loadVolume(otherFilenames{i},sysNum);
  end

  % display the volume
  refreshDisplay(sysNum);
  
  % display controls
  if gSystem{sysNum}.displayControls,displayControls(sysNum);end
else
  % command argument
  switch (filename)
   case 1
    mouseDownHandler(varargin{1});
   case 2
    mouseMoveHandler(varargin{1});
   case 3
    textHandler(varargin{1},varargin{2});
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

uniqueFigs = unique(gSystem{sysNum}.fig);
% close the figures
for i = 1:length(uniqueFigs)
  % get the location of the figure - note that if
  % we ever have multiple figures than we will need to change
  % this code here and in initSystem to have more figlocs saved
  mrSetFigLoc('mlrVol',get(uniqueFigs(i),'Position'));

  % close the figure
  delete(uniqueFigs(i));
end

% remove the variable
gSystem{sysNum} = [];


%%%%%%%%%%%%%%%%%%%%%%%
%%   buttonHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function buttonHandler(textNum,sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% display controls
if textNum < 0
  switch textNum
   case -1
    closeHandler(sysNum);
   case -2
    displayControls(sysNum);
  end
  return
end

% if we are already running the animation, then turn it off
if gSystem{sysNum}.animating == textNum
  gSystem{sysNum}.animating = 0;
  return
end

% get number of dimensions
nDim = vol.h.nDim;

% start the animation
gSystem{sysNum}.animating = textNum;
coord = vol.coord;

% loop that runs the animation
while(gSystem{sysNum}.animating)
  % increment coord, go forward for an axis that moves in
  % the positive direction 
  if vol.axisDir(textNum) == 1
    coord(textNum) = coord(textNum)+1;
    if coord(textNum) > vol.h.dim(textNum)
      coord(textNum) = 1;
    end
    % and negatively otherwise
  else
    coord(textNum) = coord(textNum)-1;
    if coord(textNum) < 1
      coord(textNum) = vol.h.dim(textNum);
    end
  end
  % set the volume coord
  setVolCoord(sysNum,1,coord);
  % set the text boxes
  for iCoord = 1:nDim
    set(gSystem{sysNum}.hCoordTextbox(iCoord),'String',coord(iCoord));
  end
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(coord(1),coord(2),coord(3),coord(4),coord(5)));
  % and refresh
  refreshDisplay(sysNum);
end
%%%%%%%%%%%%%%%%%%%%%%%
%%   incdecHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function incdecHandler(textNum,incdec,sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% turn off animation
gSystem{sysNum}.animating = false;

% get the textbox we are updating. If textNum is less
% than zero we are updating the extra coord textboxes
if textNum > 0
  hTextbox = gSystem{sysNum}.hCoordTextbox(textNum);
else
  hTextbox = gSystem{sysNum}.hExtraCoordTextbox(-textNum);
end
  
% inc or dec the text box
val = str2num(get(hTextbox,'String'));
val = val+incdec;

% if we are doing an extra coordinate than transform
if (textNum < 0)
  set(hTextbox,'String',val);
  textHandler(sysNum,-1);
elseif (val>=1) && (val<=vol.h.dim(textNum))
  % if we have made a valid change then set it
  set(hTextbox,'String',val);
  % now refresh
  textHandler(sysNum,1);
end

%%%%%%%%%%%%%%%%%%%%%
%%   textHandler   %%
%%%%%%%%%%%%%%%%%%%%%
function textHandler(sysNum,textboxNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% turn off animation
gSystem{sysNum}.animating = 0;

% textboxNums that are negative are the extra coordinates
if textboxNum < 0
  % set coordinate according to extra
  setCoordToExtra(sysNum);
end

% get number of dimensions
nDim = gSystem{sysNum}.vols(1).h.nDim;

% get coordinate
for iCoord = 1:nDim
  coord(iCoord) = str2num(get(gSystem{sysNum}.hCoordTextbox(iCoord),'String'));
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
setVolCoord(sysNum,1,coord);

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
  updateExtraCoords(sysNum);
  % nDims hard coded to 5 here
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(vol.coord(1),vol.coord(2),vol.coord(3),vol.coord(4),vol.coord(5)));
else
  set(gcf,'Pointer','fullcrosshair');
  set(gSystem{sysNum}.hCoordTextbox(1),'String',coord(1));
  set(gSystem{sysNum}.hCoordTextbox(2),'String',coord(2));
  set(gSystem{sysNum}.hCoordTextbox(3),'String',coord(3));
  updateExtraCoords(sysNum);
  % nDims hard coded to 5 here
  set(gSystem{sysNum}.hValTextbox(1),'String',vol.data(coord(1),coord(2),coord(3),coord(4),coord(5)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   updateExtraCoords   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateExtraCoords(sysNum)

global gSystem
if gSystem{sysNum}.vols(1).altXforms.n > 0
  % get the current coord
  coord(1) = str2num(get(gSystem{sysNum}.hCoordTextbox(1),'String'));
  coord(2) = str2num(get(gSystem{sysNum}.hCoordTextbox(2),'String'));
  coord(3) = str2num(get(gSystem{sysNum}.hCoordTextbox(3),'String'));
  coord(4) = 1;
  coord = coord(:);
  
  % get the extra coordinate xform
  xform = gSystem{sysNum}.vols(1).altXforms.xforms{gSystem{sysNum}.vols(1).altXforms.currentXform};
  
  % convert
  coord = xform * coord;
  
  % and display
  set(gSystem{sysNum}.hExtraCoordTextbox(1),'String',sprintf('%0.1f',coord(1)));
  set(gSystem{sysNum}.hExtraCoordTextbox(2),'String',sprintf('%0.1f',coord(2)));
  set(gSystem{sysNum}.hExtraCoordTextbox(3),'String',sprintf('%0.1f',coord(3)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   setCoordToExtra   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function setCoordToExtra(sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(1);

% get coordinates form extra textbox
coord(1) = str2num(get(gSystem{sysNum}.hExtraCoordTextbox(1),'String'));
coord(2) = str2num(get(gSystem{sysNum}.hExtraCoordTextbox(2),'String'));
coord(3) = str2num(get(gSystem{sysNum}.hExtraCoordTextbox(3),'String'));
coord(4) = 1;
coord = coord(:);

% get the extra coordinate xform
xform = gSystem{sysNum}.vols(1).altXforms.xforms{gSystem{sysNum}.vols(1).altXforms.currentXform};
  
% convert
coord = round(inv(xform) * coord);

% check to see if we have valid coordinates
if all(coord(1:3) > 0) && all(coord(1:3) <= vol.h.dim(1:3))
  % if so, then display
  set(gSystem{sysNum}.hCoordTextbox(1),'String',sprintf('%0.0f',coord(1)));
  set(gSystem{sysNum}.hCoordTextbox(2),'String',sprintf('%0.0f',coord(2)));
  set(gSystem{sysNum}.hCoordTextbox(3),'String',sprintf('%0.0f',coord(3)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mouseDownHandler    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(sysNum)

set(gcf,'Pointer','arrow');
global gSystem;
vol = gSystem{sysNum}.vols(1);

% stop any ongoing animation
gSystem{sysNum}.animating = 0;

% get the mouse pos
coord = getMouseCoord(sysNum);
if isempty(coord),return,end

% and set the volume coordinate
setVolCoord(sysNum,1,coord);

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
  a = subplot(gSystem{sysNum}.subplotRows(1),gSystem{sysNum}.subplotCols(1),subplotNum);
else
  return
end
viewNum = find(gSystem{sysNum}.a(1,1:3) == a);
if isempty(viewNum),return,end

% get pointer loc
pointerLoc = get(a,'CurrentPoint');
pointerLoc = round(pointerLoc(1,2:-1:1));

% check transpose to see which coordinate is which
% note that matlab displays matrices in a transposed fashion
if ~vol.transpose(viewNum)
  pointerX = pointerLoc(1);
  pointerY = pointerLoc(2);
else
  pointerX = pointerLoc(2);
  pointerY = pointerLoc(1);
end

% check image boundaries 
if (pointerX < 1) | (pointerX > vol.h.dim(vol.xDim(viewNum))),return,end
if (pointerY < 1) | (pointerY > vol.h.dim(vol.yDim(viewNum))),return,end


% get the coordinate. Note that we have to be careful both of whether the axis
% is flipped AND whether there was a transpose (since the y-axis goes in
% the opposite direction as x for matlab displayed images
if ((vol.axisDir(vol.xDim(viewNum)) == 1) && (vol.transpose(viewNum))) || ((vol.axisDir(vol.xDim(viewNum)) == -1) && (~vol.transpose(viewNum)))
  coord(vol.xDim(viewNum)) = pointerX;
else
  coord(vol.xDim(viewNum)) = vol.h.dim(vol.xDim(viewNum)) - pointerX + 1;
end
% note that for the y-dimension, matlab's axis are flipped
if ((vol.axisDir(vol.yDim(viewNum)) == -1) && (vol.transpose(viewNum))) || ((vol.axisDir(vol.yDim(viewNum)) == 1) && (~vol.transpose(viewNum)))
  coord(vol.yDim(viewNum)) = pointerY;
else
  coord(vol.yDim(viewNum)) = vol.h.dim(vol.yDim(viewNum)) - pointerY + 1;
end
coord(vol.viewDim(viewNum)) = vol.coord(vol.viewDim(viewNum));

% nDims hard coded to 5 here
coord(end+1:5) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%    refreshDisplay    %
%%%%%%%%%%%%%%%%%%%%%%%%
function refreshDisplay(sysNum)

global gSystem;

f = gcf;
% set the figure
figure(gSystem{sysNum}.vols(1).fig(1));

% all other volumes are tethered to the first volume
for iVol = 1:gSystem{sysNum}.n
  % if the volume is tethered then set its coordinates according
  % to whatever volume it is tethered to.
  if gSystem{sysNum}.vols(iVol).tethered
    setTetheredCoord(sysNum,iVol,gSystem{sysNum}.vols(iVol).tethered);
  end

  % display the coordinate
  dispVolume(iVol,sysNum)

  % set that we have displayed this coordinate
  gSystem{sysNum}.vols(iVol).curCoord = gSystem{sysNum}.vols(iVol).coord;
end

updateExtraCoords(sysNum);

drawnow;
figure(f);


%%%%%%%%%%%%%%%%%%%%
%    dispVolume    %
%%%%%%%%%%%%%%%%%%%%
function dispVolume(iVol,sysNum)

global gSystem;
vol = gSystem{sysNum}.vols(iVol);

% view labels
viewLabel = {'Sagittal','Coronal','Axial'};

for iView = 1:3
  % see if we need to redisplay, first condition is whether we have changed the volume or other dimenson
  % the second condition is for primary volumes wether we have updated the coordinates being displayed
  % in the view, and the third condition is for tethered volumes - whether we have change the coordinates
  % for the volume we are tethered to
  if (~isequal(vol.curCoord(4:end),vol.coord(4:end)) || ...
      (~vol.tethered && ~isequal(vol.curCoord(vol.viewDim(iView)),vol.coord(vol.viewDim(iView)))) || ...
      (vol.tethered && ~isequal(vol.curCoord(gSystem{sysNum}.vols(vol.tethered).viewDim(iView)),vol.coord(gSystem{sysNum}.vols(vol.tethered).viewDim(iView)))))
    % get the slice
    % nDims hard coded to 5 here
    if vol.tethered
      % tethered volumes have their display slices precomputed (by interpolation) in 
      % setTetheredCoord
      dispSlice = vol.dispSlice{iView};

      % get the axis that we are tethered to (so that we can draw the overlay)
      aTether = subplot(gSystem{sysNum}.subplotRows(iView),gSystem{sysNum}.subplotCols(iView),gSystem{sysNum}.vols(vol.tethered).subplotNum(iView));

      % the transpose and axis directions need to be taken from the volume this is tethered to
      % prepare image for display
      [dispOverlaySlice xLabelStr yLabelStr] = prepareImageForDisplay(dispSlice,vol,iView,gSystem{sysNum}.vols(vol.tethered).transpose,gSystem{sysNum}.vols(vol.tethered).axisDir);
      
      % if we are not displaying the interpolated image in the
      % second row, then we have to prepare the image that
      % corresponds to the same location
      if gSystem{sysNum}.displayInterpolated
	dispSlice = dispOverlaySlice;
      else
	% need to get the coordinate of the tethered to volume
	% in these coordinates
	dispSlice = getMatchingSlice(sysNum,vol,iView);
	[dispSlice xLabelStr yLabelStr] = prepareImageForDisplay(dispSlice,vol,iView,vol.transpose,vol.axisDir);
      end
    else
      % otherwise, grab the data for this image
      dispSlice = squeeze(vol.data(vol.viewIndexes{iView,1},vol.viewIndexes{iView,2},vol.viewIndexes{iView,3},vol.coord(4),vol.coord(5)));

      % prepare image for display
      [dispSlice xLabelStr yLabelStr] = prepareImageForDisplay(dispSlice,vol,iView,vol.transpose,vol.axisDir);
    end

    % get the correct axis to draw into
    a = subplot(gSystem{sysNum}.subplotRows(iView),gSystem{sysNum}.subplotCols(iView),vol.subplotNum(iView));

    % and display the image
    cla(a);
    if isempty(dispSlice)
      axis off;
    else
      subimage(dispSlice,gray(256));
      % turn off labels
      set(a,'XTickLabel','');
      set(a,'YTickLabel','');
      % and put on what axis we have
      xlabel(a,xLabelStr);
      ylabel(a,yLabelStr);
    end
    % set title
    titleStr = sprintf('%s (%s)',viewLabel{iView},vol.h.filename);
    h = title(titleStr,'Interpreter','none');
    
    % and display the overlay
    if vol.tethered
      axes(aTether);
      hold on;
      gSystem{sysNum}.vols(iVol).overlay(iView) = subimage(dispOverlaySlice,hot(256));
      gSystem{sysNum}.vols(iVol).overlayMask{iView} = ~isnan(dispOverlaySlice(:));
      alphaData = zeros(size(dispOverlaySlice));
      if gSystem{sysNum}.overlayToggleState
	alphaData(gSystem{sysNum}.vols(iVol).overlayMask{iView}) = gSystem{sysNum}.overlayAlpha;
      end
      set(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData',alphaData);
      hold off;
    end
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getMatchingSlice   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispSlice = getMatchingSlice(sysNum,vol,iView)

global gSystem;

% get coordinated of volume we are tethered to
coord = gSystem{sysNum}.vols(vol.tethered).coord;
coord(4) = 1;coord = coord(1:4);coord = coord(:);

% convert to this volume coordinates
coord = round(applySystemXform(sysNum,vol.xform)*coord);
coord = coord(1:3);

% FIX: this hardcodes the volume number to 1
coord(end+1:5) = 1;

% check to see if coordinates are in volume
if any(coord(1:3)<1) || any(coord(1:3)>vol.h.dim(1:3))
  dispSlice = [];
  return
end

% if not set the volume coordinates
setVolCoord(sysNum,vol.volnum,coord);
vol = gSystem{sysNum}.vols(vol.volnum);

% and get display slice
dispSlice = squeeze(vol.data(vol.viewIndexes{iView,1},vol.viewIndexes{iView,2},vol.viewIndexes{iView,3},vol.coord(4),vol.coord(5)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   prepareImageForDisplay   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dispSlice xLabelStr yLabelStr] = prepareImageForDisplay(dispSlice,vol,iView,transpose,axisDir)

if isempty(dispSlice)
  xLabelStr = '';yLabelStr = '';
  return
end
% clip
%  dispSlice = clipImage(dispSlice);

% make into image with index values
minDispSlice = min(dispSlice(:));
maxDispSlice = max(dispSlice(:));
dispSlice = ceil(256*(dispSlice-minDispSlice)/(maxDispSlice-minDispSlice));

% do transpose if necessary
if transpose(iView)
  dispSlice = dispSlice';
  xLabelStr = vol.dispAxisLabels{vol.xDim(iView)};
  yLabelStr = vol.dispAxisLabels{vol.yDim(iView)};
else
  xLabelStr = vol.dispAxisLabels{vol.yDim(iView)};
  yLabelStr = vol.dispAxisLabels{vol.xDim(iView)};
end

% flip axis if necessary (note that Matlab shows the x axis as - to +
% and the y -axis in the opposite orientation, so we treat the x and y differently)
if axisDir(vol.xDim(iView)) == -1,dispSlice = fliplr(dispSlice);end
if axisDir(vol.yDim(iView)) == 1,dispSlice = flipud(dispSlice);end


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
function tf = loadVolume(filename,sysNum)

tf = false;
global gSystem;

[d h] = mlrLoadImage(filename);
if isempty(d),return,end

% updated number of volumes
gSystem{sysNum}.n = gSystem{sysNum}.n+1;
n = gSystem{sysNum}.n;

% set volume number
gSystem{sysNum}.vols(n).volnum = n;

% update the fields in vols
gSystem{sysNum}.vols(n).data = d;
gSystem{sysNum}.vols(n).h = h;

% set which figure numbers this volume will display into
gSystem{sysNum}.vols(n).fig(1:3) = gSystem{sysNum}.fig(1:3);

% set which subplot
gSystem{sysNum}.vols(n).subplotNum(1:3) = (1:3)+(n-1)*3;

% compute magnet directions of each axis based on qform and then
% choose which axis will be displayed in what figure based on this
% information
if h.qform_code && ~gSystem{sysNum}.imageOrientation
  [gSystem{sysNum}.vols(n).axisLabels gSystem{sysNum}.vols(n).axisMapping gSystem{sysNum}.vols(n).axisDir] = getAxisLabels(h.qform44);
  gSystem{sysNum}.vols(n).viewDim(1:3) = gSystem{sysNum}.vols(n).axisMapping;
else
  gSystem{sysNum}.vols(n).axisLabels = [];
  % default to assuming LPI orientation
  gSystem{sysNum}.vols(n).axisMapping = [1 2 3];
  gSystem{sysNum}.vols(n).axisDir = [1 1 1];
  % no information about what axis is what, so just show each axis in order
  gSystem{sysNum}.vols(n).viewDim(1:3) = 1:3;
end

% for convenience set which dimension is x and y for each of the views
% given the viewDim info set above and also info about transposing and flipping
desiredAxis = {[2 3],[1 3],[1 2]};
for iView = 1:3
  % get what the other axis are
  otherDims = setdiff(1:3,gSystem{sysNum}.vols(n).viewDim(iView));
  gSystem{sysNum}.vols(n).xDim(iView) = otherDims(1);
  gSystem{sysNum}.vols(n).yDim(iView) = otherDims(2);
  % next figure out the apporpriate transpose and flips needed to
  % show the images in standard LPI (i.e. the first view will
  % be sagittal with nose to right, second view willbe coronal and
  % third view will be axial) all with left is left, right is right
  if gSystem{sysNum}.vols(n).axisMapping(otherDims(1)) == desiredAxis{iView}(1)
    % transpose when the desired axis is the same (since
    % images are displayed as y/x by matlab)
    gSystem{sysNum}.vols(n).transpose(iView) = 1;
  else
    gSystem{sysNum}.vols(n).transpose(iView) = 0;
  end
  % now make axis labels that can be used for displaying. Note that
  % if the axis are flipped, dispSlice will unflip them so we
  % have to reverse the order of the labels provided by getAxisLabels
  axisLabel = {'X','Y','Z'};
  if isempty(gSystem{sysNum}.vols(n).axisLabels)
    % no transform, so just use simple X,Y,Z labels
    gSystem{sysNum}.vols(n).dispAxisLabels{iView} = axisLabel{iView};
  else
    % check axis direction
    if gSystem{sysNum}.vols(n).axisDir == -1
      % and make reversed labels
      gSystem{sysNum}.vols(n).dispAxisLabels{iView} = sprintf('%s <- %s -> %s',gSystem{sysNum}.vols(n).axisLabels{iView}{2},axisLabel{iView},gSystem{sysNum}.vols(n).axisLabels{iView}{1});
    else
      % and make normal labels
      gSystem{sysNum}.vols(n).dispAxisLabels{iView} = sprintf('%s <- %s -> %s',gSystem{sysNum}.vols(n).axisLabels{iView}{1},axisLabel{iView},gSystem{sysNum}.vols(n).axisLabels{iView}{2});
    end
  end
end

% update the coordinate (start displaying in middle of volume)
coord = round(h.dim/2);
% nDims hard coded to 5 here
coord(end+1:5) = 1;
setVolCoord(sysNum,n,coord);

% now this sets the indexes from the volume for which the
% image will be displayed
for iView = 1:3
  for jAxis = 1:3
    if jAxis ~= gSystem{sysNum}.vols(n).viewDim(iView)
      gSystem{sysNum}.vols(n).viewIndexes{iView,jAxis} = 1:gSystem{sysNum}.vols(n).h.dim(jAxis);
    else
      gSystem{sysNum}.vols(n).viewIndexes{iView,jAxis} = gSystem{sysNum}.vols(n).coord(jAxis);
    end
  end
end

% set the current coordinates for the first time
gSystem{sysNum}.vols(n).curCoord = [nan nan nan];

% print out the header of the image
dispHeaderInfo(sysNum,gSystem{sysNum}.vols(n));

% see what other transforms we have
gSystem{sysNum}.vols(n).altXforms.shortNames = {};
gSystem{sysNum}.vols(n).altXforms.names = {};
gSystem{sysNum}.vols(n).altXforms.xforms = {};
gSystem{sysNum}.vols(n).altXforms.n = 0;
gSystem{sysNum}.vols(n).altXforms.currentXform = [];
if ~isempty(h.vol2tal)
  gSystem{sysNum}.vols(n).altXforms.names{end+1} = 'Talairach';
  gSystem{sysNum}.vols(n).altXforms.shortNames{end+1} = 'Tal';
  gSystem{sysNum}.vols(n).altXforms.xforms{end+1} = h.vol2tal;
  gSystem{sysNum}.vols(n).altXforms.n = gSystem{sysNum}.vols(n).altXforms.n+1;
end
if ~isempty(h.vol2mag)
  gSystem{sysNum}.vols(n).altXforms.names{end+1} = 'Canonical';
  gSystem{sysNum}.vols(n).altXforms.shortNames{end+1} = 'Base';
  gSystem{sysNum}.vols(n).altXforms.xforms{end+1} = h.vol2mag;
  gSystem{sysNum}.vols(n).altXforms.n = gSystem{sysNum}.vols(n).altXforms.n+1;
end
if h.qform_code
  gSystem{sysNum}.vols(n).altXforms.names{end+1} = 'Qform';
  gSystem{sysNum}.vols(n).altXforms.shortNames{end+1} = 'Base';
  gSystem{sysNum}.vols(n).altXforms.xforms{end+1} = h.qform44;
  gSystem{sysNum}.vols(n).altXforms.n = gSystem{sysNum}.vols(n).altXforms.n+1;
end

% set the text boxes
if n == 1
  names = {'X','Y','Z','Slice','Receiver'};
  for i = 1:h.nDim
    % make the inc/dec textboxes
    gSystem{sysNum}.hCoordTextbox(i) = makeTextboxIncDec(sysNum,coord(i),i,1,i);
    % get name for button
    if i < length(names)
      name = names{i};
    else 
      name = sprintf('dim %i',i);
    end
    % make the button
    gSystem{sysNum}.hButton(i) = makeButton(sysNum,name,i,2,i);
  end
  % make a text box for the value of the voxel
  makeTextbox(sysNum,'Value',2,h.nDim+1);
  gSystem{sysNum}.hValTextbox(1) = makeTextbox(sysNum,'',1,h.nDim+1);

  % make the extra coordinates text boxes
  if gSystem{sysNum}.vols(n).altXforms.n > 0
    gSystem{sysNum}.vols(n).altXforms.currentXform = 1;
    shortName = gSystem{sysNum}.vols(n).altXforms.shortNames{gSystem{sysNum}.vols(n).altXforms.currentXform};
    gSystem{sysNum}.hExtraCoordTitle(1) = makeTextbox(sysNum,sprintf('%s X',shortName),2,h.nDim+2);
    gSystem{sysNum}.hExtraCoordTextbox(1) = makeTextboxIncDec(sysNum,'',-1,1,h.nDim+2);
    gSystem{sysNum}.hExtraCoordTitle(2) = makeTextbox(sysNum,sprintf('%s Y',shortName),2,h.nDim+3);
    gSystem{sysNum}.hExtraCoordTextbox(2) = makeTextboxIncDec(sysNum,'',-2,1,h.nDim+3);
    gSystem{sysNum}.hExtraCoordTitle(3) = makeTextbox(sysNum,sprintf('%s Z',shortName),2,h.nDim+4);
    gSystem{sysNum}.hExtraCoordTextbox(3) = makeTextboxIncDec(sysNum,'',-3,1,h.nDim+4);
  end
end

% add another row for displaying the image
if n > 1
  makeButton(sysNum,'Controls',-2,4,1);
  % mark that the volume display is "tethered" to the first volume
  gSystem{sysNum}.vols(n).tethered = 1;
  % make a cache for storing images
  gSystem{sysNum}.vols(n).c = mrCache('init',2*max(gSystem{sysNum}.vols(n).h.dim(1:3)));
  % set the initial transform
  gSystem{sysNum}.vols(n).xform = getVol2vol(sysNum,gSystem{sysNum}.vols(n),gSystem{sysNum}.vols(1));
  % setup subplotRows and subplotCols
  gSystem{sysNum}.subplotRows = [n n n];
  gSystem{sysNum}.subplotCols = [3 3 3];
  % redo the subplots
  for i = 1:n
    for j = 1:3
      gSystem{sysNum}.a(i,j) = subplot(n,3,j+(i-1)*3);
      cla;
      axis off;
    end
  end
  % set to display controls
%  gSystem{sysNum}.displayControls = true;
else
  % first volume is displayed independently
  gSystem{sysNum}.vols(n).tethered = 0;
  % make close button
  makeButton(sysNum,'Close',-1,3,1);
end

tf = true;

%%%%%%%%%%%%%%%%%%%%%%%
%    getAxisLabels    %
%%%%%%%%%%%%%%%%%%%%%%%
function [axisLabels axisMapping axisDirection] = getAxisLabels(qform44)

% This function uses the qform to determine in which direction
% each axis of the image goes. axisLabels are readable labels 
% and axisDirs are a vector for the closest pointing direction
% for each axis

for axisNum = 1:3
  % get the vector of the axis in the image that we want to label
  axisVector = zeros(3,1);
  axisVector(axisNum) = 1;

  % find out which direction in magnet space the axis vector goes
  axisDirVector = qform44(1:3,1:3)*axisVector;

  % normalize to unit length
  axisDirVector = axisDirVector ./ sqrt(sum(axisDirVector.^2));

  % these are the magnet cardinal axis and their names
  cardinalAxisDirs = {[1 0 0],[0 1 0],[0 0 1],[-1 0 0],[0 -1 0],[0 0 -1]};
  cardinalAxisLabels = {'right','anterior','superior','left','posterior','inferior'};

  % now get the angle of this axisDirVector with
  % each of the magnet cardinal axis. Note that
  % we are computing for both directions of the axis
  % for conveinence (i.e. left and right)
  for i = 1:length(cardinalAxisDirs)
    angles(i) = r2d(acos(dot(axisDirVector,cardinalAxisDirs{i})));
  end

  % sort the angles (remembering which axis they originally came from). Thus
  % sortedAxisNum contains an ordered list of which axis the vector is closest
  % too. The closest axis is sortedAxisNum(1) and the farthest axis is sortedAxisNum(6)
  [angles sortedAxisNum] = sort(angles);

  % get the closest axis direction (i.e. the one with the smallest angle which
  % is the first in the list of sortedAxisNum)
  axisDirs(axisNum,:) = cardinalAxisDirs{sortedAxisNum(1)}';
  
  % if the closest angle is less than an arbitrary value than we will consider
  % the axis to be a pure direction - if not, we will label
  % with a combination of the two closest axis.
  if angles(1) < 5
    axisLabels{axisNum} = {cardinalAxisLabels{sortedAxisNum(6)} cardinalAxisLabels{sortedAxisNum(1)}};
  else
    axisLabels{axisNum} = {sprintf('%s/%s',cardinalAxisLabels{sortedAxisNum(6)},cardinalAxisLabels{sortedAxisNum(5)}) sprintf('%s/%s',cardinalAxisLabels{sortedAxisNum(1)},cardinalAxisLabels{sortedAxisNum(2)})};
  end

end

% convert the axisDirs to axisMapping (i.e. what axis in the magnet each axis in the image
% corresponds to and in which direction it points
[axisMapping row axisDirection] = find(axisDirs);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applySystemXform    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function xform = applySystemXform(sysNum,xform)

global gSystem;

xform = inv(shiftOriginXform) * xform * gSystem{sysNum}.xform * shiftOriginXform;

%%%%%%%%%%%%%%%%%%%%
%    getVol2vol    %
%%%%%%%%%%%%%%%%%%%%
function vol2vol = getVol2vol(sysNum,vol1,vol2)

global gSystem;
verbose = gSystem{sysNum}.verbose;

if vol1.h.sform_code && vol2.h.sform_code
  vol2vol = inv(vol1.h.sform44) * vol2.h.sform44;
  if verbose
    dispHeader('Aliging using sform');
    disp(sprintf('%s',mrnum2str(vol2vol,'compact=0')))
    dispHeader;
  end
%  vol2vol = inv(shiftOriginXform) * xform * shiftOriginXform;
elseif vol1.h.qform_code && vol2.h.qform_code
  vol2vol = inv(vol1.h.qform44) * vol2.h.qform44;
  if verbose
    dispHeader('Aliging using qform');
    disp(sprintf('%s',mrnum2str(vol2vol,'compact=0')))
    dispHeader;
  end
else
  vol2vol = eye(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispHeaderInfo   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function dispHeaderInfo(sysNum,vol)

global gSystem;

if gSystem{sysNum}.verbose
  dispHeader(vol.h.filename);
  disp(sprintf('type: %s (%s)',vol.h.type,vol.h.ext));
  disp(sprintf('dim: [%s]',mrnum2str(vol.h.dim(:)','sigfigs=0')));
  disp(sprintf('pixdim: [%s]',mrnum2str(vol.h.pixdim(:)')));
  disp(sprintf('qform_code: %i',vol.h.qform_code));
  disp(sprintf('qform:'));
  disp(sprintf('%s',mrnum2str(vol.h.qform44,'compact=0','sigfigs=-1')));
  disp(sprintf('sform_code: %i',vol.h.sform_code));
  disp(sprintf('sform:'));
  disp(sprintf('%s',mrnum2str(vol.h.sform44,'compact=0','sigfigs=-1')));

  % display axis information
  if ~isempty(vol.axisLabels)
    cardinalAxisLabels = {'X','Y','Z'};
    disp(sprintf('Volume orientation is: %s%s%s',upper(vol.axisLabels{1}{1}(1)),upper(vol.axisLabels{2}{1}(1)),upper(vol.axisLabels{3}{1}(1))));
    for axisNum = 1:3
      disp(sprintf('Axis %s goes from %s to %s',cardinalAxisLabels{axisNum},vol.axisLabels{axisNum}{1},vol.axisLabels{axisNum}{2}));
    end
  end
  
  % if there is a talInfo field, display that
  if isfield(vol.h,'base') && isfield(vol.h.base,'talInfo')
    disp(sprintf('AC: [%s]',mrnum2str(vol.h.base.talInfo.AC,'compact=1','sigfigs=0')));
    disp(sprintf('PC: [%s]',mrnum2str(vol.h.base.talInfo.PC,'compact=1','sigfigs=0')));
    disp(sprintf('SAC: [%s]',mrnum2str(vol.h.base.talInfo.SAC,'compact=1','sigfigs=0')));
    disp(sprintf('IAC: [%s]',mrnum2str(vol.h.base.talInfo.IAC,'compact=1','sigfigs=0')));
    disp(sprintf('PPC: [%s]',mrnum2str(vol.h.base.talInfo.PPC,'compact=1','sigfigs=0')));
    disp(sprintf('AAC: [%s]',mrnum2str(vol.h.base.talInfo.AAC,'compact=1','sigfigs=0')));
    disp(sprintf('LAC: [%s]',mrnum2str(vol.h.base.talInfo.LAC,'compact=1','sigfigs=0')));
    disp(sprintf('RAC: [%s]',mrnum2str(vol.h.base.talInfo.RAC,'compact=1','sigfigs=0')));
  end
end

% display detailed header information
if gSystem{sysNum}.verbose>1
  hdrFields = fieldnames(vol.h.hdr);
  for iField = 1:length(hdrFields)
    val = vol.h.hdr.(hdrFields{iField});
    if isnumeric(val)
      if (size(val,1) == 1) && (size(val,2) == 1)
	disp(sprintf('%s: %s',hdrFields{iField},mrnum2str(val)));
      elseif size(val,1) == 1
	disp(sprintf('%s: [%s]',hdrFields{iField},mrnum2str(val)));
      elseif size(val,2) == 1
	disp(sprintf('%s: [%s]',hdrFields{iField},mrnum2str(val')));
      else
	disp(sprintf('%s:\n%s',hdrFields{iField},mrnum2str(val,'compact=0')));
      end
    elseif isstr(val)
      disp(sprintf('%s: %s',hdrFields{iField},val));
    elseif isempty(val)
      disp(sprintf('%s: []',hdrFields{iField}));
    elseif isstruct(val)
      disp(sprintf('%s: struct',hdrFields{iField}));
    else
      disp(sprintf('%s: Unknown type',hdrFields{iField}));
    end      
  end
end

if gSystem{sysNum}.verbose==1,dispHeader;end


%%%%%%%%%%%%%%%%%%%
%    mrnum2str    %
%%%%%%%%%%%%%%%%%%%
function str = mrnum2str(num,arg1,arg2)

% just a wrapper function so that we can use mynum2str (which is in my matlab directory not in mrTools)
if exist('mynum2str')==2
  switch (nargin)
   case 1
    str = mynum2str(num);
   case 2
    str = mynum2str(num,arg1);
   case 3
    str = mynum2str(num,arg1,arg2);
  end
else
  % otherwise use matlabs normal function
  str = num2str(num(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setTetheredCoord    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function setTetheredCoord(sysNum,iVol,iRefvol)

global gSystem;

% get coord we are tethering to
coord = gSystem{sysNum}.vols(iRefvol).coord;

% see if we have already updated
if isequal(gSystem{sysNum}.vols(iVol).curCoord,coord)
  return
end

% look in cache for precalculate
dispSlice = mrCache('find',gSystem{sysNum}.vols(iVol).c,mrnum2str(coord));
if ~isempty(dispSlice)
  gSystem{sysNum}.vols(iVol).dispSlice = dispSlice;
else
  % create an image
  for iView = 1:3
    % get the coordinates of the reference volume
    [x y z] = ndgrid(gSystem{sysNum}.vols(iRefvol).viewIndexes{iView,1},gSystem{sysNum}.vols(iRefvol).viewIndexes{iView,2},gSystem{sysNum}.vols(iRefvol).viewIndexes{iView,3});
    s = size(x);
    % make into coords for multiplying
    coords = [x(:) y(:) z(:)]';
    coords(4,:) = 1;
    % convert from the reference volume coordinates to our coordinates
    coords = applySystemXform(sysNum,gSystem{sysNum}.vols(iVol).xform) * coords;
    x = reshape(coords(1,:),s);
    y = reshape(coords(2,:),s);
    z = reshape(coords(3,:),s);
    % get the interpolated image (note that interp3 needs to have y and x swaped to work correctly here)
    gSystem{sysNum}.vols(iVol).dispSlice{iView} = squeeze(interp3(gSystem{sysNum}.vols(iVol).data,y,x,z,gSystem{sysNum}.interpMethod,nan));
  end
  % save in cache
  gSystem{sysNum}.vols(iVol).c = mrCache('add',gSystem{sysNum}.vols(iVol).c,mrnum2str(coord),gSystem{sysNum}.vols(iVol).dispSlice);
end

% set the coordinate so that we update correctly
gSystem{sysNum}.vols(iVol).coord = coord;
%%%%%%%%%%%%%%%%%%%%%
%    setVolCoord    %
%%%%%%%%%%%%%%%%%%%%%
function setVolCoord(sysNum,iVol,coord)

global gSystem;

% set the current x,y,z coordinate
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
function [sysNum otherFilenames] = initSystem(args)

global gSystem;
if isempty(gSystem)
  sysNum = 1;
else
  sysNum = length(gSystem)+1;
end

% get any aditional filenames
otherFilenames = {};
otherArgs = {};
for i = 1:length(args)
  if isstr(args{i}) && isempty(strfind(args{i},'='))
    otherFilenames{end+1} = args{i};
    otherArgs = {args{i+1:end}};
  else
    otherArgs = {args{i:end}};
    break;
  end
end

% number of loaded volumes
gSystem{sysNum}.n = 0;

% parse args here when we have settings
imageOrientation = [];verbose = [];
getArgs(otherArgs,{'imageOrientation=0','verbose=1'});
gSystem{sysNum}.imageOrientation = imageOrientation;
gSystem{sysNum}.verbose = verbose;


% defaults for button sizes
gSystem{sysNum}.buttonWidth = 100;
gSystem{sysNum}.buttonWidthMargin = 20;
gSystem{sysNum}.buttonHeightMargin = 2;
gSystem{sysNum}.buttonHeight = 25;
gSystem{sysNum}.buttonLeftMargin = 10;
gSystem{sysNum}.buttonBottomMargin = 10;

% get location of figure
figloc = mrGetFigLoc('mlrVol');

% open the fig
gSystem{sysNum}.fig(1) = figure;
if ~isempty(figloc)
  set(gSystem{sysNum}.fig(1),'Position',figloc);
end
gSystem{sysNum}.fig(2) = gSystem{sysNum}.fig(1);
gSystem{sysNum}.fig(3) = gSystem{sysNum}.fig(1);
clf;

% set the mouse functions
set(gSystem{sysNum}.fig(1),'WindowButtonDownFcn',sprintf('mlrVol(1,%i)',sysNum));
set(gSystem{sysNum}.fig(1),'WindowButtonMotionFcn',sprintf('mlrVol(2,%i)',sysNum));
set(gSystem{sysNum}.fig(1),'CloseRequestFcn',sprintf('mlrVol(6,%i)',sysNum));

% set up the axis
gSystem{sysNum}.a(1) = subplot(1,3,1);cla;axis off;
gSystem{sysNum}.a(2) = subplot(1,3,2);cla;axis off;
gSystem{sysNum}.a(3) = subplot(1,3,3);cla;axis off;
gSystem{sysNum}.subplotRows = [1 1 1];
gSystem{sysNum}.subplotCols = [3 3 3];

% no animation is running
gSystem{sysNum}.animating = false;

% get interp method
gSystem{sysNum}.interpMethod = mrGetPref('interpMethod');

% alpha for overlay
gSystem{sysNum}.overlayAlpha = 1;

% default not to display controls
gSystem{sysNum}.displayControls = false;

% default system transform
gSystem{sysNum}.xform = eye(4);

% overlay toggle state starts as on
gSystem{sysNum}.overlayToggleState = 1;

% set up system params
gSystem{sysNum}.xformParams.shiftX = 0;
gSystem{sysNum}.xformParams.shiftY = 0;
gSystem{sysNum}.xformParams.shiftZ = 0;
gSystem{sysNum}.xformParams.rotateXY = 0;
gSystem{sysNum}.xformParams.rotateXZ = 0;
gSystem{sysNum}.xformParams.rotateYZ = 0;

% display the tethered volume interpolated to
% match the primary volume display
gSystem{sysNum}.displayInterpolated = true;

% update display
drawnow

%%%%%%%%%%%%%%%%%%%%%
%%   makeTextbox   %%
%%%%%%%%%%%%%%%%%%%%%
function h = makeTextbox(sysNum,displayString,rownum,colnum);

h = uicontrol('Style','text','String',displayString,'Position',getUIControlPos(sysNum,rownum,colnum,1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(3,%i)',sysNum));

%%%%%%%%%%%%%%%%%%%%%%%
%%   makePushbuton   %%
%%%%%%%%%%%%%%%%%%%%%%%
function h = makeButton(sysNum,displayString,textboxNum,rownum,colnum)

h = uicontrol('Style','pushbutton','String',displayString,'Position',getUIControlPos(sysNum,rownum,colnum,1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(5,%i,%i)',textboxNum,sysNum));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeTextboxIncDec   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h hIncDec] = makeTextboxIncDec(sysNum,displayString,textboxNum,rownum,colnum)

% make textbox
h = uicontrol('Style','edit','String',displayString,'Position',getUIControlPos(sysNum,rownum,colnum+.125,0.75),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(3,%i,%i)',sysNum,textboxNum));

% put up incdec buttons
hIncDec(1) = uicontrol('Style','pushbutton','String','<','Position',getUIControlPos(sysNum,rownum,colnum,.1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(4,%i,-1,%i)',textboxNum,sysNum));
hIncDec(2) = uicontrol('Style','pushbutton','String','>','Position',getUIControlPos(sysNum,rownum,colnum+.9,.1),'FontSize',10,'FontName','Helvetica','HorizontalAlignment','Center','Callback',sprintf('mlrVol(4,%i,1,%i)',textboxNum,sysNum));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(sysNum,rownum,colnum,uisize)

global gSystem;

% set this buttons width
thisButtonWidth = gSystem{sysNum}.buttonWidth*uisize;

% set the position for the button
pos(1) = (gSystem{sysNum}.buttonWidth+gSystem{sysNum}.buttonWidthMargin)*(floor(colnum)-1)+gSystem{sysNum}.buttonLeftMargin+(colnum-floor(colnum))*gSystem{sysNum}.buttonWidth;
pos(2) = gSystem{sysNum}.buttonBottomMargin + (gSystem{sysNum}.buttonHeight+gSystem{sysNum}.buttonHeightMargin)*(rownum-1);
pos(3) = thisButtonWidth;
pos(4) = gSystem{sysNum}.buttonHeight;

%%%%%%%%%%%%%%%%%%%%
%    dispHeader    %
%%%%%%%%%%%%%%%%%%%%
function retval = dispHeader(header,len,c)

% check arguments
if ~any(nargin == [0 1 2 3])
  help dispHeader
  return
end

% default header is just a full line
if nargin < 1, header = '';end

% default length
if (nargin < 2) || isempty(len),len = 60;end

% default separator character
if (nargin < 3) || isempty(c),c = '=';end

% get length of texgt
headerLen = length(header);

% if it is longer than the desired header length, then
% display two lines of separators one above and below the header
if (headerLen+2) >= len
  disp(repmat(c,1,len));
  disp(header)
  disp(repmat(c,1,len));
elseif headerLen == 0
  % if the header is empty, just display a full line
  disp(repmat(c,1,len));
else
  % otherwise put header inside separator characters
  fillerLen = ((len-(headerLen+2))/2);
  
  % display the first part
  disp(sprintf('%s %s %s',repmat(c,1,floor(fillerLen)),header,repmat(c,1,ceil(fillerLen))));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   displayInterpolated   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayInterpolated(sysNum,params)

global gSystem;
gSystem{sysNum}.displayInterpolated = params.displayInterpolated;
gSystem{sysNum}.vols(2).curCoord = nan;
gSystem{sysNum}.vols(2).coord = nan;
refreshDisplay(sysNum);

%%%%%%%%%%%%%%%%%%%%%%%%%
%    displayControls    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function displayControls(sysNum)

global gSystem;

paramsInfo = {...
    {'toggleOverlay',0,'type=pushbutton','callback',@toggleOverlay,'buttonString','Toggle overlay','callbackArg',sysNum,'Toggle display the overlay'}...
    {'overlayAlpha',gSystem{sysNum}.overlayAlpha,'incdec=[-0.2 0.2]','minmax=[0 1]','callback',@overlayAlpha,'callbackArg',sysNum,'passParams=1','Change the alpha of the overlay to make it more or less transparent'}...
    {'displayInterpolated',1,'type=checkbox','callback',@displayInterpolated,'callbackArg',sysNum,'Display image interpolated to match the primary volume'}...
    {'initFromHeader',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Init from header','callbackArg',{sysNum 'initFromHeader'},'passParams=1','Reinit the alignment using the qform/sform info from the headers'}...
    {'setToIdentity',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Set to identity','callbackArg',{sysNum 'setToIdentity'},'passParams=1','Set the alignment to identity'}...
    {'swapXY',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Swap XY','callbackArg',{sysNum 'swapXY'},'passParams=1','Swap XY in the alignment'}...
    {'swapXZ',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Swap XZ','callbackArg',{sysNum 'swapXZ'},'passParams=1','Swap XZ in the alignment'}...
    {'swapYZ',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Swap YZ','callbackArg',{sysNum 'swapYZ'},'passParams=1','Swap YZ in the alignment'}...
    {'flipX',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Flip X','callbackArg',{sysNum 'flipX'},'passParams=1','Flip X axis in alignment'}...
    {'flipY',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Flip Y','callbackArg',{sysNum 'flipY'},'passParams=1','Flip Y axis in alignment'}...
    {'flipZ',0,'type=pushbutton','callback',@adjustAlignment,'buttonString','Flip Z','callbackArg',{sysNum 'flipZ'},'passParams=1','Flip Z axis in alignment'}...
    {'shiftX',gSystem{sysNum}.xformParams.shiftX,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'shiftX'},'passParams=1','Shift X axis in alignment in units of voxels'}...
    {'shiftY',gSystem{sysNum}.xformParams.shiftY,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'shiftX'},'passParams=1','Shift Y axis in alignment in units of voxels'}...
    {'shiftZ',gSystem{sysNum}.xformParams.shiftZ,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'shiftX'},'passParams=1','Shift Z axis in alignment in units of voxels'}...
    {'rotateXY',gSystem{sysNum}.xformParams.rotateXY,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'rotateXY'},'passParams=1','Rotate in XY plane in units of degrees'}...
    {'rotateXZ',gSystem{sysNum}.xformParams.rotateXZ,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'rotateXZ'},'passParams=1','Rotate in XZ plane in units of degrees'}...
    {'rotateYZ',gSystem{sysNum}.xformParams.rotateYZ,'incdec=[-1 1]','callback',@adjustAlignment,'callbackArg',{sysNum 'rotateYZ'},'passParams=1','Rotate in YZ plane in units of degrees'}...
    {'vol2vol',gSystem{sysNum}.vols(2).xform,'callback',@adjustAlignment,'callbackArg',{sysNum 'vol2vol'},'passParams=1','Directly set the alignment transform - this gets composited with the shift and rotate paramters above'}
	     };


%mrParamsDialog(paramsInfo,'mlrVol Controls',[],@controlsCallback);
mrParamsDialog(paramsInfo,'mlrVol Controls');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    controlsCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function controlsCallback(params)


%%%%%%%%%%%%%%%%%%%%%%%%%
%    adjustAlignment    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = adjustAlignment(args,params)

% some variables
global gSystem;
retval = [];
sysNum = args{1};
command = args{2};
replaceXform = false;
changeXform = true;
verbose = gSystem{sysNum}.verbose;

% get system xform
gSystem{sysNum}.xformParams = params;
%gSystem{sysNum}.xform = [1 0 0 params.shiftX;0 1 0 params.shiftY;0 0 1 params.shiftZ; 0 0 0 1];
% make rotation matrix, but need to rotate around center coordinates
dim = gSystem{sysNum}.vols(2).h.dim;
shiftToCenterOfVol = makeRotMatrix3D(0,0,0,-[dim(1)/2 dim(2)/2 dim(3)/2]);
gSystem{sysNum}.xform = inv(shiftToCenterOfVol) * makeRotMatrix3D(params.rotateXZ,params.rotateYZ,params.rotateXY,[params.shiftX params.shiftY params.shiftZ],1)*shiftToCenterOfVol;


% get the necessary xform
switch args{2}
 case {'swapXY'}
  xform = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];
 case {'swapXZ'}
  xform = [0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1];
 case {'swapYZ'}
  xform = [1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
 case {'flipX'}
  % set the xform - the nan will get set to the image size below
  xform = [-1 0 0 nan;0 1 0 0;0 0 1 0; 0 0 0 1];
 case {'flipY'}
  xform = [1 0 0 0;0 -1 0 nan;0 0 1 0; 0 0 0 1];
 case {'flipZ'}
  xform = [1 0 0 0;0 1 0 0;0 0 -1 nan; 0 0 0 1];
 case {'shiftX','shiftY','shiftZ','rotateXY','rotateXZ','rotateYZ'}
  changeXform = false;
 case {'initFromHeader'}
  xform = getVol2vol(sysNum,gSystem{sysNum}.vols(2),gSystem{sysNum}.vols(1));
  replaceXform = true;
 case {'vol2vol'}
  xform = params.vol2vol;
  replaceXform = true;
 case {'setToIdentity'}
  xform = eye(4);
  replaceXform = true;
end

% set the transform
for iVol = 1:gSystem{sysNum}.n
  if verbose,dispHeader;end
  if gSystem{sysNum}.vols(iVol).tethered
    if changeXform
      % see if there is a nan that needs to be replaced
      [row col] = find(isnan(xform));
      if ~isempty(row)
	% replace from the appropriate coordinate
	val = gSystem{sysNum}.vols(iVol).h.dim(row);
	thisXform = xform;
	thisXform(row,col) = val;
      else
	thisXform = xform;
      end
      if ~replaceXform
	% display what we are doing
	if verbose,disp(sprintf('(mlrVol) Compositing xform\n%s',mrnum2str(thisXform,'compact=0','sigfigs=-1')));end
        % now set the xform
	gSystem{sysNum}.vols(iVol).xform = gSystem{sysNum}.vols(iVol).xform*thisXform;
      else
	gSystem{sysNum}.vols(iVol).xform = xform;
      end
    end
    if verbose
      % display what we are doing
      disp(sprintf('(mlrVol) xform\n%s',mrnum2str(gSystem{sysNum}.vols(iVol).xform,'compact=0','sigfigs=-1')));
      % display system xform
      disp(sprintf('(mlrVol) System xform\n%s',mrnum2str(gSystem{sysNum}.xform,'compact=0','sigfigs=-1')));
      % display complete xform
      disp(sprintf('(mlrVol) xform after compositing system\n%s',mrnum2str(shiftOriginXform*applySystemXform(sysNum,gSystem{sysNum}.vols(iVol).xform)*inv(shiftOriginXform),'compact=0','sigfigs=-1')));
    end
    % clear cache
    gSystem{sysNum}.vols(iVol).c = mrCache('init',2*max(gSystem{sysNum}.vols(iVol).h.dim(1:3)));
  end
  % set all images to redisplay
  gSystem{sysNum}.vols(iVol).curCoord = nan;
end

% redisplay
refreshDisplay(sysNum);

%%%%%%%%%%%%%%%%%%%%%%
%    overlayAlpha    %
%%%%%%%%%%%%%%%%%%%%%%
function retval = overlayAlpha(sysNum,params)

global gSystem;
gSystem{sysNum}.overlayAlpha = params.overlayAlpha;

% go through volumes
for iVol = 1:gSystem{sysNum}.n
  % for each tethered
  if gSystem{sysNum}.vols(iVol).tethered
    % go through each view and toggle alpha
    for iView = 1:3
      alphaData = get(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData');
      mask = gSystem{sysNum}.vols(iVol).overlayMask{iView};
      alphaData(mask) = gSystem{sysNum}.overlayAlpha;
      set(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData',alphaData);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    toggleOverlay    %
%%%%%%%%%%%%%%%%%%%%%%%
function retval = toggleOverlay(sysNum)

retval = [];
global gSystem;

% go through volumes
for iVol = 1:gSystem{sysNum}.n
  % for each tethered
  if gSystem{sysNum}.vols(iVol).tethered
    % go through each view and toggle alpha
    for iView = 1:3
      alphaData = get(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData');
      mask = gSystem{sysNum}.vols(iVol).overlayMask{iView};
      if alphaData(first(find(mask)))==0
	gSystem{sysNum}.overlayToggleState = 1;
	alphaData(mask) = gSystem{sysNum}.overlayAlpha;
	set(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData',alphaData);
      else
	alphaData(mask) = 0;
	gSystem{sysNum}.overlayToggleState = 0;
	set(gSystem{sysNum}.vols(iVol).overlay(iView),'AlphaData',alphaData);
      end	  
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   showSecondRowControls   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showSecondRowControls(sysNum,show)

if show
  % move all hCoordTextboxs up one row
  if any(coord<1) || any(coord(:)>vol.h.dim(1:nDim))
    for iCoord = 1:nDim
      set(gSystem{sysNum}.hCoordTextbox(1),'String',vol.coord(iCoord));
    end
  end
end
