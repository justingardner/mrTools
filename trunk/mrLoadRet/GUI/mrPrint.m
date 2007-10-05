% mrPrint.m
%
%        $Id$
%      usage: mrPrint(v)
%         by: justin gardner
%       date: 10/04/07
%    purpose: puts a printable version of the data into the graph win
%
function retval = mrPrint(v)

% check arguments
if ~any(nargin == [1])
  help mrPrint
  return
end

mrGlobals;

% see if this a flat or not
if isempty(viewGet(v,'baseCoordMap'))
  flatAnat = 0;
else
  flatAnat = 1;
end

% grab the image
disppercent(-inf,'(mrPrint) Rerendering image');
[img base roi] = refreshMLRDisplay(viewGet(v,'viewNum'));
disppercent(inf);

% first get parameters that the user wants to display
paramsInfo = {};
paramsInfo{end+1} = {'title',sprintf('%s: %s',getLastDir(MLR.homeDir),viewGet(v,'description')),'Title of figure'};
paramsInfo{end+1} = {'backgroundColor',{'white','black'},'Background color, either white or black'};
paramsInfo{end+1} = {'colorbarLoc',{'SouthOutside','NorthOutside','EastOutside','WestOutside','None'},'Location of colorbar, select None if you do not want a colorbar'};
if flatAnat
  paramsInfo{end+1} = {'maskType',{'Circular','Remove black','None'},'Masks out anatomy image. Circular finds the largest circular aperture to view the anatomy through. Remove black removes all pixels that are black as defined by blackValue'};
end
if ~isempty(roi)
  paramsInfo{end+1} = {'roiLineWidth',1,'incdec=[-1 1]','minmax=[0 inf]','Line width for drawing ROIs. Set to 0 if you don''t want to display ROIs. Note that ROIs will only draw if you are drawing ROI perimeters.'};
  paramsInfo{end+1} = {'roiColor',{'default','yellow','magenta','cyan','red','green','blue','white','black'},'Color to use for drawing ROIs. Select default to use the color currently being displayed.'};
paramsInfo{end+1} = {'colorbarTitle',viewGet(v,'overlayName'),'Title of the colorbar'};
end

params = mrParamsDialog(paramsInfo,'Print figure options');;
if isempty(params),return,end

% get the gui, so that we can extract colorbar
fig = viewGet(v,'figNum');
gui = guidata(fig);

% grab the colorbar data
H = get(gui.colorbar,'children');
cmap = squeeze(get(H(end),'CData'));

% display in graph window
f = selectGraphWin;
keyboard
clf(f);drawnow;

set(f,'Name','Print figure');
set(f,'NumberTitle','off');
set(f,'color',params.backgroundColor)

% value to consider to be "black" in image
blackValue = 0;

if ~flatAnat,params.maskType = 'None';end

% get the mask
if strcmp(params.maskType,'None')
  mask = ones(size(img));
elseif strcmp(params.maskType,'Remove black')
  mask(:,:,1) = (base.im<=blackValue);
  mask(:,:,2) = (base.im<=blackValue);
  mask(:,:,3) = (base.im<=blackValue);
elseif strcmp(params.maskType,'Circular')
  % find the largest circular aperture
  xCenter = (size(base.im,2)/2);
  yCenter = (size(base.im,1)/2);
  x = (1:size(base.im,2))-xCenter;
  y = (1:size(base.im,1))-yCenter;
  [x y] = meshgrid(x,y);
  % now compute the distance from the center for
  % every point
  d = sqrt(x.^2+y.^2);
  % now go an find the largest distance for which there
  % are no black voxels
  maxd = max(d(:));
  mind = min(d(:));
  for circd = maxd:-1:mind
    if ~any(base.im(d<=circd)<=blackValue)
      break;
    end
  end
  mask(:,:,1) = (d>circd);
  mask(:,:,2) = (d>circd);
  mask(:,:,3) = (d>circd);
end


% get foregroundColor
if strcmp(params.backgroundColor,'white')
  foregroundColor = [0 0 0];
else
  foregroundColor = [1 1 1];
end
  
% mask out the image
if strcmp(params.backgroundColor,'white')
  img(mask) = 1;
else
  img(mask) = 0;
end
  
% display the image
colormap(cmap);
image(img);
axis equal; axis off;axis tight;hold on

% display the colormap
if ~strcmp(params.colorbarLoc,'None')
  H = colorbar(params.colorbarLoc);
  % set the colorbar ticks, making sure to switch
  % them if we have a vertical as opposed to horizontal bar
  if ismember(params.colorbarLoc,{'EastOutside','WestOutside'})
    set(H,'XTick',get(gui.colorbar,'YTick'));
    set(H,'Ytick',get(gui.colorbar,'XTick'));
    set(H,'YTickLabel',get(gui.colorbar,'XTicklabel'));
  else
    set(H,'YTick',get(gui.colorbar,'YTick'));
    set(H,'Xtick',get(gui.colorbar,'XTick'));
    set(H,'XTickLabel',get(gui.colorbar,'XTicklabel'));
  end    
  set(H,'XColor',foregroundColor);
  set(H,'YColor',foregroundColor);
  set(get(H,'Title'),'String',params.colorbarTitle);
  set(get(H,'Title'),'Color',foregroundColor);
  set(get(H,'Title'),'FontSize',14);
end

% create a title
H = title(params.title);
set(H,'Color',foregroundColor);
set(H,'FontSize',16);

drawnow;

% draw the roi
sliceNum = viewGet(v,'currentSlice');
disppercent(-inf,'(mrPrint) Rendering ROIs');
for rnum = 1:length(roi)
  % check for lines
  if params.roiLineWidth > 0
    if ~isempty(roi{rnum})
      if isfield(roi{rnum},'lines')
	if ~isempty(roi{rnum}.lines.x)
	  % if we have a circular apertuer then we need to
	  % fix all the x and y points so they don't go off the end
	  if strcmp(params.maskType,'Circular')
	    % get the distance
	    x = roi{rnum}.lines.x-xCenter;
	    y = roi{rnum}.lines.y-yCenter;
	    d = sqrt(x.^2+y.^2);
	    % find the angle of all points
	    ang = atan(y./x);
	    ysign = (y > 0)*2-1;
	    xsign = (x > 0)*2-1;
	    newx = circd*cos(ang);
	    newy = circd*sin(ang);
	    % now reset all points past the maximum radius
	    % with values at the outermost edge of the aperture
	    x(d>circd) = newx(d>circd);
	    y(d>circd) = newy(d>circd);
	    % set them back in the strucutre
	    roi{rnum}.lines.x = xsign.*abs(x)+xCenter;
	    roi{rnum}.lines.y = ysign.*abs(y)+yCenter;
	  end
	  if strcmp(params.roiColor,'default')
	    color = roi{rnum}.color;
	  else
	    color = params.roiColor;
	  end
	  % draw the lines
	  line(roi{rnum}.lines.x,roi{rnum}.lines.y,'Color',color,'LineWidth',params.roiLineWidth);
	end
      end
    end
  end
end
disppercent(inf);

% bring up print dialog
global mrPrintWarning
if isempty(mrPrintWarning)
  mrWarnDlg('(mrPrint) If you are having trouble getting colors to print correctly, try exporting the figure to an eps (use its File/Export Setup menu) and printing that');
  mrPrintWarning = 1;
end
printdlg(f);
%printpreview(f);

