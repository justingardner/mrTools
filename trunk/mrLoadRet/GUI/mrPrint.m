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
paramsInfo{end+1} = {'colorbarTitle',viewGet(v,'overlayName'),'Title of the colorbar'};
if flatAnat
  paramsInfo{end+1} = {'maskType',{'Circular','Remove black','None'},'Masks out anatomy image. Circular finds the largest circular aperture to view the anatomy through. Remove black keeps the patch the same shape, but removes pixels at the edge that are black.'};
end
if ~isempty(roi)
  paramsInfo{end+1} = {'roiLineWidth',1,'incdec=[-1 1]','minmax=[0 inf]','Line width for drawing ROIs. Set to 0 if you don''t want to display ROIs.'};
  paramsInfo{end+1} = {'roiColor',{'default','yellow','magenta','cyan','red','green','blue','white','black'},'Color to use for drawing ROIs. Select default to use the color currently being displayed.'};
  paramsInfo{end+1} = {'roiOutOfBoundsMethod',{'Remove','Max radius'},'If there is an ROI that extends beyond the circular aperture, you can either not draw the lines (Remove) or draw them at the edge of the circular aperture (Max radius). This is only important if you are using a circular aperture.'};
  paramsInfo{end+1} = {'roiLabels',1,'type=checkbox','Print ROI name at center coordinate of ROI'};
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
  % find the largest circular aperture the fits the image
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

% calcuate directions
params.plotDirections = 0;
if params.plotDirections
  % calculate gradient on baseCoords
  baseCoords = viewGet(v,'cursliceBaseCoords');
  baseCoords(baseCoords==0) = nan;

  fxx = baseCoords(round(end/2),:,1);
  fxx = fxx(~isnan(fxx));
  fxx = fxx(end)-fxx(1);
  fxy = baseCoords(:,round(end/2),1);
  fxy = fxy(~isnan(fxy));
  fxy = fxy(end)-fxy(1);

  fyx = baseCoords(round(end/2),:,2);
  fyx = fyx(~isnan(fyx));
  fyx = fyx(end)-fyx(1);
  fyy = baseCoords(:,round(end/2),2);
  fyy = fyy(~isnan(fyy));
  fyy = fyy(end)-fyy(1);

  fzx = baseCoords(round(end/2),:,3);
  fzx = fzx(~isnan(fzx));
  fzx = fzx(end)-fzx(1);
  fzy = baseCoords(:,round(end/2),3);
  fzy = fzy(~isnan(fzy));
  fzy = fzy(end)-fzy(1);

%samplingSize = 4;
%[fxx fxy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,1));
%[fyx fyy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,2));
%[fzx fzy] = gradient(baseCoords(1:samplingSize:end,1:samplingSize:end,3));

% get mean direction 
%fxx = mean(fxx(~isnan(fxx(:))));fxy = mean(fxy(~isnan(fxy)));
%fyx = mean(fyx(~isnan(fyx(:))));fyy = mean(fyy(~isnan(fyy)));
%fzx = mean(fzx(~isnan(fzx(:))));fzy = mean(fzy(~isnan(fzy)));

  startx = 0.15;starty = 0.8;maxlength = 0.075;
  scale = maxlength/max(abs([fxx fxy fyx fyy fzx fzy]));

  annotation('textarrow',startx+[scale*fzx 0],starty+[scale*fzy 0],'String','Left','HeadStyle','none');
  annotation('arrow',startx+[0 scale*fzx],starty+[0 scale*fzy]);
  annotation('textarrow',startx+[-scale*fxx 0],starty+[-scale*fxy 0],'String','Dorsal','HeadStyle','none');
  annotation('arrow',startx+[0 -scale*fxx],starty+[0 -scale*fxy]);
  annotation('textarrow',startx+[scale*fyx 0],starty+[scale*fyy 0],'String','Anterior','HeadStyle','none');
  annotation('arrow',startx+[0 scale*fyx],starty+[0 scale*fyy]);
end

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
label = {};
disppercent(-inf,'(mrPrint) Rendering ROIs');
for rnum = 1:length(roi)
  % check for lines
  if params.roiLineWidth > 0
    if ~isempty(roi{rnum})
      if isfield(roi{rnum},'lines')
	if ~isempty(roi{rnum}.lines.x)
	  % get color
	  if strcmp(params.roiColor,'default')
	    color = roi{rnum}.color;
	  else
	    color = params.roiColor;
	  end
	  % labels for rois, just create here
	  % and draw later so they are always on top
	  if params.roiLabels
	    x = roi{rnum}.lines.x;
	    y = roi{rnum}.lines.y;
	    label{end+1}.x = median(x(~isnan(x)));
	    label{end}.y = median(y(~isnan(y)));
	    label{end}.str = viewGet(v,'roiName',rnum);
	    label{end}.color = color;
	  end
	  % if we have a circular apertuer then we need to
	  % fix all the x and y points so they don't go off the end
	  if strcmp(params.maskType,'Circular')
	    % get the distance from center
	    x = roi{rnum}.lines.x-xCenter;
	    y = roi{rnum}.lines.y-yCenter;
	    d = sqrt(x.^2+y.^2);
	    if strcmp(params.roiOutOfBoundsMethod,'Max radius')
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
	      % set them back in the structure
	      roi{rnum}.lines.x = xsign.*abs(x)+xCenter;
	      roi{rnum}.lines.y = ysign.*abs(y)+yCenter;
	    else
	      % set all values greater than the radius to nan
	      x(d>circd) = nan;
	      y(d>circd) = nan;
	      
	      % set them back in the strucutre
	      roi{rnum}.lines.x = x+xCenter;
	      roi{rnum}.lines.y = y+yCenter;
	    end
	  end
	  % draw the lines
	  line(roi{rnum}.lines.x,roi{rnum}.lines.y,'Color',color,'LineWidth',params.roiLineWidth);
	end
      end
    end
  end
end
for i = 1:length(label)
  h = text(label{i}.x,label{i}.y,label{i}.str);
  set(h,'Color',foregroundColor);
  set(h,'Interpreter','None');
  set(h,'EdgeColor',label{i}.color);
  set(h,'BackgroundColor',params.backgroundColor);
  set(h,'FontSize',10);
end
disppercent(inf);

% bring up print dialog
global mrPrintWarning
if isempty(mrPrintWarning)
  mrWarnDlg('(mrPrint) If you are having trouble getting colors to print correctly, try exporting the figure to an eps (use its File/Export Setup menu) and printing that');
  mrPrintWarning = 1;
end
%printdlg(f);
%printpreview(f);

