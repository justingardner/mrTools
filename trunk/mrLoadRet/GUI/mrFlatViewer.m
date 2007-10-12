% mrFlatViewer.m
%
%      usage: mrFlatViewer(flatname,GM,WM,curv,anat)
%         by: modified by jg from surfViewer by eli merriam
%       date: 10/09/07
%    purpose: 
%
function retval = mrFlatViewer(flat,GM,WM,curv,anat)

% check arguments
if ~any(nargin == [1 5])
  help mrFlatViewer
  return
end

% if passed in a string check to see if
% it needs an extension
if isstr(flat)
  if isfile(sprintf('%s.off',stripext(flat)));
    flat = sprintf('%s.off',stripext(flat));
  end
end

% see how we are being called
if (nargin == 1) && isstr(flat) && ~isfile(flat)
  event = flat;
else
  event = 'init';
  % set defaults
  if ieNotDefined('GM'),GM = {};end
  if ieNotDefined('WM'),WM = {};end
  if ieNotDefined('curv'),curv = {};end
  if ieNotDefined('anat'),anat = {};end
  % make everybody a cell array
  flat = cellArray(flat);
  GM = cellArray(GM);
  WM = cellArray(WM);
  curv = cellArray(curv);
  anat = cellArray(anat);
end

switch (event)
 case 'init'
  initHandler(flat,GM,WM,curv,anat);
 case {'vSlider','hSlider'}
  sliderHandler;
 case {'edit'}
  editHandler;
 otherwise
  disp(sprintf('(mrFlatViewer) Could not find flat file %s',event));
end

%%%%%%%%%%%%%%%%%%%%%%
%%   init handler   %%
%%%%%%%%%%%%%%%%%%%%%%
function initHandler(flat,GM,WM,curv,anat)

global gFlatViewer;

disppercent(-inf,'(mrFlatView) Loading surfaces');
% load the flat
gFlatViewer.flat = loadSurfOFF(sprintf('%s.off',stripext(flat{1})));
if isempty(gFlatViewer.flat) || ~isfield(gFlatViewer.flat,'parentSurfaceName');
  disp(sprintf('(mrFlatViewer) %s is not a flat file',flat{1}));
  return
end

% remove any paths
gFlatViewer.flat.parentSurfaceName = getLastDir(gFlatViewer.flat.parentSurfaceName);

% load up the surfaces
if isempty(WM)
  WM{1} = gFlatViewer.flat.parentSurfaceName;
end
gFlatViewer.surfaces.WM = loadSurfOFF(WM{1});

% load the gray
if isempty(GM)
  GM{1} = sprintf('%sGM.off',strtok(stripext(WM{1}),'WM'));
end
gFlatViewer.surfaces.GM = loadSurfOFF(GM{1});

% load the curvature
if isempty(curv)
  curv{1} = sprintf('%s_Curv.vff',stripext(WM{1}));
end
gFlatViewer.curv = loadVFF(curv{1})';
disppercent(inf);

% load the volume
if isempty(anat)
  anat{1} = 'jg041001.hdr';
end
if isfile(anat{1})
  [gFlatViewer.anat.data gFlatViewer.anat.hdr] = cbiReadNifti(anat{1});
else
  gFlatViewer.anat = [];
end

% select the window
gFlatViewer.f = selectGraphWin;

% positions on figure
figLeft = 10;figBottom = 10;
sliderWidth = 20;sliderLength = 200;spacer = 10;
editWidth = 40;editHeight = 20;

% set up horizontal and vertical slider
gFlatViewer.hSliders.v = uicontrol('Style','slider','Position',[figLeft figBottom+sliderWidth sliderWidth sliderLength],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''vSlider'')','TooltipString','Rotate around y-axis');
gFlatViewer.hSliders.vText = uicontrol('Style','Edit','Position',[figLeft figBottom+sliderWidth+sliderLength+spacer editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0','HorizontalAlignment','Center');
gFlatViewer.hSliders.h = uicontrol('Style','slider','Position',[figLeft+sliderWidth figBottom sliderLength sliderWidth],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''hSlider'')','TooltipString','Rotate around z-axis');
gFlatViewer.hSliders.hText = uicontrol('Style','Edit','Position',[figLeft+sliderLength+3*spacer figBottom editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0');


% set they we are viewing white matter
gFlatViewer.whichSurface = 1;
gFlatViewer.patchColoring = 1;
% and display surface
dispSurface;
setViewAngle(0,0);

editable = 0;

% set up the parameters
paramsInfo = {};
if ~editable && (length(flat) == 1)
  paramsInfo{end+1} = {'flatFile',flat{1},'editable=0','The flat patch file'};
else
  paramsInfo{end+1} = {'flatFile',flat,'The flat patch file'};
end
if ~editable && (length(GM) == 1)
  paramsInfo{end+1} = {'GMSurface',GM{1},'editable=0','The gray matter file'};
else
  paramsInfo{end+1} = {'GMSurface',GM,'The gray matter file'};
end
if ~editable && (length(WM) == 1)
  paramsInfo{end+1} = {'WMSurface',WM{1},'editable=0','The white matter file'};
else
  paramsInfo{end+1} = {'WMSurface',WM,'The white matter file'};
end
if ~editable && (length(curv) == 1)
  paramsInfo{end+1} = {'curvature',curv{1},'editable=0','The curvature file'};
else
  paramsInfo{end+1} = {'curvature',curv,'The curvature file'};
end
if ~editable && (length(anat) == 1)
  paramsInfo{end+1} = {'anatomy',anat{1},'editable=0','The 3D anatomy file'};
else
  paramsInfo{end+1} = {'anatomy',anat,'The 3D anatomy file'};
end
% Now give choice of viewing gray or white
paramsInfo{end+1} = {'whichSurface',{'Gray matter','White matter','3D Anatomy','Patch'},'callback',@whichSurfaceCallback,'Choose which surface to view the patch on'};
gFlatViewer.patchColoringTypes = {'Uniform','Rostral in red','Right in red','Dorsal in red','Positive curvature in red','Negative curvature in red','Compressed areas in red','Stretched areas in red','High areal distortion in red'};
paramsInfo{end+1} = {'patchColoring',gFlatViewer.patchColoringTypes,'Choose how to color the patch','callback',@patchColoringCallback};
params = mrParamsDialog(paramsInfo,'View flat patch location on surface');

if isempty(params)
  close(gFlatViewer.f);
else
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   patchColoringCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function patchColoringCallback(params)

global gFlatViewer
gFlatViewer.patchColoring = find(strcmp(params.patchColoring,gFlatViewer.patchColoringTypes));
if gFlatViewer.whichSurface == 3
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  dispVolume(3,hPos);
else
  dispSurface;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   whichSurfaceCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichSurfaceCallback(params)

global gFlatViewer;

% get which surface to draw
lastWhichSurface = gFlatViewer.whichSurface;
whichSurface = find(strcmp(params.whichSurface,{'Gray matter','White matter','3D Anatomy','Patch'}));
if whichSurface ~= lastWhichSurface
  % set which surface and display
  gFlatViewer.whichSurface = whichSurface;
  % 1,2 are surfaces
  if whichSurface <= 2
    % if we are displaying the 3D anatomy, 
    % then switch to the surface view
    if lastWhichSurface > 2
      switchToSurface;
    else
      dispSurface;
    end
  % 3 is the volume
  elseif whichSurface == 3
    % switch to the volume view
    switchToVolume;
  else
    % switch to the volume view
    switchToFlat;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function switchToSurface

global gFlatViewer;
set(gFlatViewer.hSliders.v,'Visible','on');
set(gFlatViewer.hSliders.vText,'Visible','on');
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.h,'SliderStep',[15 45]./360);
set(gFlatViewer.hSliders.h,'Value',0);
set(gFlatViewer.hSliders.h,'Min',-180);
set(gFlatViewer.hSliders.h,'Max',180);
set(gFlatViewer.hSliders.h,'TooltipString','Rotate around z-axis');
set(gFlatViewer.hSliders.v,'Value',0);
set(gFlatViewer.hSliders.vText,'String',0);
set(gFlatViewer.hSliders.hText,'String',0);
dispSurface;
setViewAngle(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToVolume   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function switchToVolume

global gFlatViewer;
initSlice = 127;
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.v,'Visible','off');
set(gFlatViewer.hSliders.vText,'Visible','off');
set(gFlatViewer.hSliders.h,'SliderStep',[1 16]./256);
set(gFlatViewer.hSliders.h,'Value',initSlice);
set(gFlatViewer.hSliders.h,'Min',1);
set(gFlatViewer.hSliders.h,'Max',256);
set(gFlatViewer.hSliders.h,'TooltipString','Change viewing slice');
set(gFlatViewer.hSliders.hText,'String',num2str(initSlice));
dispVolume(3,initSlice);

%%%%%%%%%%%%%%%%%%%%%%
%%   switchToFlat   %%
%%%%%%%%%%%%%%%%%%%%%%
function switchToFlat

global gFlatViewer;
initSlice = 127;
set(gFlatViewer.hSliders.v,'Visible','off');
set(gFlatViewer.hSliders.vText,'Visible','off');
set(gFlatViewer.hSliders.h,'Visible','off');
set(gFlatViewer.hSliders.hText,'Visible','off');
dispSurface;

%%%%%%%%%%%%%%%%%%%%%%%
%%   sliderHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function sliderHandler

global gFlatViewer;

% get slider position
hPos = round(get(gFlatViewer.hSliders.h,'Value'));
vPos = round(get(gFlatViewer.hSliders.v,'Value'));

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if gFlatViewer.whichSurface <= 2
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   editHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function editHandler

global gFlatViewer;

hPos = str2num(get(gFlatViewer.hSliders.hText,'String'));
vPos = str2num(get(gFlatViewer.hSliders.vText,'String'));

% make it fit into -180:180
hPos = round(mod(hPos+180,360)-180);
vPos = round(mod(vPos+180,360)-180);

% set slider position
set(gFlatViewer.hSliders.h,'Value',hPos);
set(gFlatViewer.hSliders.v,'Value',vPos);

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if gFlatViewer.whichSurface <= 2
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   setViewAngle   %%
%%%%%%%%%%%%%%%%%%%%%%
function setViewAngle(hPos,vPos)

% flip the sign to make rotations go in the "right" direction
hPos = -hPos;vPos = -vPos;
% somehow 90 and 180 are a problem for matlab
if abs(vPos) == 90,vPos = sign(vPos)*91;,end
if abs(hPos) == 90,hPos = sign(hPos)*91;,end
if abs(hPos) == 179,hPos = sign(hPos)*179;,end

% set the camera taret to center
camtarget([0 0 0]);

% set the view angle
view(hPos,vPos);

% change the camera position to avoid the volume
% flipping back and forth, another starnge matlab thing
if (vPos >= 90) || (vPos < -90)
  camup([0 0 -1]);
else
  camup([0 0 1]);
end

%%%%%%%%%%%%%%%%%%%%%
%%   dispSurface   %%
%%%%%%%%%%%%%%%%%%%%%
function dispSurface

global gFlatViewer;
figure(gFlatViewer.f);

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

% get the vertexes/triangles and curvature
if gFlatViewer.whichSurface == 1
  vtcs = gFlatViewer.surfaces.WM.vtcs;
  tris = gFlatViewer.surfaces.WM.tris;
  c = gFlatViewer.curv;
elseif gFlatViewer.whichSurface == 2
  vtcs = gFlatViewer.surfaces.GM.vtcs;
  tris = gFlatViewer.surfaces.GM.tris;
  c = gFlatViewer.curv;
else
  vtcs = gFlatViewer.flat.vtcs;
  tris = gFlatViewer.flat.tris;
  c = gFlatViewer.curv(patchVtcs);
%  c = (c-min(c))./((max(c)-min(c)))>0.5;
  view([0 90]);
end  

% move vertices into center
vtcs(:,1) = vtcs(:,1)-mean(vtcs(:,1));
vtcs(:,2) = vtcs(:,2)-mean(vtcs(:,2));
vtcs(:,3) = vtcs(:,3)-mean(vtcs(:,3));

% clear the axis
cla;

% not sure why, but this is necessary to set up
% the axis so that right is right...
imagesc(0);

% get the colors that we want to show for that patch
co = getPatchColoring;

% now set the overlay
if gFlatViewer.whichSurface <= 2
  overlay = NaN(length(c),3);
  overlay(patchVtcs,1) = co(:,1);
  overlay(patchVtcs,2:3) = co(:,2:3);
else
  overlay(:,1) = co(:,1);
  overlay(:,2:3) = co(:,2:3);
end

% draw the surface and the overlay
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', c, ...
      'facecolor', 'interp', ...
      'edgecolor', 'none');
patch('vertices', vtcs, 'faces', tris, ...
      'FaceVertexCData', overlay, 'FaceVertexAlphaData', overlay(:,1)*.1, ...
      'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',.6);

% set axis stuff
axis off;axis equal;colormap(gray);axis tight;
camup('manual');
set(gca,'CLim',[-1.2 1.2]);

if gFlatViewer.whichSurface <= 2
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  vPos = round(get(gFlatViewer.hSliders.v,'Value'));
  setViewAngle(hPos,vPos);
end
%%%%%%%%%%%%%%
% dispVolume
%%%%%%%%%%%%%%
function dispVolume(sliceIndex,slice)

global gFlatViewer;
figure(gFlatViewer.f);
cla;

% display a slice of the anatomy image
switch sliceIndex
  case {1}
   img = gFlatViewer.anat.data(slice,:);
  case {2}
   img = gFlatViewer.anat.data(:,slice,:);
  case {3}
   img = gFlatViewer.anat.data(:,:,slice);
end
imagesc(img);
colormap(gray);
axis image;
axis off;
hold on


set(gca,'CLim',[min(img(:)) max(img(:))]);
% display patch and white matter/gray matter
whichInx = gFlatViewer.flat.patch2parent(:,2);
wmPatchNodes = gFlatViewer.surfaces.WM.vtcs(whichInx,:);
gmPatchNodes = gFlatViewer.surfaces.GM.vtcs(whichInx,:);

% get full white matter/gray matter nodes
wmNodes = gFlatViewer.surfaces.WM.vtcs;
gmNodes = gFlatViewer.surfaces.GM.vtcs;

% Plot the nodes for the gray/white matter surfaces
wmNodes = wmNodes( find( round(wmNodes(:,sliceIndex))==slice), : );
plot(wmNodes(:,1), wmNodes(:,2), 'w.', 'markersize', 1);

gmNodes = gmNodes( find( round(gmNodes(:,sliceIndex))==slice), : );
plot(gmNodes(:,1), gmNodes(:,2), 'y.', 'markersize', 1);


% plot the patch nodes, displaying both deep and superficial surfaces
co = getPatchColoring;
wmco = co(find( round(wmPatchNodes(:,sliceIndex))==slice),:);
% make into magenta vs blue
i = wmco(:,1);
wmco(:,1) = i;
wmco(:,2) = 0;
wmco(:,3) = max(i,1-i);

wmPatchNodes = wmPatchNodes( find( round(wmPatchNodes(:,sliceIndex))==slice), : );

if gFlatViewer.patchColoring == 1
  plot(wmPatchNodes(:,1), wmPatchNodes(:,2), '.', 'markersize', 1,'Color',[1 0 1]);
% otherwise each pixel has to be set
else
  for i = 1:length(wmPatchNodes(:,1))
    plot(wmPatchNodes(i,1), wmPatchNodes(i,2), '.', 'markersize', 1,'Color',wmco(i,:)');
  end
end

gmco = co(find( round(gmPatchNodes(:,sliceIndex))==slice),:);
gmPatchNodes = gmPatchNodes( find( round(gmPatchNodes(:,sliceIndex))==slice), : );
% uniform patch coloring
if gFlatViewer.patchColoring == 1
  plot(gmPatchNodes(:,1), gmPatchNodes(:,2), '.', 'markersize', 1,'Color',co(1,:));
% otherwise each pixel has to be set
else
  for i = 1:length(gmPatchNodes(:,1))
    plot(gmPatchNodes(i,1), gmPatchNodes(i,2), '.', 'markersize', 1,'Color',gmco(i,:));
  end
end

view([0 90]);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getPatchColoring   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function co = getPatchColoring

global gFlatViewer;

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

title('');
% one is uniform, 2-4 are red/blue
switch gFlatViewer.patchColoring
 % uniform
 case 1
  % make everybody red
  co = ones(1,gFlatViewer.flat.Nvtcs);
  % anatomical directions
 case {2,3,4}
  co = gFlatViewer.surfaces.GM.vtcs(patchVtcs,gFlatViewer.patchColoring-1)';
  co = (co-min(co))./(max(co)-min(co));
  % curvature
 case {5,6}
  % get curvature
  curv = gFlatViewer.curv(patchVtcs)';
  curv = (curv-min(curv))./(max(curv)-min(curv));
  if gFlatViewer.patchColoring == 6
    curv = 1-curv;
  end
  % flatten distribution
  co = flattenDistribution(curv);
 % Areal distortion
 case {7,8,9}
  tris = gFlatViewer.flat.tris;
  % get area in volume
  volumeArea = getTriangleArea(gFlatViewer.surfaces.GM.vtcs(patchVtcs,:),tris);
  patchArea = getTriangleArea(gFlatViewer.flat.vtcs,tris);
  % get the distortion as the ratio of the area in the
  % patch to the volume
  trisDistortion = patchArea./volumeArea;
  % now convert this into a color for each vertex
  distortion = ones(1,length(patchVtcs));
  % note that this is a shortcut, it just
  % sets each vertex to one value of the distortion
  % even though the vertex may belong to many triangles
  distortion(tris(:,1)) = trisDistortion;
  distortion(tris(:,2)) = trisDistortion;
  distortion(tris(:,3)) = trisDistortion;
  % normalize
  if gFlatViewer.patchColoring < 9
    co = (distortion-min(distortion))./(max(distortion)-min(distortion));
    co = flattenDistribution(co);
    % invert colors
    distortion = log10(distortion);
    if gFlatViewer.patchColoring == 8
      co = 1-co;
      title(sprintf('Stretch: max=%0.2fx median=%0.2fx mean=%0.2fx',10^max(distortion),10^median(distortion(distortion>0)),10^mean(distortion(distortion>0))));
    else
      title(sprintf('Compression: max=%0.2fx median=%0.2fx mean=%0.2fx',10^abs(min(distortion)),10^abs(median(distortion(distortion<0))),10^abs(mean(distortion(distortion<0)))));
    end
  else
    % set the colors to the absolute value
    % of the log of the distortion. This scales
    % so that doubling or halving the area from
    % the volume to the patch give you the same
    % number. Anything above 10x distortion is 1
    distortion = abs(log10(distortion));
    co = distortion;
    co(co>1) = 1;
    % print out some statistics
    distortion = 10.^distortion;
    title(sprintf('Distortion: max=%0.2fx median=%0.2fx mean=%0.2fx',max(distortion),median(distortion),mean(distortion)));
  end
end

% intermediate values turn to gray, so avoid them
co((co>0.3)&(co<0.5)) = 0.3;
co((co>0.5)&(co<0.7)) = 0.7;

% make into RGB
co(2:3,:) = [1-co;1-co];
co = co';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getTriangleArea   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function area = getTriangleArea(vtcs,tris)

% get length of a each side of triangles
a = sqrt(sum((vtcs(tris(:,1),:)-vtcs(tris(:,2),:))'.^2));
b = sqrt(sum((vtcs(tris(:,2),:)-vtcs(tris(:,3),:))'.^2));
c = sqrt(sum((vtcs(tris(:,3),:)-vtcs(tris(:,1),:))'.^2));

% get semiperimeter (i.e. 1/2 perimeter)
p = (a+b+c)/2;

% use Heron's formula for size of triangle given side lengths
area = sqrt(p.*(p-a).*(p-b).*(p-c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   flattenDistribution   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function co = flattenDistribution(co,nbins)

if ~exist('nbins','var'),nbins = 10;end

% sort the values
[co sortIndex] = sort(co);
% get binsize
binSize = floor(length(co)/nbins);
% and even number of values back into each bin
for i = 1:nbins
  co(sortIndex((i-1)*binSize+1:min(length(co),i*binSize))) = i/nbins;
end
