% renderRois3D - displays 3D rendering of currently loaded ROIs
%
%        $Id$
%      usage: [  ] = renderRois3D(thisView,overlayNum,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2010-02-15
%     inputs: 
%    outputs: 
%
%    purpose: displays 3D rendering of currently loaded ROIs
%
% planned improvements: take a subvolume (slightly larger than ROIs) from the start to reduce computing time, especially
%                       when resampling to base anatomy.

function renderRois3D(thisView,overlayNum,scanNum,x,y,z,roi)

computeContours = 0;

if isempty(thisView.ROIs)
   mrWarnDlg('(renderRois3D) Please load at least one ROI');
   return;
end

fprintf(1,'Loading data...')
tic
baseNum = viewGet(thisView,'currentBase');
basevoxelsize = viewGet(thisView,'basevoxelsize');
alphaOverlay = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay'));
d.alpha = viewGet(thisView,'alpha');
base_size = viewGet(thisView,'basedims',baseNum);

% get the analysis structure
analysis = viewGet(thisView,'analysis');
fprintf(1,'Done\n')
toc


%--------------------------------Mask and Alpha overlay data -------------------------------
if ~isempty(analysis)
   
   alphaOverlayNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay'));
   
   %First, get a mask of non-zero voxel representing the current overlay display
   %This is taken from computeOverlay.m
   fprintf(1,'Computing overlay mask...')
   tic
   [mask, overlayData]= maskOverlay(thisView,[overlayNum alphaOverlayNum],scanNum);
   if iscell(mask)
      mask = mask{1};
   end
   fprintf(1,'Done\n')
   toc

   baseMaskData = getBaseSpaceOverlay(thisView, double(mask), scanNum, baseNum,'nearest');
   baseMaskData(isnan(baseMaskData))=0;
   baseMaskData = logical(baseMaskData);
   
   %make alpha map
   if isempty(alphaOverlayNum)
      baseAlphaData = d.alpha*ones(size(baseMaskData));
   else
      fprintf(1,'Computing alpha overlay...')
      tic
      alphaOverlayData = overlayData{2};
      overlayData = overlayData{1};
      
      baseAlphaData = getBaseSpaceOverlay(thisView, alphaOverlayData, scanNum, baseNum);
      % get the range of the alpha overlay
      range = viewGet(thisView,'overlayRange',alphaOverlay);
      % handle setRangeToMax
      if strcmp(viewGet(thisView,'overlayCtype',alphaOverlay),'setRangeToMax')
       maxRange = max(clip(1),min(baseAlphaData(mask)));
       if ~isempty(maxRange),range(1) = maxRange;end
       minRange = min(max(baseAlphaData(mask)),clip(2));
       if ~isempty(minRange),range(2) = minRange;end
      end
      % now compute the alphaOverlay as a number from
      % 0 to 1 of the range
      baseAlphaData = d.alpha*((baseAlphaData-range(1))./diff(range));
      baseAlphaData(baseAlphaData>d.alpha) = d.alpha;
      baseAlphaData(baseAlphaData<0) = 0;
      alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
      if alphaOverlayExponent<0   %JB if the alpha overlay exponent is negative, set it positive and take the complementary values to 1 for the alpha map
        alphaOverlayExponent = -alphaOverlayExponent;
        baseAlphaData = 1-baseAlphaData;
      end
      baseAlphaData = baseAlphaData.^alphaOverlayExponent;
      baseAlphaData(isnan(baseAlphaData)) = 0;
   end   
   fprintf(1,'Done\n')
   toc

else
   baseMaskData = [];
   baseAlphaData = [];
end



%------------------------------------ Construct ROI objcts ------------------------------------------------%
fprintf(1,'Constructing 3D ROIs...')
tic
edgeX = [-.5 -.5; -.5 -.5; -.5 -.5; -.5 -.5; -.5 .5; -.5 .5; -.5 .5; -.5 .5; .5 .5; .5 .5; .5 .5; .5 .5];
edgeY = [-.5 .5; -.5 .5; -.5 -.5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5];
edgeZ = [-.5 -.5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; -.5 .5; -.5 .5];


allRoisBaseCoords = [];
for i_roi = 1:length(thisView.ROIs)

   xform = viewGet(thisView, 'base2roi',i_roi);
   roivoxelsize = viewGet(thisView,'roivoxelsize',i_roi);
   roiBaseCoords = xformROIcoords(thisView.ROIs(i_roi).coords,inv(xform),roivoxelsize,basevoxelsize);
   
   if ~isempty(roiBaseCoords)
      % we need to remove any coordinate that might fall outside the base anatomy
      outside_voxels = find(roiBaseCoords(1,:)<1 | roiBaseCoords(1,:)>base_size(1) |...
                           roiBaseCoords(2,:)<1 | roiBaseCoords(2,:)>base_size(2) |...
                           roiBaseCoords(3,:)<1 | roiBaseCoords(3,:)>base_size(3) );
      roiBaseCoords(:,outside_voxels) = [];
      d.roiSize(i_roi) = size(roiBaseCoords,2);
      %remember which voxels are in which rois
      d.roiDataIndex{i_roi} = size(allRoisBaseCoords,2)+(1:d.roiSize(i_roi));
      %put Roi coords together to simplify overlay computation
      allRoisBaseCoords = [allRoisBaseCoords roiBaseCoords];
      
      
      %Edges (grid and perimeter)
      edgeXcoords = [];
      edgeYcoords = [];
      edgeZcoords = [];

      for i_voxel = 1:size(roiBaseCoords,2)
         edgeXcoords = [edgeXcoords; roiBaseCoords(1,i_voxel)+edgeX];
         edgeYcoords = [edgeYcoords; roiBaseCoords(2,i_voxel)+edgeY];
         edgeZcoords = [edgeZcoords; roiBaseCoords(3,i_voxel)+edgeZ];
      end
      edgesCoords = [edgeXcoords edgeYcoords edgeZcoords];
      %find unique edges around  this roi
      [gridCoords, unique_indices] = unique(edgesCoords,'rows');
      double_indices = setdiff(1:size(edgesCoords,1),unique_indices);

      unique_edgesCoords = setdiff(edgesCoords,edgesCoords(double_indices,:),'rows');
      
      roiD.GridXcoords{i_roi} = gridCoords(:,1:2)';
      roiD.GridYcoords{i_roi} = gridCoords(:,3:4)';
      roiD.GridZcoords{i_roi} = gridCoords(:,5:6)';
      roiD.EdgeXcoords{i_roi} = unique_edgesCoords(:,1:2)';
      roiD.EdgeYcoords{i_roi} = unique_edgesCoords(:,3:4)';
      roiD.EdgeZcoords{i_roi} = unique_edgesCoords(:,5:6)';
      
      %Faces (opaque ROi)
      %coordinates of a cube around each voxel
      [roiD.FacesXcoords{i_roi}, roiD.FacesYcoords{i_roi}, roiD.FacesZcoords{i_roi}] = makeCubeFaces(roiBaseCoords,[],[],1);

   else
      roiD.GridXcoords{i_roi} = [];
      roiD.GridYcoords{i_roi} = [];
      roiD.GridZcoords{i_roi} = [];
      roiD.EdgeXcoords{i_roi} = [];
      roiD.EdgeYcoords{i_roi} = [];
      roiD.EdgeZcoords{i_roi} = [];
      roiD.FacesXcoords{i_roi} = [];
      roiD.FacesYcoords{i_roi} = [];
      roiD.FacesZcoords{i_roi} = [];
   end
end
fprintf(1,'Done\n')
toc


%------------------------------------ Construct 3D Overlay  ------------------------------------------------%
if ~isempty(analysis)
   baseOverlayData = getBaseSpaceOverlay(thisView, overlayData, scanNum, baseNum);

   %find indices of unique voxels in all rois
   [allRoisBaseCoords, dump, duplicateVoxels] = unique(allRoisBaseCoords','rows');
   allRoisBaseCoords = allRoisBaseCoords';
   for i_roi = 1:length(thisView.ROIs)
      d.roiDataIndex{i_roi} = duplicateVoxels(d.roiDataIndex{i_roi});
   end
   
%%%%%%%%%%%%%%%%%%% Compute cube faces
   fprintf(1,'Constructing overlay data cubes...')
   tic
   d.cubesCData = baseOverlayData;
   %mask the data
   d.cubesCData(~baseMaskData)=NaN;
   %keep only voxels in Rois
   cubesOverlayDataIndex = sub2ind(size(d.cubesCData), allRoisBaseCoords(1,:)', allRoisBaseCoords(2,:)', allRoisBaseCoords(3,:)' );
   d.cubesCData = d.cubesCData(cubesOverlayDataIndex);
   d.cubesAlphaData = baseAlphaData(cubesOverlayDataIndex);
   
   %exclude Nan overlay values
   voxelsToKeepIndices = find(~isnan(d.cubesCData));
   d.cubesCData = d.cubesCData(voxelsToKeepIndices);
   
   %sort faces by overlay data value, that will make things faster later
   [d.cubesCData, sortIndex] = sort(d.cubesCData);
   
   %apply sorting and excluding to other data arrays
   voxelsToKeepIndices = voxelsToKeepIndices(sortIndex);
   d.cubesAlphaData = d.cubesAlphaData(voxelsToKeepIndices);
   d.cubesBaseCoords = allRoisBaseCoords(:,voxelsToKeepIndices);
   for i_roi = 1:length(thisView.ROIs)
      [dump, d.roiDataIndex{i_roi}] = intersect(voxelsToKeepIndices,d.roiDataIndex{i_roi});
   end

   fprintf(1,'Done\n')
   toc
      
      
      
%%%%%%%%%%%%%%%%%%% Compute upsampled smooth volumes for blobs
   if computeContours
      fprintf(1,'Smoothing and upsampling overlay data...')
      %reduce volume
      minX = min(allRoisBaseCoords(1,:));
      maxX = max(allRoisBaseCoords(1,:));
      minY = min(allRoisBaseCoords(2,:));
      maxY = max(allRoisBaseCoords(2,:));
      minZ = min(allRoisBaseCoords(3,:));
      maxZ = max(allRoisBaseCoords(3,:));

      %volume slightly larger for smoothing
      margin = 3;
      minXmargin = max(minX-margin,1);
      maxXmargin = min(maxX+margin,size(baseOverlayData,1)); 
      minYmargin = max(minY-margin,1);                         
      maxYmargin = min(maxY+margin,size(baseOverlayData,2));   
      minZmargin = max(minZ-margin,1);
      maxZmargin = min(maxZ+margin,size(baseOverlayData,3));

      %reduce volume
      [Ycoords,Xcoords,Zcoords,d.subVolOverlayData] = subvolume(baseOverlayData,[minYmargin maxYmargin minXmargin maxXmargin minZmargin maxZmargin]); %subvolume swaps X and Y axes...
      [Ycoords,Xcoords,Zcoords,subMaskData] = subvolume(baseMaskData,[minYmargin maxYmargin minXmargin maxXmargin minZmargin maxZmargin]); %subvolume swaps X and Y axes...
      %upsample
      upsamplingFactor = 5;
      %input coordinates (watch out, X and Y must be inverted to use some volume function)
      %[y,x,z] = meshgrid(minYmargin:maxYmargin,minXmargin:maxXmargin,minZmargin:maxZmargin);
      %add one around the min and max to interpolate around the actual data (integer coordinates)
      [yi,xi,zi] = meshgrid(  (minYmargin:1/upsamplingFactor:maxYmargin),...
                              (minXmargin:1/upsamplingFactor:maxXmargin),...
                              (minZmargin:1/upsamplingFactor:maxZmargin));
      d.subVolOverlayData = interp3(Ycoords,Xcoords,Zcoords,d.subVolOverlayData,yi,xi,zi,'linear');
      subMaskData = interp3(Ycoords,Xcoords,Zcoords,subMaskData+0,yi,xi,zi,'nearest');
      d.subVolOverlayData = smooth3(d.subVolOverlayData);

      %mask the data with -Infs
      d.subVolOverlayData(~subMaskData)= NaN; 

      %keep only ROI voxels
      %get voxel coordinates in upsampled volume
      xform = diag([[1 1 1]*upsamplingFactor 1]);
      baseVoxelSize = viewGet(thisView,'basevoxelsize',baseNum);
      blobsRoisBaseCoords = xformROIcoords(allRoisBaseCoords,xform,baseVoxelSize,baseVoxelSize/upsamplingFactor);
      blobsRoisBaseCoords = blobsRoisBaseCoords(1:3,:) - repmat(upsamplingFactor*[minXmargin;minYmargin;minZmargin],1,size(blobsRoisBaseCoords,2))+1;
      blobsRoisBaseCoordsIndex = sub2ind(size(d.subVolOverlayData),blobsRoisBaseCoords(1,:),blobsRoisBaseCoords(2,:),blobsRoisBaseCoords(3,:));
      d.subVolOverlayDataRois = NaN(size(d.subVolOverlayData));
      d.subVolOverlayDataRois(blobsRoisBaseCoordsIndex) = d.subVolOverlayData(blobsRoisBaseCoordsIndex);
      [d.blobsYcoords,d.blobsXcoords,d.blobsZcoords,d.subVolOverlayDataRois] = subvolume(yi, xi, zi, d.subVolOverlayDataRois,[minY maxY minX maxX minZ maxZ]); %subvolume swaps X and Y axes...
      [d.blobsYcoords,d.blobsXcoords,d.blobsZcoords,d.subVolOverlayData] = subvolume(yi, xi, zi, d.subVolOverlayData,[minY maxY minX maxX minZ maxZ]); %subvolume swaps X and Y axes...

      fprintf(1,'Done\n')
   end
end
      
      

%---------------------- construct figure -----------------%
% selectGraphWin;
% global MLR;
% fignum = MLR.graphFigure;
fignum = figure;
% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','RenderRois3D');
set(fignum,'Toolbar','Figure');
roisButtonsBottom = .3;
overlayButtonsBottom = .1;
roisCheckHeight = min(.05,(.9 - roisButtonsBottom)/length(thisView.ROIs));
overlaySlidersLeft=.02;
rotateViewLeft = .6;

for i_roi = 1:length(thisView.ROIs)
   d.hRoiVoxelNum(i_roi) = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',...
     [.04 roisButtonsBottom-.01+roisCheckHeight*(length(thisView.ROIs)-i_roi+1) .08 roisCheckHeight/2]);
end

%Overlay control buttons
if ~isempty(analysis)
   set(fignum,'Name',['RenderRois3D - ' viewGet(thisView,'homeDir') ' - ' analysis.name ' - ' analysis.overlays(overlayNum).name ]);
   d.colorScale = analysis.overlays(overlayNum).range;
   d.colorMap = analysis.overlays(overlayNum).colormap;
   

   if computeContours
      d.hOverlayRestrictCheck = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','checkbox',...
        'Position',[.02 overlayButtonsBottom+.16 .15 .03],...
         'Min', 0,'Max',1, 'Value',1, 'String', 'Restrict Contour to ROIs');
      overlayList = {'Cubes','IsoContour'};
   else
      overlayList = {'Cubes'};
   end
   d.hOverlayList = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.02 overlayButtonsBottom .1 .03],...
      'Value',1, 'String', overlayList);
   d.hOverlayAlphaCheck = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','checkbox', 'Position',[.02 overlayButtonsBottom+.08 .1 .03],...
      'Min', 0,'Max',1, 'Value',0, 'String','Use Alpha');
   
   %Overlay sliders
   uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[overlaySlidersLeft .04 .13 .03],'String','Min Threshold');
   uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[overlaySlidersLeft .01 .13 .03],'String','Max Threshold');
   d.hOverlayEditMax = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[overlaySlidersLeft+.4 .01 .1 .03],...
      'Min', d.colorScale(1),'Max',d.colorScale(2), 'String',num2str(d.colorScale(2)));
   d.hOverlayEditMin = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[overlaySlidersLeft+.4 .04 .1 .03],...
      'Min', d.colorScale(1),'Max',d.colorScale(2), 'String',num2str(d.colorScale(1)));
   d.hOverlaySliderMin = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[overlaySlidersLeft+.1 .04 .3 .03],...
      'Min', d.colorScale(1),'Max',d.colorScale(2), 'value',d.colorScale(1));
   d.hOverlaySliderMax = [];
   d.hOverlaySliderMax = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[overlaySlidersLeft+.1 .01 .3 .03],...
      'Min', d.colorScale(1),'Max',d.colorScale(2), 'value',d.colorScale(2),'CreateFcn',{@changeOverlay,d});
   set(d.hOverlaySliderMin,'Callback',{@changeOverlay,d});
   set(d.hOverlayEditMin,'Callback',{@changeOverlay,d});
   set(d.hOverlaySliderMax,'Callback',{@changeOverlay,d});
   set(d.hOverlayEditMax,'Callback',{@changeOverlay,d});
   set(d.hOverlayList,'Callback',{@changeOverlay,d});
   set(d.hOverlayAlphaCheck,'Callback',{@changeOverlay,d});
   uicontrol('Parent',fignum, 'Unit','normalized', 'Style','Checkbox', 'Position',[.02 overlayButtonsBottom+.04 .20 .03],...
         'Callback',{@makeVisible,get(fignum,'userdata')}, 'String','Show Overlay','value',1);
   hOverlay = get(fignum,'userdata');
   set(hOverlay,'ambientStrength',.5);  %not too bright, not too dark
   set(hOverlay,'specularStrength',.3); %little reflectance
   set(gca,'Alim',[0 1])
   if computeContours
      set(d.hOverlayRestrictCheck,'Callback',{@changeOverlay,d});
   end
end

%Rois plot and control buttons
hold on
for i_roi = 1:length(thisView.ROIs)
   roiColor = color2RGB(thisView.ROIs(i_roi).color);
   temp = patch(roiD.EdgeXcoords{i_roi},roiD.EdgeYcoords{i_roi},roiD.EdgeZcoords{i_roi},roiColor,'edgecolor',roiColor,'LineWidth',2);
   if ~isempty(temp)
      roiD.hRoi(i_roi) = temp;
      uicontrol('Parent',fignum, 'Unit','normalized', 'Style','Checkbox',...
        'Position',[.02 roisButtonsBottom+roisCheckHeight*(length(thisView.ROIs) - i_roi+1) .15 roisCheckHeight],...
         'value',1, 'String', thisView.ROIs(i_roi).name, 'Callback', {@makeVisible,roiD.hRoi(i_roi)});
   end
end
uicontrol('Parent',fignum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.02 roisButtonsBottom .1 .03],...
      'Value',1, 'String', {'Perimeter','Grid','Opaque'},'Callback',{@drawRois,roiD});


daspect(1./basevoxelsize)
view(3);
axis tight;
grid on;
if ~isempty(analysis)
   colormap(d.colorMap);
   caxis(d.colorScale);
   colorbar('Units','normalized', 'OuterPosition',[.91 .5 .08 .45 ]);
end

xlabel('X');
ylabel('Y');
zlabel('Z');
set(gca,'Xdir','reverse'); %the X axis goes from right to left(?)
uicontrol('Parent',fignum, 'Unit','normalized', 'Style','checkbox', 'Position',[rotateViewLeft+.1 .04 .15 .03],...
         'Min', 0,'Max',1, 'Value',1, 'String', 'Reverse X Axis','Callback',{@reverseX});

%rotation control
hRotate = rotate3d;
set(hRotate,'Enable','on');
[azimuth,elevation] = view;
rotateViewSliderLength = .90-rotateViewLeft;
uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[rotateViewLeft .03 .05 .03],'String', 'Azimuth');
hAzimuthEdit = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[rotateViewLeft .01 .05 .03],...
      'String',num2str(azimuth));
hAzimuthSlider = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[rotateViewLeft+.05 .01 rotateViewSliderLength .03],...
      'Min', -360,'Max',360, 'value',azimuth);
uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[.94 .07+rotateViewSliderLength .05 .03],'String', 'Elevation');
hElevationEdit = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[.94 .05+rotateViewSliderLength .05 .03 ],...
      'String',num2str(elevation));
hElevationSlider = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[.95 .05 .03 rotateViewSliderLength],...
      'Min', -90,'Max',90, 'value',elevation);
set(hAzimuthSlider,'Callback',{@changeView,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit});
set(hElevationSlider,'Callback',{@changeView,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit});
set(hAzimuthEdit,'Callback',{@changeView,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit});
set(hElevationEdit,'Callback',{@changeView,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit});
set(hRotate,'ActionPostCallback',{@changeView,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit});

%set 2 lights at opposite sides of the object
h_light1 = camlight(90,45);
h_light2 = camlight(-90,-45);
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------- change and/or diplays the elevation and azimuth after mouse rotate (rotate3d)
function changeView(handle,eventdata,hAzimuthSlider,hElevationSlider,hAzimuthEdit,hElevationEdit)

if isnumeric(handle)
   azimuth = round(get(hAzimuthSlider,'Value'));
   elevation = round(get(hElevationSlider,'Value'));
   switch(handle)
      case hAzimuthSlider
         azimuth = round(get(handle,'Value'));
      case hElevationSlider
         elevation = round(get(handle,'Value'));
      case hAzimuthEdit
          azimuth = round(str2num(get(handle,'String')));
      case hElevationEdit
         elevation= round(str2num(get(handle,'String')));
   end
   view(azimuth,elevation);
else
   [azimuth,elevation] = view;
   azimuth = round(azimuth);
   elevation = round(elevation);
end

set(hAzimuthSlider,'Value',azimuth);
set(hElevationSlider,'Value',elevation);
set(hAzimuthEdit,'String',num2str(azimuth));
set(hElevationEdit,'String',num2str(elevation));

return;


%---------------------------------------------------- makes objects visible or invisible
function reverseX(handle,eventdata)

if get(handle,'value')
   set(gca,'xdir','reverse');
else
   set(gca,'xdir','normal');
end
return;


%---------------------------------------------------- makes objects visible or invisible
function makeVisible(handle,eventdata,h_axis)

if get(handle,'value')
   set(h_axis,'visible','on');
else
   set(h_axis,'visible','off');
end
return;


%---------------------------------------------------- draws Rois
function drawRois(handle,eventdata,roiD)

list = get(handle,'String');
switch(list{get(handle,'Value')})
   case 'Perimeter'
      Xcoords = roiD.EdgeXcoords;
      Ycoords = roiD.EdgeYcoords;
      Zcoords = roiD.EdgeZcoords;
   case 'Grid'
      Xcoords = roiD.GridXcoords;
      Ycoords = roiD.GridYcoords;
      Zcoords = roiD.GridZcoords;
   case 'Opaque'
      Xcoords = roiD.FacesXcoords;
      Ycoords = roiD.FacesYcoords;
      Zcoords = roiD.FacesZcoords;
end

for i_roi = 1:length(Xcoords)
   if ~isempty(Xcoords{i_roi})
      set(roiD.hRoi(i_roi),'Xdata',Xcoords{i_roi},'Ydata',Ycoords{i_roi},'Zdata',Zcoords{i_roi});
   end
end

      

%---------------------------------------------------- recomputes overlay object
function changeOverlay(handle,eventdata,d)

minThreshold = str2num(get(d.hOverlayEditMin,'string'));
maxThreshold = str2num(get(d.hOverlayEditMax,'string'));
if handle == d.hOverlaySliderMin
   minThreshold = get(handle,'value');
elseif ~isempty(d.hOverlaySliderMax) && handle == d.hOverlaySliderMax
   maxThreshold = get(handle,'value');
end

h_figure = get(handle,'parent');
figure(h_figure);
hOverlay = get(h_figure,'userdata');
useAlpha = get(d.hOverlayAlphaCheck,'Value');

list=get(d.hOverlayList,'string');
switch(list{get(d.hOverlayList,'Value')})
   
   case 'IsoContour'

      if get(d.hOverlayRestrictCheck,'Value')
         data = d.subVolOverlayDataRois;
      else
         data = d.subVolOverlayData;
      end
      faceColor = d.colorMap(round(size(d.colorMap,1)*(maxThreshold-d.colorScale(1))/diff(d.colorScale)),:);
      patchStruct = isosurface(d.blobsXcoords, d.blobsYcoords, d.blobsZcoords, data, maxThreshold); 

      if ~isempty(hOverlay)
         set(hOverlay,'vertices',patchStruct.vertices,'faces',patchStruct.faces,'FaceColor',faceColor,'FaceAlpha',1);
      else
         hOverlay = patch(patchStruct,'edgecolor','none','FaceColor',faceColor);
      end

   case 'Cubes'
      
      firstIndex = find(d.cubesCData>=minThreshold,1,'first');
      lastIndex = find(d.cubesCData<=maxThreshold,1,'last');
      if ~isempty(firstIndex) && ~isempty(lastIndex) && lastIndex>=firstIndex
         faceCData = d.cubesCData(firstIndex:lastIndex);
         faceAlphaData = d.cubesAlphaData(firstIndex:lastIndex);
         cubesBaseCoords = d.cubesBaseCoords(:,firstIndex:lastIndex);
         [faceXcoords, faceYcoords, faceZcoords,faceCData,faceAlphaData] = ...
               makeCubeFaces(cubesBaseCoords,faceCData,faceAlphaData,useAlpha == 0 || (d.alpha == 1 && all(faceAlphaData(:)==1)));
         %to use transparency, I need to convert X,Y,Z data to face/vertex data
         %patchStruct = surf2patch(faceXcoords,faceYcoords,faceZcoords,faceCData); %surf2patch does not handle color data correctly
         patchStruct.vertices = [reshape(faceXcoords,numel(faceXcoords),1) reshape(faceYcoords,numel(faceYcoords),1) reshape(faceZcoords,numel(faceZcoords),1)];
         patchStruct.faces = reshape((1:size(patchStruct.vertices,1)),[4 size(faceXcoords,2)])';
      else
         patchStruct.vertices = [];
         patchStruct.faces = [];
         faceCData = 1;
         faceAlphaData = 1;
         firstIndex = Inf;
         lastIndex = -Inf;
      end
      
      if ~isempty(hOverlay)
         set(hOverlay,'vertices',patchStruct.vertices,'faces',patchStruct.faces,'FaceVertexCData',faceCData','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat');
      else
         patchStruct.facevertexcdata = faceCData';
         hOverlay = patch(patchStruct,'edgecolor','none','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat');
      end
      
      %count number of visible voxels in each ROI
      for i_roi = 1:length(d.roiSize)
         set(d.hRoiVoxelNum(i_roi),'String',[num2str(sum(d.roiDataIndex{i_roi}>=firstIndex & d.roiDataIndex{i_roi}<=lastIndex)) '/' num2str(d.roiSize(i_roi))]);
      end
      
end

set(get(handle,'parent'),'userdata',hOverlay);
if ~isempty(d.hOverlaySliderMax) %this is empty when we call this callback by CreateFcn
   set(d.hOverlaySliderMax,'value',maxThreshold);
end
set(d.hOverlayEditMax,'String',num2str(maxThreshold));
set(d.hOverlaySliderMin,'Value',minThreshold);
set(d.hOverlayEditMin,'String',num2str(minThreshold));

return;


%--------------------------------------returns the coordinates of the faces of cubes centered on voxel coordinates
function [facesXcoords, facesYcoords, facesZcoords, colorData, alphaData] = makeCubeFaces(voxelCoords,colorData,alphaData,opaque)

cubeX = [-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5]';
cubeY = [-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5]';
cubeZ = [-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5]';


%coordinates of a cube around each remaining voxel
voxelsNumber = size(voxelCoords,2);
cubesXcoords = repmat(shiftdim(voxelCoords(1,:),-1),[4 6 1]) + repmat(cubeX,[1 1 voxelsNumber]);
cubesYcoords = repmat(shiftdim(voxelCoords(2,:),-1),[4 6 1]) + repmat(cubeY,[1 1 voxelsNumber]);
cubesZcoords = repmat(shiftdim(voxelCoords(3,:),-1),[4 6 1]) + repmat(cubeZ,[1 1 voxelsNumber]);

facesXcoords = reshape(cubesXcoords,4,6*voxelsNumber);
facesYcoords = reshape(cubesYcoords,4,6*voxelsNumber);
facesZcoords = reshape(cubesZcoords,4,6*voxelsNumber);

if ~ieNotDefined('colorData')
   if size(colorData,2)==1
      colorData = colorData';
   end
   colorData = reshape(repmat(colorData,[6 1]),1,6*voxelsNumber);
end
if ~ieNotDefined('alphaData')
   if size(alphaData,2)==1
      alphaData = alphaData';
   end
   alphaData = reshape(repmat(alphaData',[6 1]),1,6*voxelsNumber);
end

if opaque   %if everything is opaque, then we can just keep the outside faces
   facesCoords = [facesXcoords; facesYcoords; facesZcoords]';
   [dump, unique_indices] = unique(facesCoords,'rows');
   double_indices = setdiff(1:size(facesCoords,1),unique_indices);

   [unique_facesCoords, unique_indices] = setdiff(facesCoords,facesCoords(double_indices,:),'rows');
   if ~ieNotDefined('colorData')
      colorData = colorData(unique_indices);
   end
   if ~ieNotDefined('alphaData')
      alphaData = ones(size(colorData));
   end

   if ~isempty(unique_facesCoords)
      facesXcoords = unique_facesCoords(:,1:4)';
      facesYcoords = unique_facesCoords(:,5:8)';
      facesZcoords = unique_facesCoords(:,9:12)';
   end
end


