% renderRois3D - displays 3D rendering of currently loaded ROIs
%
%        $Id$
%      usage: [  ] = renderRois3D(thisView,overlayList,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2010-02-15
%     inputs: 
%    outputs: 
%
%    purpose: displays 3D rendering of currently loaded ROIs
%
% planned improvements: take a subvolume (slightly larger than ROIs) from the start to reduce computing time, especially
%                       when resampling to base anatomy.

function renderRois3D(thisView,overlayList,scanNum,x,y,z,roi)

computeContours = 0;

roiList = viewGet(thisView,'visibleRois');
nRois = length(roiList);

if ~nRois
   mrWarnDlg('(renderRois3D) Please make at least one ROI visible');
   return;
end

% if length(overlayList)>1
%    mrWarnDlg('(renderRois3D) This is not implemented when several overlays are displayed');
%    return;
% end

fprintf(1,'Loading data...');
tic
baseNum = viewGet(thisView,'currentBase');
basevoxelsize = viewGet(thisView,'basevoxelsize');
d.alpha = viewGet(thisView,'alpha');
base_size = viewGet(thisView,'basedims',baseNum);

% get the analysis structure
analysis = viewGet(thisView,'analysis');
currentAnalysis = viewGet(thisView,'currentAnalysis');
fprintf(1,'Done\n');
toc


%--------------------------------Mask and Alpha overlay data -------------------------------
if ~isempty(currentAnalysis)
   
  alphaOverlayNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay'));

  %First, get a mask of non-zero voxel representing the current overlay display
  %This is taken from computeOverlay.m
  fprintf(1,'Computing overlay mask...');
  [maskData, overlayData]= maskOverlay(thisView,[overlayList alphaOverlayNum],scanNum);
  maskData = maskData{1};
  overlayData = overlayData{1};
  fprintf(1,'Done\n');
  d.nOverlays = length(overlayList);

%    overlayData = getBaseSpaceOverlay(thisView, overlayData, scanNum, baseNum);
% 
%    maskData = getBaseSpaceOverlay(thisView, double(maskData), scanNum, baseNum,'nearest');
  d.base2scan = viewGet(thisView,'base2scan');
  maskData(isnan(maskData))=0;
  maskData = logical(maskData);

   %make alpha map
   if isempty(alphaOverlayNum)
      alphaOverlayData = d.alpha*ones(size(maskData));
   else
      fprintf(1,'Computing alpha overlay...');
      
      alphaOverlayData = overlayData(:,:,:,d.nOverlays+1:end);
      overlayData = overlayData(:,:,:,1:d.nOverlays);
      
      alphaMaskData = maskData(:,:,:,d.nOverlays+1:end);
      maskData = maskData(:,:,:,1:d.nOverlays);
      
      alphaOverlayData(~alphaMaskData)=NaN;
      
      % get the range of the alpha overlay 
      cOverlay=0;
      for iOverlay = alphaOverlayNum;
        cOverlay=cOverlay+1;
        range = viewGet(thisView,'overlayColorRange',iOverlay);
        % handle setRangeToMax (to debug)
        if strcmp(viewGet(thisView,'overlayCtype',alphaOverlayNum),'setRangeToMax')
          clip = viewGet(thisView,'overlayClip',iOverlay);
         maxRange = max(clip(1),min(alphaOverlayData(alphaMaskData)));
         if ~isempty(maxRange),range(1) = maxRange;end
         minRange = min(max(alphaOverlayData(alphaMaskData)),clip(2));
         if ~isempty(minRange),range(2) = minRange;end
        end
        % now compute the alphaOverlay as a number from
        % 0 to 1 of the range
        thisAlphaOverlayData = d.alpha*((alphaOverlayData(:,:,:,cOverlay)-range(1))./diff(range));
        thisAlphaOverlayData(thisAlphaOverlayData>d.alpha) = d.alpha;
        thisAlphaOverlayData(thisAlphaOverlayData<0) = 0;
        alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
        if alphaOverlayExponent<0   %JB if the alpha overlay exponent is negative, set it positive and take the complementary values to 1 for the alpha map
          alphaOverlayExponent = -alphaOverlayExponent;
          thisAlphaOverlayData = 1-thisAlphaOverlayData;
        end
        alphaOverlayData(:,:,:,cOverlay) = thisAlphaOverlayData.^alphaOverlayExponent;
      end
      alphaOverlayData(isnan(alphaOverlayData)) = 0;
   end   
   fprintf(1,'Done\n');

else
   maskData = [];
   alphaOverlayData = [];
end



%------------------------------------ Construct ROI objects ------------------------------------------------%
fprintf(1,'Constructing 3D ROIs...');
tic
edgeX = [-.5 -.5; -.5 -.5; -.5 -.5; -.5 -.5; -.5 .5; -.5 .5; -.5 .5; -.5 .5; .5 .5; .5 .5; .5 .5; .5 .5];
edgeY = [-.5 .5; -.5 .5; -.5 -.5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5];
edgeZ = [-.5 -.5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; -.5 .5; -.5 .5];
cubeEdgesCoords = [edgeX edgeY edgeZ];

allRoisCoords = [];
for iRoi = 1:nRois

  roi = viewGet(thisView, 'roi',roiList(iRoi));
%   base2roi = viewGet(thisView, 'base2roi',roiList(iRoi));
%   roiCoords = xformROIcoords(roi.coords,inv(base2roi),roi.voxelSize,basevoxelsize);
%   %past 2 lines equivalent to:
%   roiCoords = getROICoordinates(thisView,roiList(iRoi),0,[],baseNum);
  roiCoords = getROICoordinates(thisView,roiList(iRoi));
  roiColor(iRoi,:) = color2RGB(roi.color);
  roiName{iRoi} = roi.name;
  
  if ~isempty(roiCoords)
    % we need to remove any coordinate that might fall outside the base anatomy
    outside_voxels = roiCoords(1,:)<1 | roiCoords(1,:)>base_size(1) |...
                         roiCoords(2,:)<1 | roiCoords(2,:)>base_size(2) |...
                         roiCoords(3,:)<1 | roiCoords(3,:)>base_size(3) ;
    roiCoords(:,outside_voxels) = [];
  end
  if ~isempty(roiCoords)
    d.roiSize(iRoi) = size(roiCoords,2);
    %remember which voxels are in which rois
    d.roiDataIndex{iRoi} = size(allRoisCoords,2)+(1:d.roiSize(iRoi));
    %put Roi coords together to simplify overlay computation
    allRoisCoords = [allRoisCoords roiCoords];

    %compute coordinates of all cube edges around each voxel
    edgesCoords = roiCoords([1 1 2 2 3 3],:);
    edgesCoords = reshape(repmat(edgesCoords,12,1),6,12*size(roiCoords,2))';
    edgesCoords = edgesCoords+repmat(cubeEdgesCoords,size(roiCoords,2),1);
    %find unique edges around  this roi
    [gridCoords, uniqueIndices] = unique(edgesCoords,'rows');
    doubleIndices = setdiff(1:size(edgesCoords,1),uniqueIndices);

    uniqueEdgeCoords = setdiff(edgesCoords,edgesCoords(doubleIndices,:),'rows');

    roiD.GridXcoords{iRoi} = gridCoords(:,1:2)';
    roiD.GridYcoords{iRoi} = gridCoords(:,3:4)';
    roiD.GridZcoords{iRoi} = gridCoords(:,5:6)';
    roiD.EdgeXcoords{iRoi} = uniqueEdgeCoords(:,1:2)';
    roiD.EdgeYcoords{iRoi} = uniqueEdgeCoords(:,3:4)';
    roiD.EdgeZcoords{iRoi} = uniqueEdgeCoords(:,5:6)';
    
    [roiD.GridXcoords{iRoi},roiD.GridYcoords{iRoi},roiD.GridZcoords{iRoi}] = convert2baseCoords(...
      roiD.GridXcoords{iRoi},roiD.GridYcoords{iRoi},roiD.GridZcoords{iRoi},d.base2scan);
    [roiD.EdgeXcoords{iRoi},roiD.EdgeYcoords{iRoi},roiD.EdgeZcoords{iRoi}] = convert2baseCoords(...
      roiD.EdgeXcoords{iRoi},roiD.EdgeYcoords{iRoi},roiD.EdgeZcoords{iRoi},d.base2scan);
    
    %Faces (opaque ROi)
    %coordinates of a cube around each voxel
    [roiD.FacesXcoords{iRoi}, roiD.FacesYcoords{iRoi}, roiD.FacesZcoords{iRoi}] = makeCubeFaces(roiCoords,[],[],1,d.base2scan);
    %Faces (opaque ROi));

   else
      roiD.GridXcoords{iRoi} = [];
      roiD.GridYcoords{iRoi} = [];
      roiD.GridZcoords{iRoi} = [];
      roiD.EdgeXcoords{iRoi} = [];
      roiD.EdgeYcoords{iRoi} = [];
      roiD.EdgeZcoords{iRoi} = [];
      roiD.FacesXcoords{iRoi} = [];
      roiD.FacesYcoords{iRoi} = [];
      roiD.FacesZcoords{iRoi} = [];
      d.roiSize(iRoi) = 0;
      d.roiDataIndex{iRoi}=[];
   end
end
fprintf(1,'Done\n');
toc


%------------------------------------ Construct 3D Overlay  ------------------------------------------------%
if ~isempty(currentAnalysis)

  %find indices of unique voxels in all rois
  [allRoisCoords, dump, duplicateVoxels] = unique(allRoisCoords','rows');
  allRoisCoords = allRoisCoords';
  for iRoi = 1:nRois
    d.roiDataIndex{iRoi} = duplicateVoxels(d.roiDataIndex{iRoi});
  end

  %%%%%%%%%%%%%%%%%%% Compute cube faces
  fprintf(1,'Constructing overlay data cubes...');
  tic
  cOverlay = 0;
  for iOverlay=overlayList
    cOverlay = cOverlay+1;
    %d.colorMap = viewGet(thisView,'colormap',iOverlay);
    colorMap = analysis.overlays(iOverlay).colormap;
    colorRange = str2num(num2str(viewGet(thisView,'overlayColorRange',iOverlay))); %str2num(num2str()) is because the edit box rounds the values  
                                                                           %and causes problems if the min/max of the sliders is not rounded
    RGB = rescale2rgb(overlayData(:,:,:,cOverlay),colorMap,colorRange);
    if d.nOverlays>1
      thisAlphaData = alphaOverlayData(:,:,:,cOverlay);
      % apply the mask to alpha data
      thisAlphaData(~maskData(:,:,:,cOverlay))=NaN;
    % 1) pre-multiply colors by alpha
      RGB = repmat(thisAlphaData,[1 1 1 3]).*RGB;
      if cOverlay==1
        d.cubesRGB=RGB;
        d.cubesAlphaData = thisAlphaData;
      else
      % 2) add pre-multiplied colormaps (and alpha channels)
        d.cubesRGB = d.cubesRGB.*(1-RGB)+ RGB;
        d.cubesAlphaData = d.cubesAlphaData .* (1-thisAlphaData)+thisAlphaData;
      end
    else
      d.cubesRGB = RGB;
      d.cubesAlphaData = alphaOverlayData(:,:,:,cOverlay);
    end
  end
  d.cubesRGB = reshape(d.cubesRGB,[numel(d.cubesAlphaData),3]);
  
  %keep only voxels in Rois
  cubesOverlayDataIndex = sub2ind(size(d.cubesAlphaData), allRoisCoords(1,:)', allRoisCoords(2,:)', allRoisCoords(3,:)' );
  d.cubesRGB = d.cubesRGB(cubesOverlayDataIndex,:);
  d.cubesAlphaData = d.cubesAlphaData(cubesOverlayDataIndex);

  if d.nOverlays==1
    clipRange = str2num(num2str(viewGet(thisView,'overlayClip'))); 
    d.cubesColorIndex = overlayData(cubesOverlayDataIndex);
    %mask the data
    maskData = maskData(cubesOverlayDataIndex);
    d.cubesAlphaData(~maskData)=NaN;
    %exclude Nan overlay values
    voxelsToKeepIndices = find(~isnan(d.cubesAlphaData));
    d.cubesColorIndex = d.cubesColorIndex(voxelsToKeepIndices);
    %sort faces by overlay data value, that will make things faster later
    [d.cubesColorIndex, sortIndex] = sort(d.cubesColorIndex);
    voxelsToKeepIndices = voxelsToKeepIndices(sortIndex);
    d.cubesAlphaData = d.cubesAlphaData(voxelsToKeepIndices);
  else
    %here we divide by the computed alpha, because there is no base to multiply it with
    d.cubesRGB = d.cubesRGB./repmat(d.cubesAlphaData,[1 3]);
    %exclude Nan overlay values
    voxelsToKeepIndices = find(any(~isnan(d.cubesRGB),2));
    %sort according to alpha values, this will be important when removing duplicate cube faces
    d.cubesAlphaData = d.cubesAlphaData(voxelsToKeepIndices);
    [d.cubesAlphaData, sortIndex] = sort(d.cubesAlphaData);
    voxelsToKeepIndices = voxelsToKeepIndices(sortIndex);
  end

  %apply excluding and/or sorting to other RGB and alpha data arrays
  d.cubesRGB = d.cubesRGB(voxelsToKeepIndices,:);
  d.cubesBaseCoords = allRoisCoords(:,voxelsToKeepIndices);
  for iRoi = 1:nRois
    [dump, d.roiDataIndex{iRoi}] = intersect(voxelsToKeepIndices,d.roiDataIndex{iRoi});
  end

  fprintf(1,'Done\n');
  toc
      
%%%%%%%%%%%%%%%%%%% Compute upsampled smooth volumes for blobs (hasn't been debugged after important changes in this function)
   if computeContours
      fprintf(1,'Smoothing and upsampling overlay data...');
      %reduce volume
      minX = min(allRoisCoords(1,:));
      maxX = max(allRoisCoords(1,:));
      minY = min(allRoisCoords(2,:));
      maxY = max(allRoisCoords(2,:));
      minZ = min(allRoisCoords(3,:));
      maxZ = max(allRoisCoords(3,:));

      %volume slightly larger for smoothing
      margin = 3;
      minXmargin = max(minX-margin,1);
      maxXmargin = min(maxX+margin,size(overlayData,1)); 
      minYmargin = max(minY-margin,1);                         
      maxYmargin = min(maxY+margin,size(overlayData,2));   
      minZmargin = max(minZ-margin,1);
      maxZmargin = min(maxZ+margin,size(overlayData,3));

      %reduce volume
      [Ycoords,Xcoords,Zcoords,d.subVolOverlayData] = subvolume(overlayData,[minYmargin maxYmargin minXmargin maxXmargin minZmargin maxZmargin]); %subvolume swaps X and Y axes...
      [Ycoords,Xcoords,Zcoords,subMaskData] = subvolume(maskData,[minYmargin maxYmargin minXmargin maxXmargin minZmargin maxZmargin]); %subvolume swaps X and Y axes...
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
      blobsRoisBaseCoords = xformROIcoords(allRoisCoords,xform,baseVoxelSize,baseVoxelSize/upsamplingFactor);
      blobsRoisBaseCoords = blobsRoisBaseCoords(1:3,:) - repmat(upsamplingFactor*[minXmargin;minYmargin;minZmargin],1,size(blobsRoisBaseCoords,2))+1;
      blobsRoisBaseCoordsIndex = sub2ind(size(d.subVolOverlayData),blobsRoisBaseCoords(1,:),blobsRoisBaseCoords(2,:),blobsRoisBaseCoords(3,:));
      d.subVolOverlayDataRois = NaN(size(d.subVolOverlayData));
      d.subVolOverlayDataRois(blobsRoisBaseCoordsIndex) = d.subVolOverlayData(blobsRoisBaseCoordsIndex);
      [d.blobsYcoords,d.blobsXcoords,d.blobsZcoords,d.subVolOverlayDataRois] = subvolume(yi, xi, zi, d.subVolOverlayDataRois,[minY maxY minX maxX minZ maxZ]); %subvolume swaps X and Y axes...
      [d.blobsYcoords,d.blobsXcoords,d.blobsZcoords,d.subVolOverlayData] = subvolume(yi, xi, zi, d.subVolOverlayData,[minY maxY minX maxX minZ maxZ]); %subvolume swaps X and Y axes...

      fprintf(1,'Done\n');
   end
end
      
%------------------------------------ Get Base data
baseType = viewGet(thisView,'baseType',baseNum);
if ~baseType
  allRoisBaseCoords = d.base2scan\[allRoisCoords;ones(1,size(allRoisCoords,2))];
  minCoords = floor(min(allRoisBaseCoords,[],2)');
  maxCoords = ceil(max(allRoisBaseCoords,[],2)');
  baseDims = viewGet(thisView,'baseDims',baseNum); 
  minCoords = max(minCoords(1:3),[1 1 1]);    %this is probably not needed because voxels outside the base should have been excluded
  maxCoords = min(maxCoords(1:3),baseDims);   %(but might change if transformation to base is applied further down the line)
  baseColormap = gray(256);
  baseClip = viewGet(thisView,'baseClip',baseNum);
  baseGamma = viewGet(thisView,'baseGamma',baseNum);
  for iDim=1:3
    otherDims = setdiff(1:3,iDim);
    [baseImage,baseCoords{iDim},baseCoordsHomogeneous{iDim}] = ...
      getBaseSlice(thisView,minCoords(iDim),iDim,0,baseNum,baseType);
    baseRGB{iDim} = rescale2rgb(baseImage,baseColormap,baseClip,baseGamma);
    thisMinCoords = minCoords(otherDims);
    thisMaxCoords = maxCoords(otherDims);
    baseRGB{iDim} = baseRGB{iDim}(thisMinCoords(1):thisMaxCoords(1),thisMinCoords(2):thisMaxCoords(2),:);
    baseRGB{iDim}(:,end+1,:) = baseRGB{iDim}(:,end,:); %Add one column
    baseRGB{iDim}(end+1,:,:) = baseRGB{iDim}(end,:,:); %Add one row
    baseCoords{iDim} = baseCoords{iDim}(thisMinCoords(1):thisMaxCoords(1),thisMinCoords(2):thisMaxCoords(2),:);
    baseCoords{iDim}(:,end+1,:) = baseCoords{iDim}(:,end,:); %Add one column
    baseCoords{iDim}(:,end,otherDims(2)) = baseCoords{iDim}(:,end,otherDims(2))+1; %Add one column
    baseCoords{iDim}(end+1,:,:) = baseCoords{iDim}(end,:,:); %Add one row
    baseCoords{iDim}(end,:,otherDims(1)) = baseCoords{iDim}(end,:,otherDims(1))+1; %Add one column
    baseCoords{iDim} = baseCoords{iDim} - .5; %shift the coords by .5 so that each square will be centered on the integer coordinate
  end
end
  
  
%---------------------- construct figure -----------------%
% selectGraphWin;
% global MLR;
% fignum = MLR.graphFigure;
fignum = selectGraphWin;
% turn off menu/title etc.
set(fignum,'NumberTitle','off','Name','RenderRois3D','Toolbar','Figure','color',[.7 .7 .7]);
roisButtonsBottom = .33;
overlayButtonsBottom = .1;
roisCheckHeight = min(.05,(.9 - roisButtonsBottom)/nRois);
overlaySlidersLeft=.02;
rotateViewLeft = .6;

for iRoi = 1:nRois
   d.hRoiVoxelNum(iRoi) = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',...
     [.04 roisButtonsBottom-.01+roisCheckHeight*(nRois-iRoi+1) .08 roisCheckHeight/2]);
end

%Overlay control buttons
if ~isempty(currentAnalysis)
  overlayNames = viewGet(thisView,'overlayNames');
   set(fignum,'Name',['RenderRois3D - ' viewGet(thisView,'homeDir') ' - ' viewGet(thisView,'analysisName') ' - ' overlayNames{overlayList} ]);

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
      'Min', 0,'Max',1, 'Value',1, 'String','Use Alpha');
   
   %Overlay sliders
   if d.nOverlays==1
     uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[overlaySlidersLeft .04 .13 .03],'String','Min Threshold');
     uicontrol('Parent',fignum, 'Unit','normalized', 'Style','text', 'Position',[overlaySlidersLeft .01 .13 .03],'String','Max Threshold');
     d.hOverlayEditMax = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[overlaySlidersLeft+.4 .01 .1 .03],...
        'Min', clipRange(1),'Max',clipRange(2), 'String',num2str(clipRange(2)));
     d.hOverlayEditMin = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','edit', 'Position',[overlaySlidersLeft+.4 .04 .1 .03],...
        'Min', clipRange(1),'Max',clipRange(2), 'String',num2str(clipRange(1)));
     d.hOverlaySliderMin = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[overlaySlidersLeft+.1 .04 .3 .03],...
        'Min', clipRange(1),'Max',clipRange(2), 'value',clipRange(1));
     d.hOverlaySliderMax = [];
     d.hOverlaySliderMax = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','slider', 'Position',[overlaySlidersLeft+.1 .01 .3 .03],...
        'Min', clipRange(1),'Max',clipRange(2), 'value',clipRange(2));%,'CreateFcn',{@changeOverlay,d});
     set(d.hOverlaySliderMin,'Callback',{@changeOverlay,d});
     set(d.hOverlayEditMin,'Callback',{@changeOverlay,d});
     set(d.hOverlaySliderMax,'Callback',{@changeOverlay,d});
     set(d.hOverlayEditMax,'Callback',{@changeOverlay,d});
   end
   set(d.hOverlayAlphaCheck,'Callback',{@changeOverlay,d});
   set(d.hOverlayList,'Callback',{@changeOverlay,d});
   hOverlay = get(fignum,'userdata');
   set(hOverlay,'ambientStrength',.5);  %not too bright, not too dark
   set(hOverlay,'specularStrength',.3); %little reflectance
   if computeContours
      set(d.hOverlayRestrictCheck,'Callback',{@changeOverlay,d});
   end
   changeOverlay(d.hOverlayAlphaCheck,[],d);
   %this has to be after calling cahngeOverlay, because we want to use the hande of the overlay cubes
   uicontrol('Parent',fignum, 'Unit','normalized', 'Style','Checkbox', 'Position',[.02 overlayButtonsBottom+.04 .20 .03],...
         'Callback',{@makeVisible,get(fignum,'userdata')}, 'String','Show Overlay','value',1);

end

hAxes = gca;
set(hAxes,'Alim',[0 1],'color','none')
%Rois plot and control buttons
hold on
for iRoi = 1:nRois
   temp = patch(roiD.EdgeXcoords{iRoi},roiD.EdgeYcoords{iRoi},roiD.EdgeZcoords{iRoi},roiColor(iRoi,:),'edgecolor',roiColor(iRoi,:),'LineWidth',2,'parent',hAxes);
   if ~isempty(temp)
      roiD.hRoi(iRoi) = temp;
      uicontrol('Parent',fignum, 'Unit','normalized', 'Style','Checkbox',...
        'Position',[.02 roisButtonsBottom+roisCheckHeight*(nRois - iRoi+1) .15 roisCheckHeight],...
         'value',1, 'String', roiName{iRoi}, 'Callback', {@makeVisible,roiD.hRoi(iRoi)});
   end
end
uicontrol('Parent',fignum, 'Unit','normalized', 'Style','popupmenu', 'Position',[.02 roisButtonsBottom .1 .03],...
      'Value',1, 'String', {'Perimeter','Grid','Opaque'},'Callback',{@drawRois,roiD});

% plot Base slices
if ~baseType
  for iDim = 1:3
    hSurface(iDim) = surface(baseCoords{iDim}(:,:,1),baseCoords{iDim}(:,:,2),baseCoords{iDim}(:,:,3),baseRGB{iDim}, 'EdgeColor','none','parent',hAxes);
  end
  d.hOverlayShowBase = uicontrol('Parent',fignum, 'Unit','normalized', 'Style','checkbox',...
    'Position',[.02 overlayButtonsBottom+.18 .15 .03],...
     'Min', 0,'Max',1, 'Value',1, 'String', 'Show Base', 'callback',{@makeVisible,hSurface});
end

daspect(1./basevoxelsize)
view(3);
axis tight;
grid on;
% if ~isempty(currentAnalysis)
%    colormap(d.colorMap);
%    caxis(d.colorRange);
%    colorbar('Units','normalized', 'OuterPosition',[.91 .5 .08 .45 ]);
% end

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

for iRoi = 1:length(Xcoords)
   if ~isempty(Xcoords{iRoi})
      set(roiD.hRoi(iRoi),'Xdata',Xcoords{iRoi},'Ydata',Ycoords{iRoi},'Zdata',Zcoords{iRoi});
   end
end

      

%---------------------------------------------------- recomputes overlay object
function changeOverlay(handle,eventdata,d)

if d.nOverlays==1
  minThreshold = str2num(get(d.hOverlayEditMin,'string'));
  maxThreshold = str2num(get(d.hOverlayEditMax,'string'));
  if handle == d.hOverlaySliderMin
     minThreshold = get(handle,'value');
  elseif ~isempty(d.hOverlaySliderMax) && handle == d.hOverlaySliderMax
     maxThreshold = get(handle,'value');
  elseif handle == d.hOverlayEditMin || handle == d.hOverlayEditMax
    sliderMin = get(d.hOverlaySliderMin,'min');
    sliderMax = get(d.hOverlaySliderMax,'max');
    if  minThreshold<sliderMin
      minThreshold = sliderMin;
    end
    if  maxThreshold>sliderMax
      maxThreshold = sliderMax;
    end
  end
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
      faceColor = d.colorMap(round(size(d.colorMap,1)*(maxThreshold-d.colorRange(1))/diff(d.colorRange)),:);  
      patchStruct = isosurface(d.blobsXcoords, d.blobsYcoords, d.blobsZcoords, data, maxThreshold); 

      if ~isempty(hOverlay)
         set(hOverlay,'vertices',patchStruct.vertices,'faces',patchStruct.faces,'FaceColor',faceColor,'FaceAlpha',1);
      else
         hOverlay = patch(patchStruct,'edgecolor','none','FaceColor',faceColor);
      end

   case 'Cubes'
      if d.nOverlays==1
        firstIndex = find(d.cubesColorIndex>=minThreshold,1,'first');
        lastIndex = find(d.cubesColorIndex<=maxThreshold,1,'last');
      else
        firstIndex=1;
        lastIndex=size(d.cubesAlphaData,1);
      end
      if ~isempty(firstIndex) && ~isempty(lastIndex) && lastIndex>=firstIndex
         %faceCData = d.cubesColorIndex(firstIndex:lastIndex);
         faceCData = d.cubesRGB(firstIndex:lastIndex,:);
         faceAlphaData = d.cubesAlphaData(firstIndex:lastIndex);
         cubesBaseCoords = d.cubesBaseCoords(:,firstIndex:lastIndex);
         [faceXcoords, faceYcoords, faceZcoords,faceCData,faceAlphaData] = ...
               makeCubeFaces(cubesBaseCoords,faceCData',faceAlphaData',useAlpha == 0 || (d.alpha == 1 && all(faceAlphaData(:)==1)),d.base2scan);
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
         set(hOverlay,'vertices',patchStruct.vertices,'faces',patchStruct.faces,...
           'FaceVertexCData',faceCData','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat');
      else
         patchStruct.facevertexcdata = faceCData';
         hOverlay = patch(patchStruct,'edgecolor','none','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat');
      end
      
      %count number of visible voxels in each ROI
      for iRoi = 1:length(d.roiSize)
         set(d.hRoiVoxelNum(iRoi),'String',[num2str(sum(d.roiDataIndex{iRoi}>=firstIndex & d.roiDataIndex{iRoi}<=lastIndex)) '/' num2str(d.roiSize(iRoi))]);
      end
      
end

set(get(handle,'parent'),'userdata',hOverlay);
if d.nOverlays==1
  set(d.hOverlaySliderMax,'value',maxThreshold);
  set(d.hOverlayEditMax,'String',num2str(maxThreshold));
  set(d.hOverlaySliderMin,'Value',minThreshold);
  set(d.hOverlayEditMin,'String',num2str(minThreshold));
end

return;


%--------------------------------------returns the coordinates of the faces of cubes centered on voxel coordinates
function [facesXcoords, facesYcoords, facesZcoords, colorData, alphaData] = makeCubeFaces(voxelCoords,colorData,alphaData,opaque,base2scan)

%coordinates of a cube around each remaining voxel
cubeX = [-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5]';
cubeY = [-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5]';
cubeZ = [-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5]';

voxelsNumber = size(voxelCoords,2);
cubesXcoords = repmat(shiftdim(voxelCoords(1,:),-1),[4 6 1]) + repmat(cubeX,[1 1 voxelsNumber]);
cubesYcoords = repmat(shiftdim(voxelCoords(2,:),-1),[4 6 1]) + repmat(cubeY,[1 1 voxelsNumber]);
cubesZcoords = repmat(shiftdim(voxelCoords(3,:),-1),[4 6 1]) + repmat(cubeZ,[1 1 voxelsNumber]);

facesXcoords = reshape(cubesXcoords,4,6*voxelsNumber);
facesYcoords = reshape(cubesYcoords,4,6*voxelsNumber);
facesZcoords = reshape(cubesZcoords,4,6*voxelsNumber);

if ~ieNotDefined('colorData')
  colorData = reshape(repmat(colorData,[6 1]),size(colorData,1),6*voxelsNumber);
end
if ~ieNotDefined('alphaData')
  alphaData = reshape(repmat(alphaData,[6 1]),1,6*voxelsNumber);
end

%find duplicate faces to remove (one of each duplicate)
facesCoords = [facesXcoords; facesYcoords; facesZcoords]';
[uniqueFaceCoords, uniqueIndices] = unique(facesCoords,'rows');

if opaque   %if everything is opaque, then we can just keep the outside faces, that is remove the other set of duplicate faces
   doubleIndices = setdiff(1:size(facesCoords,1),uniqueIndices);
   [uniqueFaceCoords, uniqueIndices] = setdiff(facesCoords,facesCoords(doubleIndices,:),'rows');
end
if ~ieNotDefined('colorData')
  colorData = colorData(:,uniqueIndices);
end
if ~opaque && ~ieNotDefined('alphaData')
  alphaData = alphaData(uniqueIndices);
  %round any alpha value >.985 to 1, because of a bug in the display
  alphaData(alphaData>.985)=1;
else
  alphaData = ones(size(uniqueIndices))';
end

if ~isempty(uniqueFaceCoords)
  facesXcoords = uniqueFaceCoords(:,1:4)';
  facesYcoords = uniqueFaceCoords(:,5:8)';
  facesZcoords = uniqueFaceCoords(:,9:12)';
end

[facesXcoords,facesYcoords,facesZcoords] = convert2baseCoords(facesXcoords,facesYcoords,facesZcoords,base2scan);

%-----------------------
function [baseX,baseY,baseZ] = convert2baseCoords(X,Y,Z,base2scan)

baseX = reshape(X,1,numel(X));
baseY = reshape(Y,1,numel(Y));
baseZ = reshape(Z,1,numel(Z));

baseCoords = base2scan\[baseX;baseY;baseZ;ones(1,size(baseX,2))];

baseX = baseCoords(1,:);
baseY = baseCoords(2,:);
baseZ = baseCoords(3,:);

baseX = reshape(baseX,size(X));
baseY = reshape(baseY,size(Y));
baseZ = reshape(baseZ,size(Z));



