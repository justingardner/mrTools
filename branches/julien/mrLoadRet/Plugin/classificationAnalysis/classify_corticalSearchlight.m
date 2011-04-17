function [out,params]=classify_corticalSearchlight(thisView,d,params,scanParams)
%This fucntion does a searchlight classification with a defined radius in
%an area defined by an ROI, which is normally a mask of the skullstripped
%brain. Has to be the run from the root directory of an MLR session.

out=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the mesh and corresponding an anatomy file
if exist('Segmentation')==7
  pathname = 'Segmentation';
else
  pathname = mrGetPref('volumeDirectory');
end

baseNum=[];
while isempty(baseNum)
  if fieldIsNotDefined(params,'canonicalBaseName')
    [filename,pathname] = uigetfile([pathname '/*.hdr'], 'Choose Base Anatomy');
    if isnumeric(filename),  return; end
    params.canonicalBaseName = [pathname filename];
  end
  baseHdr = cbiReadNiftiHeader(params.canonicalBaseName);
  [~,basename,~] = fileparts(baseHdr.hdr_name);
  baseNum = viewGet(thisView,'basenum',basename);
  if isempty(baseNum)
    if strcmp(questdlg('The canonical base needs to be loaded in the view','Base is not loaded in view','Choose another file','Abort','Choose another file'),'Abort')
      return;
    else
      params.canonicalBaseName = [];
    end
  end
end

% load the appropriate surface files
if fieldIsNotDefined(params,'innerSurfaceName')
  [filename,pathname] = uigetfile([pathname '/*WM.off'], 'Choose Inner Surface');
  if isnumeric(filename),  return; end
  params.innerSurfaceName = [pathname filename];
end
surface.inner = loadSurfOFF(params.innerSurfaceName);
if isempty(surface.inner)
  mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName)));
  return
end

surface.inner = xformSurfaceWorld2Array(surface.inner, baseHdr);
if fieldIsNotDefined(params,'outerSurfaceName')
  [filename,pathname] = uigetfile([pathname '/*GM.off'], 'Choose Outer Surface');
  if isnumeric(filename),  return; end
  params.outerSurfaceName = [pathname filename];
end
surface.outer = loadSurfOFF(params.outerSurfaceName);
if isempty(surface.outer)
  mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName)));
  return
end
surface.outer = xformSurfaceWorld2Array(surface.outer, baseHdr);

% build up a mrMesh-style structure, taking account of the current corticalDepth
corticalDepth = 0.5;
m.vertices = surface.inner.vtcs+corticalDepth*(surface.outer.vtcs-surface.inner.vtcs);
m.faceIndexList  = surface.inner.tris;


% get the FNIRT warp coeff files
if exist('FNIRT')==7
  pathname = 'FNIRT';
else
  pathname = pwd;
end
if fieldIsNotDefined(params,'func2anatFnirtCoeffs')
  params.func2anatFnirtCoeffs = [];
end
[filename,pathname2] = uigetfile([pathname '/*.hdr'], 'Choose Functional to Structural FNIRT coefficient file (Press Cancel if none is required)',params.func2anatFnirtCoeffs);
if ~isnumeric(filename)
  params.func2anatFnirtCoeffs = [pathname2 filename];
  pathname = pathname2;
end
if fieldIsNotDefined(params,'anat2funcFnirtCoeffs')
  params.anat2funcFnirtCoeffs = [];
end
[filename,pathname] = uigetfile([pathname '/*.hdr'], 'Choose Structural to Functional FNIRT coefficient file (Press Cancel if none is required)',params.anat2funcFnirtCoeffs);
if ~isnumeric(filename)
  params.anat2funcFnirtCoeffs = [pathname filename];
end


%%%%% All that's before this point should outside this function because not specific to a scanNum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scanNum = scanParams.scanNum;
if ~fieldIsNotDefined(params,'func2anatFnirtCoeffs') || ~fieldIsNotDefined(params,'anat2funcFnirtCoeffs')
  scanHdr = viewGet(thisView,'niftiHdr',scanNum);
end

%This assume that the base is loaded. if not, then the current base is chosen, which will be wrong if it's neither a flat map, a surface or the canonical base
base2scan = viewGet(thisView,'base2scan',scanNum,[],baseNum);
scan_surf = round(base2scan*[m.vertices,ones(length(m.vertices),1)]');


%%%% remove vertices/faces outside ROI
margin = ceil([params.radius params.radius params.radius]./viewGet(thisView,'scanVoxelSize'));
subsetBoxScan = getRoisBox(thisView,scanNum,margin,viewGet(thisView,'roinum',params.roiMask));

scanCoords = getROICoordinates(thisView,viewGet(thisView,'roinum',params.roiMask),scanNum);
scanCoords = [scanCoords; ones(1,size(scanCoords,2))];
scanDims = viewGet(thisView,'scanDims',scanNum);
% first apply non-linear transformation if needed
if ~fieldIsNotDefined(params,'func2anatFnirtCoeffs')
  scanCoords = applyFNIRT(scanCoords,params.func2anatFnirtCoeffs,scanDims,scanHdr,1);
end

%find the coordinates of the subset box including all the voxels from the ROI 
subsetBoxBase(1,1) = min(scanCoords(1,:))-margin(1);
subsetBoxBase(1,2) = max(scanCoords(1,:))+margin(1);
subsetBoxBase(2,1) = min(scanCoords(2,:))-margin(2);
subsetBoxBase(2,2) = max(scanCoords(2,:))+margin(2);
subsetBoxBase(3,1) = min(scanCoords(3,:))-margin(3);
subsetBoxBase(3,2) = max(scanCoords(3,:))+margin(3);


%find vertices that are in the box
verticesInBox = scan_surf(1,:)>=subsetBoxBase(1,1) & scan_surf(1,:)<=subsetBoxBase(1,2) &...
                scan_surf(2,:)>=subsetBoxBase(2,1) & scan_surf(2,:)<=subsetBoxBase(2,2) &...
                scan_surf(3,:)>=subsetBoxBase(3,1) & scan_surf(3,:)<=subsetBoxBase(3,2);
%remove vertices outside ROI
m.vertices = m.vertices(verticesInBox,:);
m.innerVertices = surface.inner.vtcs(verticesInBox,:);
m.outerVertices = surface.outer.vtcs(verticesInBox,:);
%replace ones by its 'vertex in the box' number
verticesInBox = double(verticesInBox);     %need to convert to double first     
verticesInBox(verticesInBox>0)=double(1:nnz(verticesInBox));
%remove any face that involves vertices outside the box
m.faceIndexList = verticesInBox(m.faceIndexList);
m.faceIndexList = m.faceIndexList(all(m.faceIndexList>0,2),:);

% calculate the connection matrix
[m.uniqueVertices,m.vertsToUnique,m.UniqueToVerts] = unique(m.vertices,'rows');
m.uniqueFaceIndexList = findUniqueFaceIndexList(m);
m.connectionMatrix = findConnectionMatrix(m);
scaleFactor =  baseHdr.pixdim(2:4)';
D = find3DNeighbourDists(m,scaleFactor);


m.uniqueInnerVertices = m.innerVertices(m.vertsToUnique,:);
m.uniqueOuterVertices = m.outerVertices(m.vertsToUnique,:);
corticalThickness = sqrt(sum((m.uniqueInnerVertices - m.uniqueOuterVertices).^2,2));
maxCorticalThickness = max(corticalThickness);

% min_surf = round(max(min(scan_surf(1:3,:)'),[1 1 1]));
% max_surf = round(min(max(scan_surf(1:3,:)'),d.dim(1:3)));
disp('Loading data')
%tmp = loadScan(thisView,scanNum,[],[min_surf(3) max_surf(3)],[],[min_surf(1) max_surf(1)],[min_surf(2) max_surf(2)]);
tmp = loadScan(thisView,scanNum,[],[subsetBoxScan(3,1) subsetBoxScan(3,2)],[],[subsetBoxScan(1,1) subsetBoxScan(1,2)],[subsetBoxScan(2,1) subsetBoxScan(2,2)]);
d.data = tmp.data;
clear('tmp');
d.dim = size(d.data);

col_mean = repmat(mean(d.data,4),[1 1 1, d.dim(4)]);
d.data=d.data-col_mean;
d.data=100*d.data./col_mean;
clear col_mean;

disp('preparing data for classification')
cc = viewGet(thisView,'concatinfo',scanNum);

run=[];
for i=1:size(cc.runTransition,1)
  run(cc.runTransition(i,1):cc.runTransition(i,2))=i;
  %remove any event that's too close to the end of the run
  for j=1:length(d.stimvol)
    eventsToRemove = find(d.stimvol{j}>d.dim(4)-scanParams.eventLength-scanParams.hdLag+1);
    if ~isempty(eventsToRemove)
      fprintf('(roiClassification) Removing %d event(s) because they''re too close to the end of the run\n',length(eventsToRemove));
      d.stimvol{j}(eventsToRemove) = [];
    end
  end
end

lab = [];
for i=1:length(d.stimvol)
    lab(d.stimvol{i})=i;
end


idx = find(lab>0);
if scanParams.averageEvent	 %average across the TRs for each event
  m_ = NaN([size(d.data,1) size(d.data,2) size(d.data,3) size(idx,2)]);
  for i=1:size(idx,2)
    m_(:,:,:,i)=mean(d.data(:,:,:,idx(i)+scanParams.hdLag:idx(i)+scanParams.eventLength+scanParams.hdLag-1),4);
  end
  d.data=m_;clear m_
  run=run(idx);
  lab = lab(idx);
elseif scanParams.eventLength>1 %create instance from each TR in the stim duration
  for i=1:size(idx,2)
  m_(:,:,:,idx(i):idx(i)+scanParams.eventLength-1)=d.data(:,:,:,idx(i)+scanParams.hdLag:idx(i)+scanParams.eventLength+scanParams.hdLag-1);
  l_(idx(i):idx(i)+scanParams.eventLength-1)=repmat(lab(idx(i)),1,scanParams.eventLength);
  r_(idx(i):idx(i)+scanParams.eventLength-1)=repmat(run(idx(i)),1,scanParams.eventLength);
  end
  d.data=m_(l_>0);
  lab = l_(l_>0);
  run=r_(l_>0);
  clear m_ l_ r_
else %create instance from individual TRs related to each TR
  d.data = d.data(:,:,:,idx+scanParams.hdLag);
  run=run(idx);
  lab = lab(idx);
end

baseROIcoords = base2scan \  scanCoords;
[roiNearestVtcs, distances] = assignToNearest(m.uniqueVertices, baseROIcoords(1:3,:)'); 
%remove any voxel that's not within a sensible distance of the middle of the cortical sheet
roiNearestVtcs = roiNearestVtcs(distances<=maxCorticalThickness/2);
scanCoords = scanCoords(:,distances<=maxCorticalThickness/2);
baseROIcoords = baseROIcoords(:,distances<=maxCorticalThickness/2);

offsets=make_sphere(params.radius);
d.data = permute(d.data,[4 1 2 3]);

% disp('(classify_corticalSearchlight) Computing Dijkstra distances ...') %%JB: this is too long is ROI is large
% dijkstraDistances = dijkstrap( D, roiNearestVtcs);  

debug = false;
if debug
  debugFigure = figure;
  patch('vertices', m.vertices, 'faces', m.faceIndexList,'FaceVertexCData', [1 1 1]*0.5,'facecolor','flat','edgecolor','none');
  hold on
  patch('vertices', m.outerVertices, 'faces', m.faceIndexList,'FaceVertexCData', [1 1 1]*0.5,'facecolor','flat','edgecolor','none','faceAlpha',.4);
  shading flat
  lighting phong
  axis equal vis3d
  % light from the back...
  m.light_handle(1) = lightangle(-90,45);
  m.light_handle(2) = lightangle(+90,45);
  hDebug = [];
end



disppercent(-inf,'(classify_corticalSearchlight) Classifying based on spotlight....');
corr = NaN(length(roiNearestVtcs),length(run));
size_light = NaN(1,length(roiNearestVtcs));
for i=1:length(roiNearestVtcs)
  
  %find all voxels in a spherical spotlight
  spotlightBaseCoords=repmat(baseROIcoords(1:3,i),1,size(offsets,2))+offsets;
  if debug && rem(i,20)==1
    figure(debugFigure);
    delete(hDebug);
    hDebug(1) = plot3(spotlightBaseCoords(1,:),spotlightBaseCoords(2,:),spotlightBaseCoords(3,:),'.b');
  end
  
  %find middle vertices that are within radius of the center of the spotlight
  dijkstraDistances = dijkstrap( D, roiNearestVtcs(i));
  spotlightMidVertices = m.uniqueVertices(dijkstraDistances<params.radius,:);
  spotlightInnerVertices = m.uniqueInnerVertices(dijkstraDistances<params.radius,:);
  spotlightOuterVertices = m.uniqueOuterVertices(dijkstraDistances<params.radius,:);
  spotlightThickness = corticalThickness(dijkstraDistances<params.radius);
  if debug && rem(i,20)==1
    figure(debugFigure);
    hDebug(2) = plot3(spotlightMidVertices(:,1),spotlightMidVertices(:,2),spotlightMidVertices(:,3),'.y');
  end
  
  %find middle vertices that are closest to all voxels in the spotlight
  [spotlightNearestVtcs, distances] = assignToNearest(spotlightMidVertices, spotlightBaseCoords'); 
  %first exclude any that are further than half of the max cortical thickness
  spotlightNearestVtcs = spotlightNearestVtcs(distances<=maxCorticalThickness/2);
  spotlightBaseCoords = spotlightBaseCoords(:,distances<=maxCorticalThickness/2);
  
  %spotlightNearestMidVtcsCoords = spotlightMidVertices(spotlightNearestVtcs,:);
  spotlightNearestInnerVtcsCoords = spotlightInnerVertices(spotlightNearestVtcs,:);
  spotlightNearestOuterVtcsCoords = spotlightOuterVertices(spotlightNearestVtcs,:);
  spotlightThickness = spotlightThickness(spotlightNearestVtcs,:);
  
  %now only keep spotlight voxels whose coordinates project within on the inner/outer segment of the closest
  % in other words if the dot product of inner/voxel segment is less than the distance between inner and outer, but positive
  dotProduct = sum((spotlightBaseCoords' - spotlightNearestInnerVtcsCoords).* (spotlightNearestOuterVtcsCoords - spotlightNearestInnerVtcsCoords),2);
  spotlightBaseCoords = spotlightBaseCoords(:,dotProduct<=spotlightThickness.^2 & dotProduct>=0);

  if debug && rem(i,20)==1
    figure(debugFigure);
    hDebug(3) = plot3(spotlightBaseCoords(1,:),spotlightBaseCoords(2,:),spotlightBaseCoords(3,:),'or');
    drawnow;
  end
  
  
  if ~isempty(spotlightBaseCoords)
    spotlightBoxCoords=unique(round(base2scan * [spotlightBaseCoords;ones(1,size(spotlightBaseCoords,2))])','rows');
    % apply non-linear transformation if needed
    if ~fieldIsNotDefined(params,'anat2funcFnirtCoeffs')
      spotlightBoxCoords = applyFNIRT(spotlightBoxCoords',params.anat2funcFnirtCoeffs,scanDims,scanHdr,0)';
    end
    spotlightBoxCoords=spotlightBoxCoords(:,1:3) - repmat(subsetBoxScan(:,1)',size(spotlightBoxCoords,1),1) +1;
    spotlightboxIndex = mrSub2ind(d.dim(1:3),spotlightBoxCoords(:,1),spotlightBoxCoords(:,2),spotlightBoxCoords(:,3));
    %remove any voxel outside the box (that might happen if the roi is at the edge of the scan)
    spotlightboxIndex(isnan(spotlightboxIndex))=[];
    if ~isempty(spotlightboxIndex)
      patt = d.data(:,spotlightboxIndex);
      for r=1:size(cc.runTransition,1)
        try
          corr(i,run==r) = lab(run==r)'==classify(patt(run==r,:),patt(run~=r,:),lab(run~=r),'diagLinear');
        catch
          corr(i,run==r)=nan;
        end
      end
    else
      patt=[];
    end
  else
    patt=[];
  end
  size_light(i)=size(patt,2);
  disppercent(i/length(roiNearestVtcs));

end
disppercent(inf); 


precision = mrGetPref('defaultPrecision');

%calculate global and class specific accuracies
acc=mean(corr,2);
acc_Class=nan(length(roiNearestVtcs),length(d.stimvol));
for i=1:length(d.stimvol)
    acc_Class(:,i)=mean(corr(:,lab==i),2);
end

if params.sigTest
  p = binomTest(sum(corr,2),length(lab),1/length(d.stimvol),'Greater');
  p_Class=nan(size(acc_Class));
  for i=1:length(d.stimvol)
       p_Class(:,i) =  binomTest(sum(corr(:,lab==i),2),sum(lab==i),1/length(d.stimvol),'Greater');
  end

  [p fdr fwe] = transformStatistic(p,precision,params);
  [p_Class fdr_Class fwe_Class] = transformStatistic(p_Class,precision,params);
end


%Now let's put the data back in the scan volume
scanCoordsIndex = sub2ind(scanDims,scanCoords(1,:)',scanCoords(2,:)',scanCoords(3,:)');

out.accDiag = NaN(scanDims,precision);
out.pDiag = NaN(scanDims,precision);
out.fdrDiag = NaN(scanDims,precision);
out.fweDiag = NaN(scanDims,precision);
out.accDiag_Class = NaN([length(d.stimvol) scanDims],precision);
out.pDiag_Class = NaN([length(d.stimvol) scanDims],precision);
out.fdrDiag_Class = NaN([length(d.stimvol) scanDims],precision);
out.fweDiag_Class = NaN([length(d.stimvol) scanDims],precision);

out.accDiag(scanCoordsIndex) = acc;
out.pDiag(scanCoordsIndex) = p;
out.fdrDiag(scanCoordsIndex) = fdr;
out.fweDiag(scanCoordsIndex) = fwe;
out.accDiag_Class(:,scanCoordsIndex) = acc_Class';
out.pDiag_Class(:,scanCoordsIndex) = p_Class';
out.fdrDiag_Class(:,scanCoordsIndex) = fdr_Class';
out.fweDiag_Class(:,scanCoordsIndex) = fwe_Class';

out.accDiag_Class = permute(out.accDiag_Class,[2 3 4 1]);
out.pDiag_Class = permute(out.pDiag_Class,[2 3 4 1]);
out.fdrDiag_Class = permute(out.fdrDiag_Class,[2 3 4 1]);
out.fweDiag_Class = permute(out.fweDiag_Class,[2 3 4 1]);


if ~fieldIsNotDefined(params,'anat2funcFnirtCoeffs')
  out.accDiag = applyFSLWarp(out.accDiag, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.pDiag = applyFSLWarp(out.pDiag, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.fdrDiag = applyFSLWarp(out.fdrDiag, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.fweDiag = applyFSLWarp(out.fweDiag, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.accDiag_Class = applyFSLWarp(out.accDiag_Class, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.pDiag_Class = applyFSLWarp(out.pDiag_Class, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.fdrDiag_Class = applyFSLWarp(out.fdrDiag_Class, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
  out.fweDiag_Class = applyFSLWarp(out.fweDiag_Class, params.anat2funcFnirtCoeffs, 'tempInput.img', scanHdr, 'nearest', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sub functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scanCoords = applyFNIRT(scanCoords,fnirtCoeffsFile,scanDims,scanHdr,verbose)
outsideVoxels = any(scanCoords(1:3,:)< ones(3,size(scanCoords,2)) |...
scanCoords(1:3,:)> repmat(scanDims',1,size(scanCoords,2)));
if any(outsideVoxels) 
  if verbose
    mrWarnDlg('(applyFNIRT) Some voxels in are outside the scan volume, removing...');
  end
  scanCoords(:,outsideVoxels)=[];
end
warpIndices = sub2ind(scanDims(1:3),scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
roiVolume = zeros(scanDims(1:3));
roiVolume(warpIndices) = 1;
warpedRoiVolume = applyFSLWarp(roiVolume, fnirtCoeffsFile, 'tempInput.img', scanHdr, 'nearest',verbose);
warpedCoordsIndices = find(warpedRoiVolume==1);
[x,y,z] = ind2sub(scanDims(1:3),warpedCoordsIndices);
scanCoords = [x y z ones(length(warpedCoordsIndices),1)]';


