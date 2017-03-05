% remapSurfaces.m
%
%       $Id: remapSurfaces.m 1833 2010-11-13 18:37:14Z julien $	
%      usage: remapSurfaces(fssubject,fsaverage)
%         by: julien besle
%       date: 15/05/2014
%    purpose: remaps surface vertices of one Freesurfer subject to another 
%             (usually fsaverage, but not necessarily) using the spherical
%             registration output of recon-all (subj/surf/?h.sphere.reg)%        
%     output: left and right GM and WM surfaces with mesh of one subject
%             and coordinates of the other, and vice-versa

function remapSurface(fssubject,fsaverage)

if isunix || ismac
  if isempty(getenv('SUBJECTS_DIR'))
    mrErrorDlg('(remapSurfaces) FreeSurfer environment variable SUBJECTS_DIR is not set');
  % implement another way to get the subjects folder
  end
  subjPath = [getenv('SUBJECTS_DIR') '/' fssubject];
  if isempty(dir(subjPath))
    mrWarnDlg(['(remapSurfaces) Freesurfer subject ' fssubject ' does not exist']);
    return
  end
  fsaveragePath = [getenv('SUBJECTS_DIR') '/' fsaverage];
  if isempty(dir(fsaveragePath))
    mrWarnDlg(['(remapSurfaces) Freesurfer subject ' fsaverage ' does not exist']);
    return
  end
else
  mrErrorDlg('(remapSurfaces) Not implemented for platforms other than Unix or Mac');
  % implement another way to get the subjects folder
end

if isempty(dir([subjPath '/surfRelax']))
  mrWarnDlg(['(remapSurfaces) surfRelax folder does not exist in ' subjPath '. You must first run mlrImportFreesurfer.']);
  return
end

if isempty(dir([fsaveragePath '/surfRelax']))
  mrWarnDlg(['(remapSurfaces) surfRelax folder does not exist in ' fsaveragePath '. You must first run mlrImportFreesurfer.']);
  return
end
        
side = {'left','right'};
surfs = {'GM','WM'};
fsSide = {'lh','rh'};

disp('(remapSurfaces) Will process:');
for iSide=1:2
  for iSurf = 1:2
    %find all OFF files for a given side and surface (WM or GM)
    subjFiles{iSide,iSurf} = dir([subjPath '/surfRelax/' fssubject '_' side{iSide} '_' surfs{iSurf} '*.off']);
    %find OFF files to exclude (already processed or flat maps)
    toExclude = dir([subjPath '/surfRelax/' fssubject '_' side{iSide} '_' surfs{iSurf} '*Colin*.off']);
    toExclude = [toExclude; dir([subjPath '/surfRelax/' fssubject '_' side{iSide} '_' surfs{iSurf} '*MNI*.off'])];
    toExclude = [toExclude; dir([subjPath '/surfRelax/' fssubject '_' side{iSide} '_' surfs{iSurf} '*Flat*.off'])];
    
    [~,toKeep]=setdiff({subjFiles{iSide,iSurf}(:).name},{toExclude(:).name});
    subjFiles{iSide,iSurf} = {subjFiles{iSide,iSurf}(toKeep).name};
    disp(subjFiles{iSide,iSurf}');
%     if length(subjFiles{iSide,iSurf})>2
%       dir([subjPath '/surfRelax/' fssubject '_' side{iSide} '_' surfs{iSurf} '*.off'])
%       filelist = input('Which files do you wish to transform (input index vector) ?')
%       subjFiles{iSide,iSurf} = subjFiles{iSide,iSurf}(filelist);
%     end
  end
  disp([subjPath '/surfRelax/' fssubject '_' side{iSide} '_Curv.vff'])
end

for iSide=1:2
  %get reg sphere surfaces
  [vertSphereSubj, triSphereSubj] = freesurfer_read_surf([subjPath '/surf/' fsSide{iSide} '.sphere.reg']);
  [vertSphereAverage, triSphereAverage] = freesurfer_read_surf([fsaveragePath '/surf/' fsSide{iSide} '.sphere.reg']);
  
  % compute re-gridding matrix (expresses the vertices coordinates in one
  % sphere as a linear combination of face vertices in the other)
  [averageToSubj,subjToAverage] = findCorrespondingVertices(vertSphereSubj,vertSphereAverage,triSphereSubj,triSphereAverage);

  for iSurf = 1:2
    %get surfaces in OFF format
    averageSurf = loadSurfOFF([fsaveragePath '/surfRelax/' fsaverage '_' side{iSide} '_' surfs{iSurf} '.off']);
    
    for jSurf = subjFiles{iSide,iSurf}
      pattern = ['_' side{iSide} '_' surfs{iSurf}];
      fssubjectPrefix = jSurf{1}([1:strfind(jSurf{1},pattern)-1 strfind(jSurf{1},pattern)+length(pattern):end-4]);
      subjSurf = loadSurfOFF([subjPath '/surfRelax/' jSurf{1}]);

      %apply regridding matrix to surfaces
      thisAverageSurf = averageSurf;
      tmpVtcs = averageToSubj*thisAverageSurf.vtcs;
      thisAverageSurf.vtcs = subjToAverage*subjSurf.vtcs;
      subjSurf.vtcs = tmpVtcs;

      %change file name
      [path,filename,extension]=fileparts(subjSurf.filename);
      subjSurf.filename = [path '/' filename '_' fsaverage extension];
      [path,filename,extension]=fileparts(thisAverageSurf.filename);
      thisAverageSurf.filename = [path '/' filename '_' fssubjectPrefix extension];

      %save file
      disp(['(remapSurfaces) Writing ' subjSurf.filename]);
      writeOFF(subjSurf, subjSurf.filename);
      disp(['(remapSurfaces) Writing ' thisAverageSurf.filename]);
      writeOFF(thisAverageSurf, thisAverageSurf.filename);

    end
  end
  
  %interpolate curvature data
  subjCurv = loadVFF([subjPath '/surfRelax/' fssubject '_' side{iSide} '_Curv.vff'])';
  subjToAverageCurv = subjToAverage*subjCurv;
  subjToAverageCurv(subjToAverageCurv>max(subjCurv))=max(subjCurv); %clip the data to original min/max (because saveVFF normalizes them)
  subjToAverageCurv(subjToAverageCurv<min(subjCurv))=min(subjCurv);
  saveVFF([fsaveragePath '/surfRelax/' fsaverage '_' side{iSide} '_Curv_' fssubject '.vff'], subjToAverageCurv');
  averageCurv = loadVFF([fsaveragePath '/surfRelax/' fsaverage '_' side{iSide} '_Curv.vff'])';
  averageToSubjCurv = averageToSubj*averageCurv;
  averageToSubjCurv(averageToSubjCurv>max(averageCurv))=max(subjCurv);
  averageToSubjCurv(averageToSubjCurv<min(averageCurv))=min(averageCurv);
  saveVFF([subjPath '/surfRelax/' fssubject '_' side{iSide} '_Curv_' fsaverage '.vff'], averageToSubjCurv');

end
    
    

% function [vert2to1,vert1to2]=findCorrespondingVertices(vert1,vert2,tri1,tri2)
function [vert2to1,vert1to2]=findCorrespondingVertices(vert1,vert2,tri1,tri2)

nVerts1=size(vert1,1);
nVerts2=size(vert2,1);


% % %brute-force method (takes 2-3 hours for just closest vertex): 
% % %loop over each vertex of a mesh, compute its distance
% % %to entire other mesh and find closest vertex (shortest distance)
% % vert2to1=nan(nVerts1,1);
% % for iVert=1:nVerts1
% %   [value,vert2to1(iVert)]=min(sum((repmat(vert1(iVert,:),nVerts2,1)-vert2).^2,2));
% % end
% % 
% % vert1to2=nan(nVerts2,1);
% % hWaitbar = mrWaitBar(0,'(remapSurfaces) Finding corresponding vertices in both meshes...');
% % for iVert=1:nVerts2
% %   [value,vert1to2(iVert)]=min(sum((repmat(vert2(iVert,:),nVerts1,1)-vert1).^2,2));
% %   mrWaitBar(iVert/nVerts2,hWaitbar)
% % end
% % mrCloseDlg(hWaitbar);

% % %Dave Langers' method (takes about 3 hour for one remapping)
% % nTris2=size(tri2,1);
% % vert2to1 = zeros(nVerts1, 3);
% % weights2to1d = zeros(nVerts1, 3);
% % score = inf(nVerts1, 1);
% % hWaitbar = mrWaitBar(0,'(remapSurfaces) Finding corresponding vertices in both meshes (Dave)...');
% % cFace=0;
% % for fitindex = tri2'
% %   cFace=cFace+1;
% %   fitparam = vert1/vert2(fitindex, :);
% %   fitscore = sum(abs(fitparam), 2);
% %   vertset = find((fitscore < score) & (sum(fitparam, 2) > 0));
% %   vert2to1(vertset, :) = repmat(fitindex',   length(vertset), 1);
% %   weights2to1d(vertset, :) = fitparam(vertset, :)./repmat(sum(fitparam(vertset, :), 2), 1, 3);
% %   score(vertset) = fitscore(vertset);
% %   mrWaitBar(cFace/nTris2,hWaitbar)
% % end
% % mrCloseDlg(hWaitbar);

%%%%% FASTER method taking into account spatial proximity of corresponding vertices in both meshes
%%%%%  (takes about 4 min for both remappings)
hWaitbar = mrWaitBar(0,'(remapSurfaces) Finding corresponding vertices in both meshes...');
%vectors that will index one set of vertices into the other
vert2to1=nan(nVerts1,3);
vert1to2=nan(nVerts2,3);
%associated weights
weights2to1=nan(nVerts1,3);
weights1to2=nan(nVerts2,3);

%calculate azimuth and elevation of both meshes
[az1,el1]=cart2sph(vert1(:,1),vert1(:,2),vert1(:,3));
[az2,el2]=cart2sph(vert2(:,1),vert2(:,2),vert2(:,3));

%average azimuth and elevation of triangles
aztri1=(az1(tri1(:,1))+az1(tri1(:,2))+az1(tri1(:,3)))/3;
eltri1=(el1(tri1(:,1))+el1(tri1(:,2))+el1(tri1(:,3)))/3;
aztri2=(az2(tri2(:,1))+az2(tri2(:,2))+az2(tri2(:,3)))/3;
eltri2=(el2(tri2(:,1))+el2(tri2(:,2))+el2(tri2(:,3)))/3;

%average location of triangles
meanTri1=(vert1(tri1(:,1),:)+vert1(tri1(:,2),:)+vert1(tri1(:,3),:))/3;
meanTri2=(vert2(tri2(:,1),:)+vert2(tri2(:,2),:)+vert2(tri2(:,3),:))/3;

%now go in patches of approximately equal area and increasing elevation,
margin = 0.02;
nSteps=20;
elevations = linspace(-pi/2,pi/2,nSteps);
counter=0;
totalCount = 2*sum(1:nSteps/2)-nSteps/2; %not sure that works if nSteps is odd
for iEl = 1:length(elevations)-1
  elMin= elevations(iEl);
  elMax= elevations(iEl+1);
  elMinMargin= elevations(iEl)-margin;
  elMaxMargin= elevations(iEl+1)+margin;
  el1Indices = el1>=elMin  & el1<=elMax;
  el2Indices = el2>=elMin  & el2<=elMax;
  switch(iEl)
    case 1
% %       el1IndicesMargin = (el1>=elMinMargin | (el1-pi)>=elMinMargin) & el1<=elMaxMargin;
% %       el2IndicesMargin = (el2>=elMinMargin | (el2-pi)>=elMinMargin) & el2<=elMaxMargin;
      eltri1Indices = (eltri1>=elMinMargin | (eltri1-pi)>=elMinMargin)  & eltri1<=elMaxMargin;
      eltri2Indices = (eltri2>=elMinMargin | (eltri2-pi)>=elMinMargin)  & eltri2<=elMaxMargin;
    case length(elevations)-1
% %       el1IndicesMargin = el1>=elMinMargin  & ((el1<=elMaxMargin) | (el1+pi)<=elMaxMargin);
% %       el2IndicesMargin = el2>=elMinMargin  & ((el2<=elMaxMargin) | (el2+pi)<=elMaxMargin);
      eltri1Indices = eltri1>=elMinMargin  & (eltri1<=elMaxMargin | (eltri1+pi)<=elMaxMargin);
      eltri2Indices = eltri2>=elMinMargin  & (eltri2<=elMaxMargin | (eltri2+pi)<=elMaxMargin);
    otherwise
% %       el1IndicesMargin = el1>=elMinMargin & el1<=elMaxMargin;
% %       el2IndicesMargin = el2>=elMinMargin & el2<=elMaxMargin;
      eltri1Indices = eltri1>=elMinMargin  & eltri1<=elMaxMargin;
      eltri2Indices = eltri2>=elMinMargin  & eltri2<=elMaxMargin;
  end
  %rotating around the azimuth axis
  azimuths = linspace(-pi,pi,nSteps/2-abs(nSteps/2-iEl)+1); %change the step size according to elevation
  for iAz = 1:length(azimuths)-1
    counter=counter+1;
    azMin= azimuths(iAz);
    azMax= azimuths(iAz+1);
    azMinMargin= azimuths(iAz)-margin;
    azMaxMargin= azimuths(iAz+1)+margin;
    az1Indices = az1>=azMin & az1<=azMax;
    az2Indices = az2>=azMin & az2<=azMax;
    switch(iAz)
      case [1,length(azimuths)-1]
% %         az1IndicesMargin = (az1>=azMinMargin | (az1-2*pi)>=azMinMargin) & ((az1<=azMaxMargin) | (az1+2*pi)<=azMaxMargin);
% %         az2IndicesMargin = (az2>=azMinMargin | (az2-2*pi)>=azMinMargin) & ((az2<=azMaxMargin) | (az2+2*pi)<=azMaxMargin);
        aztri1Indices = (aztri1>=azMin | (aztri1-2*pi)>=azMinMargin) & (aztri1<=azMax | (aztri1+2*pi)<=azMaxMargin);
        aztri2Indices = (aztri2>=azMin | (aztri2-2*pi)>=azMinMargin) & (aztri2<=azMax | (aztri2+2*pi)<=azMaxMargin);
      otherwise
% %         az1IndicesMargin = az1>=azMinMargin & az1<=azMaxMargin;
% %         az2IndicesMargin = az2>=azMinMargin & az2<=azMaxMargin;
        aztri1Indices = aztri1>=azMin & aztri1<=azMax;
        aztri2Indices = aztri2>=azMin & aztri2<=azMax;
    end
    
    %select the patch vertices and triangles by intersection of the elevation and azimuth indices
    elaz1Indices=find(el1Indices & az1Indices); %use find here because we'll need the actual indices later
    thisVert1 = vert1(elaz1Indices,:);
    elaztri1Indices=find(eltri1Indices & aztri1Indices);
    thisTri1=tri1(elaztri1Indices,:);
    elaz2Indices=find(el2Indices & az2Indices); 
    thisVert2 = vert2(elaz2Indices,:);
    elaztri2Indices=find(eltri2Indices & aztri2Indices);
    thisTri2=tri2(elaztri2Indices,:);

    nThisTris1=size(thisTri1,1);
    nThisTris2=size(thisTri2,1);
    nThisVerts1=size(thisVert1,1);
    nThisVerts2=size(thisVert2,1);

% %     %now loop over the subsets of triangles and for each, find which
% %     %vertices in the other mesh they can best model (as a linear combination of their vertices)
% %     %(this is adapted from Dave Langers' code)
% % % hWaitbar = mrWaitBar(0,'(remapSurfaces) Finding corresponding vertices in both meshes...');
% %     thisWeights2to1=nan(nThisVerts1,3);
% %     thisWeights1to2=nan(nThisVerts2,3);
% %     thisVert2to1=nan(nThisVerts1,3);
% %     thisVert1to2=nan(nThisVerts2,3);
% %     score = inf(nThisVerts2, 1);
% %     for iTri=1:nThisTris1
% %       fitparam = thisVert2/vert1(thisTri1(iTri,:)', :);
% %       fitscore = sum(abs(fitparam), 2);
% %       vertset = find((fitscore < score) & (sum(fitparam, 2) > 0));
% %       thisVert1to2(vertset, :) = repmat(thisTri1(iTri,:),   length(vertset), 1);
% %       thisWeights1to2(vertset, :) = fitparam(vertset, :)./repmat(sum(fitparam(vertset, :), 2), 1, 3);
% %       score(vertset) = fitscore(vertset);
% % % mrWaitBar(iTri/nThisTris1,hWaitbar)
% %     end
% % % mrCloseDlg(hWaitbar);
% %     vert1to2(elaz2Indices,:)=thisVert1to2;
% %     weights1to2(elaz2Indices,:)=thisWeights1to2;
% %     
% %     score = inf(nThisVerts1, 1);
% %     for iTri=1:nThisTris2
% %       fitparam = thisVert1/vert2(thisTri2(iTri,:)', :);
% %       fitscore = sum(abs(fitparam), 2);
% %       vertset = find((fitscore < score) & (sum(fitparam, 2) > 0));
% %       thisVert2to1(vertset, :) = repmat(thisTri2(iTri,:),   length(vertset), 1);
% %       thisWeights2to1(vertset, :) = fitparam(vertset, :)./repmat(sum(fitparam(vertset, :), 2), 1, 3);
% %       score(vertset) = fitscore(vertset);
% %     end
% %     vert2to1(elaz1Indices,:)=thisVert2to1;
% %     weights2to1(elaz1Indices,:)=thisWeights2to1;
    
%%%% Slightly faster method (about 2 min): find the triangle in one mesh whose centre
%%%% (average position) is closest to a given vertex in the other and then
%%%% solve the position equation to get the weights
    thisMeanTri1=meanTri1(elaztri1Indices,:);
    thisMeanTri2=meanTri2(elaztri2Indices,:);
    triIndex=nan(nThisVerts2,1);
    weights=nan(nThisVerts2,3);
    for iVert=1:nThisVerts2
      %find the triangle closest to a given vertex
      [~,triIndex(iVert)]=min(sum((repmat(thisVert2(iVert,:),nThisTris1,1)-thisMeanTri1).^2,2));
      %find the linear combination of the triangle vertices that gives the given vertex
      weights(iVert,:) = thisVert2(iVert,:)/vert1(thisTri1(triIndex(iVert)',:)',:);
%       Y = thisVert2(iVert,:)';
%       X = vert1(thisTri1(triIndex(iVert)',:)',:)';
%       weights(iVert,:) =  ((X'*X)\X'*Y)';
% mrWaitBar(iVert/nThisVerts2,hWaitbar)
    end
% mrCloseDlg(hWaitbar);
    vert1to2(elaz2Indices,:)=thisTri1(triIndex,:);
    weights1to2(elaz2Indices,:)=weights;
  
    triIndex=nan(nThisVerts1,1);
    weights=nan(nThisVerts1,3);
    for iVert=1:nThisVerts1
      [~,triIndex(iVert)]=min(sum((repmat(thisVert1(iVert,:),nThisTris2,1)-thisMeanTri2).^2,2));
      weights(iVert,:) = thisVert1(iVert,:)/vert2(thisTri2(triIndex(iVert)',:)',:);
%       Y = thisVert1(iVert,:)';
%       X = vert2(thisTri2(triIndex(iVert)',:)',:)';
%       weights(iVert,:) =  ((X'*X)\X'*Y)';
    end
    vert2to1(elaz1Indices,:)=thisTri2(triIndex,:);
    weights2to1(elaz1Indices,:)=weights;

    
% % %older version not using the triangles
%%%%% (2 min for single closest vertex and 6 min for several closest)
% %     elaz1IndicesMargin=find(el1IndicesMargin & az1IndicesMargin);
% %     thisVert1Margin = vert1(elaz1IndicesMargin,:);
% %     elaz2IndicesMargin=find(el2IndicesMargin & az2IndicesMargin);
% %     thisVert2Margin = vert2(elaz2IndicesMargin,:);

% %     %now loop over the subsets of vertices, compute distances and select
% %     %closest  corresponding vertex 
% %     nThisVerts1=size(thisVert1,1);
% %     nThisVerts2=size(thisVert2,1);
% %     nThisVerts1Margin=size(thisVert1Margin,1);
% %     nThisVerts2Margin=size(thisVert2Margin,1);
% %     thisVert2to1=nan(nThisVerts1,3);
% %     for iVert=1:nThisVerts1
% % %       %get closet vertex
% % %       [value,thisVert2to1(iVert)]=min(sum((repmat(thisVert1(iVert,:),nThisVerts2Margin,1)-thisVert2Margin).^2,2));
% %       %get first 3 closest vertices
% %       [~,sortingIndices]=sort(sum((repmat(thisVert1(iVert,:),nThisVerts2Margin,1)-thisVert2Margin).^2,2));
% %       thisVert2to1(iVert,:)=sortingIndices(1:3)';
% %     end
% %     thisVert1to2=nan(nThisVerts2,3);
% %     for iVert=1:nThisVerts2
% % %       %get closet vertex
% % %       [value,thisVert1to2(iVert)]=min(sum((repmat(thisVert2(iVert,:),nThisVerts1Margin,1)-thisVert1Margin).^2,2));
% %       %get first 3 closest vertices
% %       [~,sortingIndices]=sort(sum((repmat(thisVert2(iVert,:),nThisVerts1Margin,1)-thisVert1Margin).^2,2));
% %       thisVert1to2(iVert,:)=sortingIndices(1:3)';
% %     end
% %     
% %     vert2to1(elaz1Indices,:)=elaz2IndicesMargin(thisVert2to1);
% %     vert1to2(elaz2Indices,:)=elaz1IndicesMargin(thisVert1to2);

    mrWaitBar(counter/totalCount,hWaitbar)

  end
  
end
mrCloseDlg(hWaitbar);

%put indices and weights in a big sparse matrix so we can apply it in one multiplication
vert1to2=sparse(repmat((1:nVerts2)', 1, 3), vert1to2, weights1to2, nVerts2, nVerts1, nVerts2*3);
vert2to1=sparse(repmat((1:nVerts1)', 1, 3), vert2to1, weights2to1, nVerts1, nVerts2, nVerts1*3);


% % function [outvert1,outvert2]=computeRemappedCoords(vert1,vert2,vert2to1,vert1to2,sphereVert1,sphereVert2)
% % 
% % outvert1=computeRemappedCoords2(vert2,vert2to1,sphereVert1,sphereVert2);
% % outvert2=computeRemappedCoords2(vert1,vert1to2,sphereVert2,sphereVert1);
% % 
% % function outvert1=computeRemappedCoords2(vert2,vert2to1,sphereVert1,sphereVert2)
% % 
% % %nearestneighbour
% % nnoutvert1 = vert2(vert2to1(:,1),:);
% % 
% % %interpolation:
% % %express each vertex of one sphere as a weighted sum of the 3 closest
% % %points on the other sphere and solve the linear equation 
% % %(this fails if any two points are close to each other)
% % outvert1=nan(size(sphereVert1));
% % % warning('off','MATLAB:nearlySingularMatrix')
% % % warning('off','MATLAB:singularMatrix')
% % for i=1:size(sphereVert1,1)
% % %   if sqrt(sum(( sphereVert1(i,:)-sphereVert2(vert2to1(i,1),:)).^2,2))>0.04
% %     Y = sphereVert1(i,:)';
% %     X = sphereVert2(vert2to1(i,1:3)',:)';
% %     weights =  ((X'*X)\X'*Y)';  %could have done this step in findCorrespondingVertices to avoid repeating it 
% %     %for each surface (GM and WM), but doesn`t take that long, so couldn`t be bothered
% % 
% %     % apply the weights to the surface:
% %     outvert1(i,:) = weights*vert2(vert2to1(i,1:3),:);
% % %   else
% % %     outvert1(i,:)=nnoutvert1(i,:);
% % %   end
% % end
% % warning('on','MATLAB:nearlySingularMatrix')
% % warning('on','MATLAB:singularMatrix')
% % 
% % %If the resulting vertices are far from the nearest neighbor approximation
% % %(>1 mm seems to be a good cutoff), this is probably because the sphere vertex was too
% % %close to one of the 3 other ones to start with. In this case, use nearest neigbour
% % cutoff=1.7;
% % toofar=sqrt(sum((outvert1-nnoutvert1).^2,2))>cutoff;
% % toto= sqrt(sum((outvert1-nnoutvert1).^2,2));
% % figure;hist(toto,100);
% % fprintf('(remapSurfaces) Using nearest neighbour interpolation for %d/%d vertices.\n',nnz(toofar),size(sphereVert1,1));
% % outvert1(toofar,:)= nnoutvert1(toofar,:);


