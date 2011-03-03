function m_d=classify_corticalSearchlight(view,d,radius,params)
%This fucntion does a searchlight classification with a defined radius in
%an area defined by an ROI, which is normally a mask of the skullstripped
%brain. Has to be the run from the root directory of an MLR session.

if 0
    %Define which scan from the concatenation group is to be analyzed
    scan=1;


    view=newView;
    view=viewSet(view,'curGroup','Concatenation');
    load('Concatenation/erAnal/erAnal.mat');
    d = erAnal.d{scan};
    if isempty(d)
      disp('No analysis');
      return
    end

    %radius defined if not already
    if ~exist('radius')
        radius=input('Define radius');
    end
end
% a function to calculate the offsets of the spotlight voxels relative to
% the centre of the sphere



baseHdr = cbiReadNiftiHeader(['Segmentation/',uigetfile('Segmentation/*.hdr')]);
% load the appropriate surface files
% disp(sprintf('Loading %s', baseCoordMap.innerCoordsFileName));
surface.inner = loadSurfOFF(fullfile('Segmentation/', uigetfile('Segmentation/*WM.off')));
if isempty(surface.inner)
mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName)));
return
end
surface.inner = xformSurfaceWorld2Array(surface.inner, baseHdr);
% disp(sprintf('Loading %s', baseCoordMap.outerCoordsFileName));
surface.outer = loadSurfOFF(fullfile('Segmentation/', uigetfile('Segmentation/*GM.off')));
if isempty(surface.outer)
mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName)));
return
end
surface.outer = xformSurfaceWorld2Array(surface.outer, baseHdr);

corticalDepth = 0.5;
%% build up a mrMesh-style structure, taking account of the current corticalDepth
m.vertices = surface.inner.vtcs+corticalDepth*(surface.outer.vtcs-surface.inner.vtcs);
m.faceIndexList  = surface.inner.tris;
%% calculate the connection matrix
[m.uniqueVertices,m.vertsToUnique,m.UniqueToVerts] = unique(m.vertices,'rows');
m.uniqueFaceIndexList = findUniqueFaceIndexList(m);
m.connectionMatrix = findConnectionMatrix(m);
scaleFactor =  baseHdr.pixdim(2:4)';
D = find3DNeighbourDists(m,scaleFactor);

[~,basename,~] = fileparts(baseHdr.hdr_name);
base2scan = viewGet(view,'base2scan',params.scanNum,[],viewGet(view,'basenum',basename));

scan_surf = base2scan*[m.vertices,ones(length(m.vertices),1)]';

min_surf = round(max(min(scan_surf(1:3,:)'),[1 1 1]));
max_surf = round(min(max(scan_surf(1:3,:)'),d.dim(1:3)));
disp('Loading data')
tic
tmp = loadScan(view,params.scanNum,[],[min_surf(3) max_surf(3)],[],[min_surf(1) max_surf(1)],[min_surf(2) max_surf(2)]);

col_mean = repmat(mean(tmp.data,4),[1 1 1, d.dim(4)]);

d.data=tmp.data-col_mean;
d.data=100*d.data./col_mean;
clear tmp;


disp('preparing data for classification')
cc = viewGet(view,'concatinfo',params.scanNum);
lab = [];
for i=1:length(d.stimvol)
    lab(d.stimvol{i})=i;
end

run=[];
for i=1:size(cc.runTransition,1)
run(cc.runTransition(i,1):cc.runTransition(i,2))=i;
end


% d.tmean=mean(d.data,4);
% col_mean = repmat(d.tmean,[1 1 1 size(d.data,4)]);
% d.data = 100*(d.data-col_mean)./col_mean;


idx = find(lab>0);
if params.averageEvent %average across the TRs for each event
for i=1:size(idx,2)
m_(:,:,:,i)=mean(d.data(:,:,:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1),4);
end
d.data=m_;clear m_
run=run(idx);
lab = lab(idx);
elseif params.eventLength>1 %create instance from each TR in the stim duration
for i=1:size(idx,2)
m_(:,:,:,idx(i):idx(i)+params.eventLength-1)=d.data(:,:,:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1);
l_(idx(i):idx(i)+params.eventLength-1)=repmat(lab(idx(i)),1,params.eventLength);
r_(idx(i):idx(i)+params.eventLength-1)=repmat(run(idx(i)),1,params.eventLength);
end
d.data=m_(l_>0);
lab = l_(l_>0);
run=r_(l_>0);
clear m_ l_ r_
else %create instance from individual TRs related to each TR
d.data = d.data(:,:,:,idx+params.hdLag);
run=run(idx);
lab = lab(idx);
end

n.data=nan([d.dim(1:3) size(d.data,4)]);
n.data(min_surf(1):max_surf(1),min_surf(2):max_surf(2),min_surf(3):max_surf(3),:)=d.data;
d.data=n.data;
clear n


disp('Classifying')
for i=1:length(m.vertices)
  dist = dijkstrap( D, m.UniqueToVerts(i));dist=dist(m.UniqueToVerts);  
base_light=[m.vertices(dist<radius,:)';ones(1,size(m.vertices(dist<radius,:),1))];
scan_light=unique(round(base2scan * base_light)','rows')';
[~,j] = find((scan_light(1,:)<=0 | scan_light(1,:)>d.dim(1)) | (scan_light(2,:)<=0 | scan_light(2,:)>d.dim(2)) | (scan_light(3,:)<=0 | scan_light(3,:)>d.dim(3)));
scan_light=scan_light(:,setdiff(1:size(scan_light,2),j));
if ~isempty(scan_light)
    for j=1:size(scan_light,2)
        patt(j,:)=d.data(scan_light(1,j),scan_light(2,j),scan_light(3,j),:);
    end
    for r=1:size(cc.runTransition,1)
        try
        acc(r)=mean(lab(run==r)'==classify(patt(:,run==r)',patt(:,run~=r)',lab(run~=r),'diagLinear'));
        catch
            acc(r)=nan;
%        s_acc(r) = mean(lab(run==r)'==classifyWithSvm_2(patt(:,run==r),patt(:,run~=r),lab(run~=r)));
        end
    end
    mean_acc(i)=mean(acc);
%     mean_s_acc(i)=mean(s_acc);
else
    mean_acc(i)=NaN;
    patt=[];
%     mean_s_acc(i)=NaN;
end
size_light(i)=size(patt,1);
patt=[];
end,toc

c = loadVFF(fullfile('Segmentation/', uigetfile('Segmentation/*vff')));%('Segmentation/ab_left_curv.vff')
cmap = gray;
limval = 1.2;
c(c>limval) = limval;
c(c<-limval) = -limval;
c = round((size(cmap,1)-1)*(c-min(c))./(max(c)-min(c))+1);
combinedOverlay(:,1) = cmap(c);
combinedOverlay(:,2) = cmap(c);
combinedOverlay(:,3) = cmap(c);
combinedOverlay(~isnan(mean_acc),:)=squeeze(cdata2rgb(mean_acc(~isnan(mean_acc)),'hot',[0 1]));

figure,patch('vertices', m.vertices, 'faces', m.faceIndexList,'FaceVertexCData',combinedOverlay,'facecolor','interp','edgecolor','none');
colormap('hot')
colorbar
axis equal
axis vis3d
axis off

keyboard
% offsets=make_sphere(radius);
% roi_name='mask';
% 
% % %function that does the classification. 
% % classify_hem(view,scan,d,offsets,radius,roi)
% % 
% % function classify_hem(view,scan,d,offsets,radius,hem)
% 
% %laods
% load(['ROIs/',roi_name],roi_name);
% roi{1}=eval(roi_name);
% 
% roi{1}.scanCoords = getROICoordinates(view,roi{1},params.scanNum);
% 
% d.roi{1} = loadROITSeries(view,roi{1},params.scanNum);





warning('off');
% idx=cell(1,length(d.roi{1}.scanCoords));
% 
% % finds out which voxels are in which spotlight, and removes any that are
% % outside the head
% disppercent(-inf,'Indexing spotlight....');
% for i_vox=1:length(d.roi{1}.scanCoords)
%  
%     idx{i_vox}=repmat(d.roi{1}.scanCoords(:,i_vox),1,size(offsets,2))+offsets;
%     [~,j]=find(idx{i_vox}(1,:)<min(d.roi{1}.scanCoords(1,:)) | idx{i_vox}(2,:)<min(d.roi{1}.scanCoords(2,:)) | idx{i_vox}(3,:)<min(d.roi{1}.scanCoords(3,:)) );
%     idx{i_vox}=idx{i_vox}(:,setdiff(1:size(idx{i_vox},2),j));
%     [~,j]=find(idx{i_vox}(1,:)>max(d.roi{1}.scanCoords(1,:)) | idx{i_vox}(2,:)>max(d.roi{1}.scanCoords(2,:)) | idx{i_vox}(3,:)>max(d.roi{1}.scanCoords(3,:)) );
%     idx{i_vox}=idx{i_vox}(:,setdiff(1:size(idx{i_vox},2),j));
% 
% 
%     disppercent(i_vox/length(d.roi{1}.scanCoords));
% end
% 
% 
% disppercent(inf);



% preprocess the functional data



% in=[];
% dir=[];
% disppercent(-inf,'Creating instances...');
% indx=cell(1,length(d.stimvol));
% indx_start=cell(1,length(d.stimvol));
% indx_end=cell(1,length(d.stimvol));
% instances=cell(1,length(d.stimvol));
% 
% for i_stim=1:length(d.stimvol)
% indx_start{i_stim} = d.stimvol{i_stim}+2;
% indx_end{i_stim} = d.stimvol{i_stim}+9;
% 
% 
% 
%     for i_trial=1:length(indx_start{i_stim})
%         instances{i_stim}(:,i_trial) = mean(d.roi{1}.tSeries(:,indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)),2);
%         indx{i_stim}=[indx{i_stim},indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)];
% 
%     end
%     in = [in,instances{i_stim}];
%     dir =[dir,i_stim*ones(1,size(instances{i_stim},2))];
%     disppercent(i_stim/length(d.stimvol));
% end
%  
% 
% test_idx=reshape(1:length(dir(:)),length(dir(:))/length(d.stimvol),length(d.stimvol));
% disppercent(inf);
   

   
% for i=1:size(test_idx,1)
%     train_idx(i,:)=setdiff(1:length(dir),test_idx(i,:));
% end

%pick out the eventstrings and average if requested
idx = find(lab>0);
if params.averageEvent %average across the TRs for each event
    for i=1:size(idx,2)
        m_(:,i)=mean(d.roi{1}.tSeries(:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1),2);
    end
    d.roi{1}.tSeries=m_;clear m_
    run=run(idx);
    lab = lab(idx);
elseif params.eventLength>1 %create instance from each TR in the stim duration
    for i=1:size(idx,2)
        m_(:,idx(i):idx(i)+params.eventLength-1)=d.roi{1}.tSeries(:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1);
        l_(idx(i):idx(i)+params.eventLength-1)=repmat(lab(idx(i)),1,params.eventLength);
        r_(idx(i):idx(i)+params.eventLength-1)=repmat(run(idx(i)),1,params.eventLength);
    end
    d.roi{1}.tSeries=m_(l_>0);
    lab = l_(l_>0);
    run=r_(l_>0);
    clear m_ l_ r_
else %create instance from individual TRs related to each TR
    d.roi{1}.tSeries = d.roi{1}.tSeries(:,idx+params.hdLag);
    run=run(idx);
    lab = lab(idx);
end


% works out which roi coords are indexed by the spotlights
disppercent(-inf, 'Creating 4D TSeries');
for i=1:length(d.roi{1}.scanCoords)
mm(d.roi{1}.scanCoords(1,i),d.roi{1}.scanCoords(2,i),d.roi{1}.scanCoords(3,i),:) = d.roi{1}.tSeries(i,:);
disppercent(i/length(d.roi{1}.scanCoords));
end
disppercent(inf);



warning('off');

%performs leave one out classifcation for each spotlight
disppercent(-inf,'Classifying based on spotlight....');
% m=cell(1,length(idx));
% m_2=cell(1,length(idx));
% s=cell(1,length(idx));
% class_corr_diag=cell(1,length(idx));
% class_corr_linear=cell(1,length(idx));
class_acc_diag=nan(1,size(d.roi{1}.scanCoords,2));
% class_acc_linear=nan(1,length(idx));

% mm= mm(:,:,:,lab>0);
% run = run(lab>0);
% lab=lab(lab>0);

for i_sphere=1:size(d.roi{1}.scanCoords,2)
    
    idx=repmat(d.roi{1}.scanCoords(:,i_sphere),1,size(offsets,2))+offsets;
    [~,j] = find((idx(1,:)<min(d.roi{1}.scanCoords(1,:)) | idx(1,:)>max(d.roi{1}.scanCoords(1,:))) | (idx(2,:)<min(d.roi{1}.scanCoords(2,:)) | idx(2,:)>max(d.roi{1}.scanCoords(2,:))) | (idx(3,:)<min(d.roi{1}.scanCoords(3,:)) | idx(3,:)>max(d.roi{1}.scanCoords(3,:))));
    idx=idx(:,setdiff(1:size(idx,2),j));


    for i_vol=1:size(idx,2)
    xxx(i_vol,:)=mm(idx(1,i_vol),idx(2,i_vol),idx(3,i_vol),:);
    end
[I,~]=find(xxx>0);
xxx=xxx(unique(I),:);

for i=1:size(cc.runTransition,1)
    acc(i)=mean(lab(run==i)'==classify(xxx(:,run==i)',xxx(:,run~=i)',lab(run~=i),'diagLinear'));
end

%     whichstim=[1 3];
% 
%     ncat=length(d.stimvol);
% 
% 
%     for t_idx=1:length(indx_start{whichstim(1)})
%         test_patt=xxx(:,test_idx(t_idx,:));
%         test_lab=dir(test_idx(t_idx,:));
%         train_patt=xxx(:,train_idx(t_idx,:));
%         train_lab=dir(train_idx(t_idx,:));
% 
%         class_lab=classify(test_patt',train_patt',train_lab,'diagLinear');
%         class_corr_diag{i_sphere}(:,t_idx)= class_lab'==test_lab;
% 
%         class_lab=classify(test_patt',train_patt',train_lab,'otherlinear');
%         class_corr_linear{i_sphere}(:,t_idx)= class_lab'==test_lab;
% 
%     end
class_acc_diag(i_sphere)=mean(acc);
%  class_acc_diag(i_sphere)=mean(class_corr_diag{i_sphere}(:));
%  class_acc_linear(i_sphere)=mean(class_corr_linear{i_sphere}(:));

 xxx=[];

disppercent(i_sphere/length(d.roi{1}.scanCoords));
end

disppercent(inf);


%reshapes the data  and saves into into an img/hdr file
m_d=nan(viewGet(view,'scandims'));
% m_l=nan(viewGet(view,'scandims'));
    for i=1:length(d.roi{1}.scanCoords)
        m_d(d.roi{1}.scanCoords(1,i),d.roi{1}.scanCoords(2,i),d.roi{1}.scanCoords(3,i))=class_acc_diag(i);
%         m_l(d.roi{1}.scanCoords(1,i),d.roi{1}.scanCoords(2,i),d.roi{1}.scanCoords(3,i))=class_acc_linear(i);
    end
   
return


m_d_f=sprintf('rad_%i_acc_diag_%s.img',radius,hem);
% m_l_f=sprintf('rad_%i_acc_linear_%s.img',radius,hem);
% keyboard
cbiWriteNifti(m_d_f,m_d);
% cbiWriteNifti(m_l_f,m_l);

function offsets = make_sphere(radius);

% sphere.
offsets = [];

% we want to allow for the possibility of non-integer
% radii. however, we need the CEIL call here to convert the
% OFFSETS offsets into integers (because we are going to
% use these to index into the MASK matrix)
ceil_radius = ceil(radius);

for x = -1*ceil_radius:ceil_radius
  for y = -1*ceil_radius:ceil_radius
    for z = -1*ceil_radius:ceil_radius
      if sqrt(x^2+y^2+z^2) <= radius
        % append the current offset-triplet to the next row
        offsets(end+1,:) = [x y z];
      end
    end
  end
end

offsets=offsets';