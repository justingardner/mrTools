function [m_d m_c m_p]=classify_searchlight(view,d,radius,roiMask,params)
%This fucntion does a searchlight classification with a defined radius in
%an area defined by an ROI, which is normally a mask of the skullstripped
%brain. Has to be the run from the root directory of an MLR session.


% a function to calculate the offsets of the spotlight voxels relative to
% the centre of the sphere
offsets=make_sphere(radius);

if roiMask

% d.roi{1}.scanCoords = getROICoordinates(view,roiMask,params.scanNum);
% d.roi{1} = loadROITSeries(view,roiMask,params.scanNum);
[xxx, roiVoxelIndices, roiCoords] = loadScanRois(view,params.scanNum,viewGet(view,'roinum',roiMask));
d.roi{1}.tSeries=squeeze(xxx.data);
d.roi{1}.scanCoords = roiCoords;
clear xxx
end



warning('off');


% preprocess the functional data

cc = viewGet(view,'concatinfo',params.scanNum);

for ii=1:size(params.stimToEVmatrix,2)
    tmp.stimvol{ii}=[d.stimvol{[find(params.stimToEVmatrix(:,ii))]}];
end
    d.stimvol=tmp.stimvol;
    d.stimNames=params.EVnames;
lab = [];

for i=1:length(d.stimvol)
    lab(d.stimvol{i})=i;
end

run=[];
for i=1:size(cc.runTransition,1)
run(cc.runTransition(i,1):cc.runTransition(i,2))=i;
end

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
for i=1:length(d.roi{1}.scanCoords{1})
mm(d.roi{1}.scanCoords{1}(1,i),d.roi{1}.scanCoords{1}(2,i),d.roi{1}.scanCoords{1}(3,i),:) = d.roi{1}.tSeries(i,:);
disppercent(i/length(d.roi{1}.scanCoords{1}));
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
class_acc=nan(length(d.stimvol),size(d.roi{1}.scanCoords{1},2));
class_acc_diag=nan(1,size(d.roi{1}.scanCoords{1},2));
% class_acc_linear=nan(1,length(idx));

% mm= mm(:,:,:,lab>0);
% run = run(lab>0);
% lab=lab(lab>0);

for i_sphere=1:size(d.roi{1}.scanCoords{1},2)
    
    idx=repmat(d.roi{1}.scanCoords{1}(:,i_sphere),1,size(offsets,2))+offsets;
    [~,j] = find((idx(1,:)<min(d.roi{1}.scanCoords{1}(1,:)) | idx(1,:)>max(d.roi{1}.scanCoords{1}(1,:))) | (idx(2,:)<min(d.roi{1}.scanCoords{1}(2,:)) | idx(2,:)>max(d.roi{1}.scanCoords{1}(2,:))) | (idx(3,:)<min(d.roi{1}.scanCoords{1}(3,:)) | idx(3,:)>max(d.roi{1}.scanCoords{1}(3,:))));
    idx=idx(:,setdiff(1:size(idx,2),j));


    for i_vol=1:size(idx,2)
    xxx(i_vol,:)=mm(idx(1,i_vol),idx(2,i_vol),idx(3,i_vol),:);
    end
[I,~]=find(xxx>0);
xxx=xxx(unique(I),:);

class_lab=nan(1,length(run));
corr=nan(1,length(run));
for i=1:size(cc.runTransition,1)
    class_lab(run==i)=classify(xxx(:,run==i)',xxx(:,run~=i)',lab(run~=i),'diagLinear')';
    corr(run==i)=class_lab(run==i)==lab(run==i);    
end
acc_p(i_sphere) = binomTest(sum(corr),length(corr),1/length(d.stimvol),'Greater');
for i=1:length(d.stimvol)
    class_acc(i,i_sphere)=mean(corr(lab==i));
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
class_acc_diag(i_sphere)=mean(class_acc(:,i_sphere));
%  class_acc_diag(i_sphere)=mean(class_corr_diag{i_sphere}(:));
%  class_acc_linear(i_sphere)=mean(class_corr_linear{i_sphere}(:));

 xxx=[];

disppercent(i_sphere/length(d.roi{1}.scanCoords{1}));
end

disppercent(inf);


%reshapes the data  and saves into into an img/hdr file
m_d=nan(viewGet(view,'scandims'));
m_p=nan(viewGet(view,'scandims'));
m_c=nan([viewGet(view,'scandims'),length(d.stimvol)]);
% m_l=nan(viewGet(view,'scandims'));
    for i=1:length(d.roi{1}.scanCoords{1})
        m_d(d.roi{1}.scanCoords{1}(1,i),d.roi{1}.scanCoords{1}(2,i),d.roi{1}.scanCoords{1}(3,i))=class_acc_diag(i);
        m_c(d.roi{1}.scanCoords{1}(1,i),d.roi{1}.scanCoords{1}(2,i),d.roi{1}.scanCoords{1}(3,i),:)=class_acc(:,i);
        m_p(d.roi{1}.scanCoords{1}(1,i),d.roi{1}.scanCoords{1}(2,i),d.roi{1}.scanCoords{1}(3,i),:)=acc_p(i);
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