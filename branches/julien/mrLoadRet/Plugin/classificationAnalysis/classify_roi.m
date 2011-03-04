function out = classify_roi(view,d,params,roi_n,select_vox,numShuff)
% Interrogator fucntion to perform Matlab and SVM classifcation on a given
% ROI. Orders the voxels by coherence for the localiser, but this can be
% changed to any other statistic by changing the code. The number of voxels
% used is also specified within the code, and can either be a single value
% or a vector that the code will cycle through.

% The code can also average voxels together based on the phase value from
% the localizer, but most people won't need to do that.


% fprintf('\n Preparing Data\n')

% for i=1:length(roi_n)
% r_coord{i}=getROICoordinates(view,roi_n(i),params.scanNum);
% end

if ~exist('numShuff','var'),numShuff=0;end

[xxx r_idx r_coord] = loadScanRois(view,params.scanNum,roi_n);
if select_vox
    overlay = selectSortOverlay(view);
    [~, thisOverlay] = isoverlay(overlay);
    curGroupName = viewGet(view,'groupName');
    curGroupNum =  viewGet(view,'groupNum',curGroupName);
    fromGroupNum = viewGet(view,'groupNum',thisOverlay.groupName);
  	scan2scan = viewGet(view,'scan2scan',1,curGroupNum,1,fromGroupNum);
    scandims = viewGet(view, 'scandims', 1);
    if ~all(all((scan2scan - eye(4))<1e-6)) %check if we're in the scan space
      %transform values in current scan space
      [Ycoords,Xcoords,Zcoords] = meshgrid(1:scandims(2),1:scandims(1),1:scandims(3));
%       for iScan = 1:length(scanList)
%          fprintf(1,['Overlay ' num2str(iOverlay) ', Scan ' num2str(iScan) '\n']);
         overlay.data = getNewSpaceOverlay(thisOverlay.data{1}, scan2scan,Xcoords,Ycoords,Zcoords);
%       end
   else
      overlay.data = thisOverlay.data;
    end
   for r = 1:size(r_coord,2)
       for v = 1:size(r_coord{r},2)
           sort_data{r}(v)=overlay.data{1}(r_coord{r}(1,v),r_coord{r}(2,v),r_coord{r}(3,v));
       end
   end
end


cc = viewGet(view,'concatinfo',params.scanNum);
lab = [];
for i=1:length(d.stimvol)
    lab(d.stimvol{i})=i;
end

run=[];
for i=1:size(cc.runTransition,1)
run(cc.runTransition(i,1):cc.runTransition(i,2))=i;
end

xxx.data = squeeze(xxx.data);

%conversion to % signal change
xxx.t_mean = mean(xxx.data,2);
xxx.data = 100*(xxx.data-repmat(xxx.t_mean,1,size(xxx.data,2)))./repmat(xxx.t_mean,1,size(xxx.data,2));

disp('(classify_roi) Creating patterns for classification')

%pick out the eventstrings and average if requested
idx = find(lab>0);
if params.averageEvent %average across the TRs for each event
    for i=1:size(idx,2)
        m_(:,i)=mean(xxx.data(:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1),2);
    end
    xxx.data=m_;clear m_
    run=run(idx);
    lab = lab(idx);
elseif params.eventLength>1 %create instance from each TR in the stim duration
    m_=[];
    l_=[];
    r_=[];
%     for i=1:size(idx,2)
%         m_(:,idx(i):idx(i)+params.eventLength-1)=xxx.data(:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1);
%         l_(idx(i):idx(i)+params.eventLength+-1)=repmat(lab(idx(i)),1,params.eventLength);
%         r_(idx(i):idx(i)+params.eventLength+-1)=repmat(run(idx(i)),1,params.eventLength);
%     end
    for i=1:size(idx,2)
        m_(:,idx(i):idx(i)+params.eventLength-1)=xxx.data(:,idx(i)+params.hdLag:idx(i)+params.eventLength+params.hdLag-1);
        l_(idx(i):idx(i)+params.eventLength+-1)=repmat(lab(idx(i)),1,params.eventLength);
        r_(idx(i):idx(i)+params.eventLength+-1)=repmat(run(idx(i)),1,params.eventLength);
    end
    xxx.data=m_(:,l_>0);
    lab = l_(l_>0);
    run=r_(l_>0);
    clear m_ l_ r_
else %create instance from individual TRs related to each TR
    xxx.data = xxx.data(:,idx+params.hdLag);
    run=run(idx);
    lab = lab(idx);
end

disp('(classify_roi) Classifying')

for r = 1:size(r_idx,2)
    patt=xxx.data(r_idx{r},:);
    
    if (select_vox==0) 
        patt_sort=1:size(patt,1);
    elseif (select_vox <  size(patt,1))
       [~, patt_sort] =sort(sort_data{r},'descend');
       patt_sort=patt_sort(1:select_vox);
    elseif (select_vox >= size(patt,1))
       [~, patt_sort] =sort(sort_data{r},'descend');   
    end
    
    class_lab=[];
    acc=[];
    svm_lab=[];
    svm_weight=[];
    svm_acc=[];
    for i=1:size(cc.runTransition,1)
    
    class_lab{r}(run==i)=classify(patt(patt_sort,run==i)',patt(patt_sort,run~=i)',lab(run~=i),'diagLinear');
    acc{r}(i)=mean(lab(run==i)==class_lab{r}(run==i));
    [svm_lab{r}(run==i) svm_weight{r}(:,:,:,i)]=classifyWithSvm(patt(patt_sort,run==i),patt(patt_sort,run~=i),lab(run~=i));
    svm_acc{r}(i)=mean(lab(run==i)==svm_lab{r}(run==i));
    end

m_acc{r}=mean(acc{r});
m_sacc{r}=mean(svm_acc{r});
svm_mean_weight{r}=squeeze(sum(svm_weight{r}(:,end,:,:),4));
[~,max_w{r}] = max(svm_mean_weight{r});
end

for r=1:size(r_idx,2)
for i=1:size(svm_mean_weight{r},1)
    for j=1:size(svm_mean_weight{r},1)
    svm_count{r}(i,j)=sum(svm_lab{r}(lab==i)==j)/sum(lab==i);
    end
end
end

  
if numShuff
    disp('(classify_roi) Runninng shuffled classification')
for r = 1:size(r_idx,2)
    patt=xxx.data(r_idx{r},:);
    
    if (select_vox==0) 
        patt_sort=1:size(patt,1);
    elseif (select_vox <  size(patt,1))
       [~, patt_sort] =sort(sort_data{r},'descend');
       patt_sort=patt_sort(1:select_vox);
    elseif (select_vox >= size(patt,1))
       [~, patt_sort] =sort(sort_data{r},'descend');   
    end
    for s=1:numShuff
    s_lab=lab(randperm(length(lab)));
    
        for i=1:size(cc.runTransition,1)
        s_acc(i)=mean(s_lab(run==i)'==classify(patt(patt_sort,run==i)',patt(patt_sort,run~=i)',lab(run~=i),'diagLinear'));
%         s_acc(i)=mean(s_lab(run==i)'==classify(xxx.data(r_idx{r},run==i)',xxx.data(r_idx{r},run~=i)',s_lab(run~=i),'diagLinear'));
        end

    sm_acc{r}(s)=mean(s_acc);
    th_95{r} = prctile(sm_acc{r},95);
    
end
end
end


figure('name','''tuning curves''')
for r=1:size(r_idx,2)
    colors = rainbow_colors(size(svm_count{r},1));
    subplot(1,size(r_idx,2),r)
    for i=1:size(svm_count{r},1)
        plot(svm_count{r}(i,:),'color',colors(i,:))
        hold on
    end
       axis([1 size(svm_count{r},1) 0 1])
       set(gca,'XTick',1:size(svm_count{r},1))
        set(gca,'XTickLabel',d.stimNames)
        legend(d.stimNames)
        line([1 size(svm_count{r},1)],[1/size(svm_count{r},1) 1/size(svm_count{r},1)],'Color','k')
        hold on
        if numShuff
        line([1 size(svm_count{r},1)],[th_95{r} th_95{r}],'Color','k','linestyle','--');
        hold on
        end
      title(sprintf('%s, Mean accuracy %f/%f, n-vox %i',viewGet(view,'roiname',r),m_acc{r},m_sacc{r},min([select_vox length(r_coord{r})])))  
end


keyboard

return
% This part decides the number of voxels and bins to use. To skip the phase
% bin averaging, set bin_vector to 0;
% n_vector = 300%[100:100:min(1000,size(roi_t,1)),size(roi_t,1)]
% bin_vector= [0,round(logspace(log10(3000),log10(3),10))]
% % 
% % 
% for i = 1:length(n_vector);
% 
% 
%     for j = 1:length(bin_vector)
%         
%         n = n_vector(i);
%         bin.n = bin_vector(j);
%         
%             roi_t_n = prepare_data(roi_t,roi_ph,n,bin);
%             [a, b, c] = classify_roi(roi_t_n,d.stimvol,d.stimNames,roi,n,bin,1);
% 
%             class(i,j)=a.acc;
%             class_o(i,j)=b.acc;
%             svm(i,j)=c.acc;
%     end
%    
% end

% 
% 
% bin.shuff=0;
% 
% roi_ph_2 = roi_ph(randperm(length(roi_ph)));
% if 1%bin.shuff;
% fprintf('\nShuffling Phase Values\n')
% 
% for i = 1:length(n_vector);
% 
%     for j = 1:length(bin_vector)
%         
%         n = n_vector(i);
%         bin.n = bin_vector(j);
%         
%             roi_t_n = prepare_data(roi_t,roi_ph_2,n,bin);
%             [a, b, c] = classify_roi(roi_t_n,d.stimvol,d.stimNames,roi,n,bin,0);
% 
%             class_shuff(i,j)=a.acc;
%             class_o_shuff(i,j)=b.acc;
%             svm_shuff(i,j)=c.acc;
%     end
%    
% end
% 
% fprintf('\nShuffling Phase Values: Run 2\n')
% roi_ph_3 = roi_ph(randperm(length(roi_ph)));
% 
% 
% for i = 1:length(n_vector);
% 
%     for j = 1:length(bin_vector)
%         
%         n = n_vector(i);
%         bin.n = bin_vector(j);
%         
%             roi_t_n = prepare_data(roi_t,roi_ph_3,n,bin);
%             [a, b, c] = classify_roi(roi_t_n,d.stimvol,d.stimNames,roi,n,bin,0);
% 
%             class_shuff_2(i,j)=a.acc;
%             class_o_shuff_2(i,j)=b.acc;
%             svm_shuff_2(i,j)=c.acc;
%     end
%    
% end
% end
% 
% save(['class_res_',date,'_',roi{1}.name],'svm*','class*','n_vector','bin_vector')
% bin.n=0;
% for n_shuff=1:1000
%     roi_t_n = prepare_data(roi_t,roi_ph,n,bin);
%     [a, b, c] = classify_roi_null(roi_t_n,d.stimvol,d.stimNames,roi,n,bin,0);
%     class_null(n_shuff)=a.acc;
%     class_o_null(n_shuff)=b.acc;
%     svm_shuff_null(n_shuff)=c.acc;
% end
% 
% bin.n=6;
% for n_shuff=1:1000
%     roi_t_n = prepare_data(roi_t,roi_ph,n,bin);
%     [a, b, c] = classify_roi_null(roi_t_n,d.stimvol,d.stimNames,roi,n,bin,0);
%     class_null_2(n_shuff)=a.acc;
%     class_o_null_2(n_shuff)=b.acc;
%     svm_shuff_null_2(n_shuff)=c.acc;
% end


keyboard

function roi_t_n = prepare_data(roi_t,roi_ph,n,bin)
if size(roi_t,1)<n, n=size(roi_t,1);end

roi_t_n=(roi_t(1:n,:));
roi_ph_n=(roi_ph(1:n));

if bin.shuff
    roi_ph_n = roi_ph_n(randperm(length(roi_ph_n)));
end

if bin.n
    [ph_n ph_x] = hist(roi_ph_n,bin.n);
    ph_x = ph_x-ph_x(1);
    for i=1:length(ph_x)-1
        tmp(i,:) = mean(roi_t_n(roi_ph_n >= ph_x(i) & roi_ph_n < ph_x(i+1),:),1);
    end
    tmp(length(ph_x),:) = mean(roi_t_n(roi_ph_n >= ph_x(length(ph_x)) & roi_ph_n < 2*pi,:),1);
    [nan_x nan_y] = find(isnan(tmp));
    tmp = tmp(setdiff(1:size(tmp,1),unique(nan_x)),:);
    roi_t_n = tmp;
    
end


% function [class class_o, svm] = classify_roi(roi_t_n,stimvol,stimNames,roi,n,bin,disp)
% 
% svm_flag=1;
% indx=[];
% for i_stim=1:length(stimvol)
% indx{i_stim}=[];
% indx_start{i_stim} = stimvol{i_stim}+2;
% indx_end{i_stim} = stimvol{i_stim}+9;
% 
% % keyboard
% 
%     for i_trial=1:length(indx_start{i_stim})
% 
% %       instances{i_stim}(:,i_trial) = mean(d.roi{1}.tSeries(1:n,indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)),2);
%         instances{i_stim}(:,i_trial) = mean(roi_t_n(:,indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)),2);
%         indx{i_stim}=[indx{i_stim},indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)];
% 
%     end
% end
% 
% warning('off')
%         
%         whichstim=[1:length(stimvol)];
%         ncat=length(whichstim);
% 
% 
%         for test_idx=1:length(indx_start{whichstim(1)})
%             train_idx=setdiff(1:length(indx_start{whichstim(1)}),test_idx);
% 
% 
%             test_patt=[];
%             train_patt=[];
%             lab=[];
%             test_lab=[];
% 
%             lin_test_patt=[];
%             lin_test_lab=[];
%             lin_train_patt=[];
%             lin_train_lab=[];
% 
%             for i_stim=1:length(whichstim)
%                 test_patt=[test_patt,instances{whichstim(i_stim)}(:,test_idx)];
%                 train_patt = [train_patt,instances{whichstim(i_stim)}(:,train_idx)];
%                 lab=[lab,whichstim(i_stim)*ones(1,length(train_idx))];
%                 test_lab=[test_lab;whichstim(i_stim)];
%             
% 
% %                 idx_temp=indx_start{i_stim}(test_idx):indx_end{i_stim}(test_idx);
% %                 lin_test_patt=[lin_test_patt,roi_t_n(1:n,idx_temp)];
% %                 lin_test_lab=[lin_test_lab,i_stim*ones(1,length([indx_start{i_stim}(test_idx):indx_end{i_stim}(test_idx)]))];
% %                 lin_train_patt=[lin_train_patt,roi_t(1:n,setdiff(indx{1},idx_temp))];
% %                 lin_train_lab=[lin_train_lab,i_stim*ones(1,length(setdiff(indx{1},idx_temp)))];
%             end
% %              test_patt=test_patt(:,[2 1]);
% %              test_lab=test_lab([2 1]);
%              
%                     % x1=train_patt(:,lab==svmstim(1));
%                     % x2=train_patt(:,lab==svmstim(2));
%                     % svm=getsvm(x1',x2');
%                     % svm_class{test_idx}=getsvm(test_patt',svm);
%                     % svm_acc(test_idx)=mean((svm_class{test_idx}>0)==(test_lab==test_lab(1))');
% 
% 
%             class_lab_o(:,test_idx)=classify_ab(test_patt',train_patt',lab,'otherLinear');
%             class_lab(:,test_idx)=classify_ab(test_patt',train_patt',lab,'diagLinear');
%             svm_lab(:,test_idx)=classifyWithSvm(test_patt,train_patt,lab);
% 
%         %     try
%         %         class_lab_lin=classify(test_patt',train_patt',lab);
%         %         acc_lin(test_idx)=mean(class_lab_lin==test_lab);
%         %     catch
%         %     end
% 
%             class_corr(:,test_idx)=class_lab(:,test_idx)==test_lab;
%             class_corr_o(:,test_idx)=class_lab_o(:,test_idx)==test_lab;
%             class_corr_svm(:,test_idx)=svm_lab(:,test_idx)==test_lab;
% 
% %             if svm_flag
% %             for j = 1:(ncat-1)
% %                 for k = (j+1):ncat
% %                   % get all instances
% %                   x1 = train_patt(:,lab==whichstim(j));
% %                   x2 = train_patt(:,lab==whichstim(k));
% %                   % remove the one instance
% %                   trainset1 = setdiff(1:size(x1,1),test_idx);
% %                   trainset2 = setdiff(1:size(x2,1),test_idx);
% %                   % create the classifer
% %                   svm{test_idx,j,k} = getsvm(x1',x2');%,kernelfun,kernelargs,C);
% %                   % store the weights
% %                   weights(test_idx,j,k,:) = svm{test_idx,j,k}.w;
% %                 end
% %             end
% %             % now test each one of the left out response
% %               for teststim = 1:ncat
% %                 % get each one of the test stimuli in turn
% %                 testvec =  test_patt(:,test_lab==whichstim(teststim));
% %                 for j = 1:ncat
% %                   for k = j:ncat
% %                     if (j == k) 
% %                       % classification against oneself we define as 0
% %                       classifierout{test_idx,teststim}(j,k) = 0;
% %                     else
% %                       % use the svm, calculated above w/out the
% %                       % test stimulus to classify the stimulus against
% %                       % each one of the other stimuli categories
% %                       classifierout{test_idx,teststim}(j,k) = getsvm(testvec',svm{test_idx,j,k});
% %                       % in the matrix of outputs, the classification of the
% %                       % other stimulus categories vs. this one, is symmetric
% %                       % i.e. if we know the classification for stimulus 1 vs
% %                       % stimulus 2, then the classification for stimulus 2 vs
% %                       % stimulus 1 is just the negative of that.
% %                       classifierout{test_idx,teststim}(k,j) = -classifierout{test_idx,teststim}(j,k);
% %                     end
% %                   end
% %                   % now calculate the sum of the classifier against all other
% %                   % stimuli types
% %                     classifierout{test_idx,teststim}(j,ncat+1) = sum(classifierout{test_idx,teststim}(j,1:ncat));
% %                 end
% %                 % make a tuning curve--this gives the summed classifier outputs
% %                 % for each possible stimulus types.
% % 
% %                 temp_tuningcurve = squeeze(classifierout{test_idx,teststim}(:,ncat+1));
% %                 tuningcurve(test_idx,teststim,:) = temp_tuningcurve;
% %                 % the maximum point on this tuning curve is the best guess at
% %                 % how to classify this instance
% %                 classification(test_idx,teststim) = find(first(max(temp_tuningcurve))==temp_tuningcurve);
% % 
% %                 svm_corr(test_idx,teststim) = classification(test_idx,teststim) == teststim;
% %               end
% %             end
% 
% 
%         end
% 
% class.acc=mean(class_corr(:));
% class.p = myBinomTest(sum(class_corr(:)),length(class_corr(:)),1/length(whichstim),'Greater');
% 
% class_o.acc=mean(class_corr_o(:));
% class_o.p = myBinomTest(sum(class_corr_o(:)),length(class_corr_o(:)),1/length(whichstim),'Greater');
% 
% if svm_flag
% svm.acc=mean(class_corr_svm(:));
% svm.p = myBinomTest(sum(class_corr_svm(:)),length(class_corr_svm(:)),1/length(whichstim),'Greater');
% else
% svm.acc=NaN;
% svm.p=NaN;
% end
% if disp
%     
%     fprintf('\nResults for ROI %s: %i voxels included.',roi{1}.name,n);
% 
%     if bin.n
%         fprintf('\n%i bins used = bin width %d degrees.',bin.n,360/bin.n);
%         if bin.shuff,fprintf(' SHUFFLED phase values!')
%     end
%     end
% 
% 
% 
%     fprintf('\n\nClassification resuts for directions: '),fprintf('%i ',whichstim);
%     fprintf('\n\nresults for matlab classification\n');
%     fprintf('\nAcc = %f, p = %.3f',class.acc,class.p);
% 
% 
%     % [class_h class_p class_ci] =ttest(class_corr(:),1/ncat,0.05)
%     fprintf('\n\nresults for other matlab classification\n')
%     fprintf('\nAcc = %f, p = %.3f',class_o.acc,class_o.p);
% 
% 
%     % [class_h_o class_p_o class_ci_o] =ttest(class_corr_o(:),1/ncat,0.05)
% 
% 
%     fprintf('\n\nresults for SVM classification\n')
%     fprintf('\nAcc = %f, p = %.3f',svm.acc,svm.p);
%     % [svm_h svm_sig svm_ci] = ttest(svm_corr(:),1/ncat,0.05)
%     fprintf('\n');
% 
% %     for i_r =1:size(classification,2)
% % 
% %         r(i_r,:) = [sum(classification(:,i_r)==1) sum(classification(:,i_r)==2) sum(classification(:,i_r)==3) sum(classification(:,i_r)==4) sum(classification(:,i_r)==5) sum(classification(:,i_r)==6) sum(classification(:,i_r)==7) sum(classification(:,i_r)==8) sum(classification(:,i_r)==1)]/i_trial;
% %     end
% 
% 
%     for i_r = 1:size(class_lab_o,1)
%         for j = 1:length(whichstim)
%         r(i_r,j)=sum(class_lab_o(i_r,:)==j)./size(class_lab_o,2);
%         end
%     end
% 
% %     t=[0:(pi/4):2*pi];
% t = [1:size(r,1)];
% 
%     colors = rainbow_colors(size(r,1));
%     figure('name',[roi{1}.name,' ''tuning curves'''])
%     for i_r=1:size(r,1)
%         subplot(3,ceil(size(r,1)/2),i_r)
%         plot(t,r(i_r,:),'color',colors(i_r,:))
%         set(gca,'XTick',[i_r],'XTickLabel',stimNames([i_r]))
%         hold on
%         line([1 size(r,1)],[1/size(r,1) 1/size(r,1)],'color','k','LineStyle','--')
%         axis([1 size(r,1) 0 1])
% % %       polar(0,1,'w')
% % 
%         title(stimNames{i_r})
% %         hold on
% %          h(i_r) = polar(t,r(i_r,:));
% %         polar(0:.1:2*pi,repmat(1/8,1,length(0:.1:2*pi)),'--k')
% %         h(i_r) = polar(t,r(i_r,:));
% %         set(h(i_r),'color',colors(i_r,:));
% %         hold off
%     end
% 
%     for j = 1:size(r,1)
%         rr(j,:) = r(j,[j:end,1:j-1]);
%     end
%     m_r = mean(rr,1);
%     subplot(3,ceil(size(r,1)/2),i_r+1:3*ceil(size(r,1)/2))
%     plot(t,m_r,'k')
%     set(gca,'XTick',t,'XTickLabel',stimNames)
%     xlabel('Offset (degrees)')
%     ylabel('Chance of classification')
%     yaxis([0 1])
%     hold on
%     line([1 size(r,1)],[1/size(r,1) 1/size(r,1)],'color','k','LineStyle','--')
% 
%     
% 
% %     tt=t(1:8);
% %     rr=r(:,1:8);
% % 
% %     for i=1:size(r,1)
% %         t_dm(i,:)=tt-repmat(tt(i),1,length(tt));
% %     end
% % 
% %     for j=1:8
% %     t_dm_w(j,:)=t_dm(j,[j:end,1:j-1]);r_w(j,:)=rr(j,[j:end,1:j-1]);
% %     end
% % 
% %     r_c=r_w(:,[5:end,1:5]);
% %     r_c_m=mean(r_c,1);
% %     % t_l=[-180 -135 -90 -45 0 45 90 135 180];
% %     % figure,plot(t_l,r_c_m);
% %     figure
% %     polar(0,0.5,'w')
% %     hold on
% %     polar(0:pi/4:2*pi,r_c_m([1:9]))
% %     % ylabel('% classified');
% %     % xlabel('direction offset from actual');
% end
% 
% function [class class_o, svm] = classify_roi_null(roi_t_n,stimvol,stimNames,roi,n,bin,disp)
% 
% svm_flag=1;
% indx=[];
% for i_stim=1:length(stimvol)
% indx{i_stim}=[];
% indx_start{i_stim} = stimvol{i_stim}+2;
% indx_end{i_stim} = stimvol{i_stim}+9;
% 
% % keyboard
% 
%     for i_trial=1:length(indx_start{i_stim})
% 
% %       instances{i_stim}(:,i_trial) = mean(d.roi{1}.tSeries(1:n,indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)),2);
%         instances{i_stim}(:,i_trial) = mean(roi_t_n(:,indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)),2);
%         indx{i_stim}=[indx{i_stim},indx_start{i_stim}(i_trial):indx_end{i_stim}(i_trial)];
% 
%     end
% end
% 
% warning('off')
%         
%         whichstim=[1:length(stimvol)];
%         ncat=length(whichstim);
%         shuff = randperm(ncat*(length(indx_start{1})-1));
% 
%         for test_idx=1:length(indx_start{whichstim(1)})
%             train_idx=setdiff(1:length(indx_start{whichstim(1)}),test_idx);
% 
% 
%             test_patt=[];
%             train_patt=[];
%             lab=[];
%             test_lab=[];
% 
%             lin_test_patt=[];
%             lin_test_lab=[];
%             lin_train_patt=[];
%             lin_train_lab=[];
% 
%             for i_stim=1:length(whichstim)
%                 test_patt=[test_patt,instances{whichstim(i_stim)}(:,test_idx)];
%                 train_patt = [train_patt,instances{whichstim(i_stim)}(:,train_idx)];
%                 lab=[lab,whichstim(i_stim)*ones(1,length(train_idx))];
%                 test_lab=[test_lab;whichstim(i_stim)];
%             
% 
% %                 idx_temp=indx_start{i_stim}(test_idx):indx_end{i_stim}(test_idx);
% %                 lin_test_patt=[lin_test_patt,roi_t_n(1:n,idx_temp)];
% %                 lin_test_lab=[lin_test_lab,i_stim*ones(1,length([indx_start{i_stim}(test_idx):indx_end{i_stim}(test_idx)]))];
% %                 lin_train_patt=[lin_train_patt,roi_t(1:n,setdiff(indx{1},idx_temp))];
% %                 lin_train_lab=[lin_train_lab,i_stim*ones(1,length(setdiff(indx{1},idx_temp)))];
%             end
%             lab = lab(randperm(length(lab)));
% %              test_patt=test_patt(:,[2 1]);
% %              test_lab=test_lab([2 1]);
%              
%                     % x1=train_patt(:,lab==svmstim(1));
%                     % x2=train_patt(:,lab==svmstim(2));
%                     % svm=getsvm(x1',x2');
%                     % svm_class{test_idx}=getsvm(test_patt',svm);
%                     % svm_acc(test_idx)=mean((svm_class{test_idx}>0)==(test_lab==test_lab(1))');
% 
% 
%             class_lab_o(:,test_idx)=classify_ab(test_patt',train_patt',lab,'otherLinear');
%             class_lab(:,test_idx)=classify_ab(test_patt',train_patt',lab,'diagLinear');
%             svm_lab(:,test_idx)=classifyWithSvm(test_patt,train_patt,lab);
% 
%         %     try
%         %         class_lab_lin=classify(test_patt',train_patt',lab);
%         %         acc_lin(test_idx)=mean(class_lab_lin==test_lab);
%         %     catch
%         %     end
% 
%             class_corr(:,test_idx)=class_lab(:,test_idx)==test_lab;
%             class_corr_o(:,test_idx)=class_lab_o(:,test_idx)==test_lab;
%             class_corr_svm(:,test_idx)=svm_lab(:,test_idx)==test_lab;
% 
% %             if svm_flag
% %             for j = 1:(ncat-1)
% %                 for k = (j+1):ncat
% %                   % get all instances
% %                   x1 = train_patt(:,lab==whichstim(j));
% %                   x2 = train_patt(:,lab==whichstim(k));
% %                   % remove the one instance
% %                   trainset1 = setdiff(1:size(x1,1),test_idx);
% %                   trainset2 = setdiff(1:size(x2,1),test_idx);
% %                   % create the classifer
% %                   svm{test_idx,j,k} = getsvm(x1',x2');%,kernelfun,kernelargs,C);
% %                   % store the weights
% %                   weights(test_idx,j,k,:) = svm{test_idx,j,k}.w;
% %                 end
% %             end
% %             % now test each one of the left out response
% %               for teststim = 1:ncat
% %                 % get each one of the test stimuli in turn
% %                 testvec =  test_patt(:,test_lab==whichstim(teststim));
% %                 for j = 1:ncat
% %                   for k = j:ncat
% %                     if (j == k) 
% %                       % classification against oneself we define as 0
% %                       classifierout{test_idx,teststim}(j,k) = 0;
% %                     else
% %                       % use the svm, calculated above w/out the
% %                       % test stimulus to classify the stimulus against
% %                       % each one of the other stimuli categories
% %                       classifierout{test_idx,teststim}(j,k) = getsvm(testvec',svm{test_idx,j,k});
% %                       % in the matrix of outputs, the classification of the
% %                       % other stimulus categories vs. this one, is symmetric
% %                       % i.e. if we know the classification for stimulus 1 vs
% %                       % stimulus 2, then the classification for stimulus 2 vs
% %                       % stimulus 1 is just the negative of that.
% %                       classifierout{test_idx,teststim}(k,j) = -classifierout{test_idx,teststim}(j,k);
% %                     end
% %                   end
% %                   % now calculate the sum of the classifier against all other
% %                   % stimuli types
% %                     classifierout{test_idx,teststim}(j,ncat+1) = sum(classifierout{test_idx,teststim}(j,1:ncat));
% %                 end
% %                 % make a tuning curve--this gives the summed classifier outputs
% %                 % for each possible stimulus types.
% % 
% %                 temp_tuningcurve = squeeze(classifierout{test_idx,teststim}(:,ncat+1));
% %                 tuningcurve(test_idx,teststim,:) = temp_tuningcurve;
% %                 % the maximum point on this tuning curve is the best guess at
% %                 % how to classify this instance
% %                 classification(test_idx,teststim) = find(first(max(temp_tuningcurve))==temp_tuningcurve);
% % 
% %                 svm_corr(test_idx,teststim) = classification(test_idx,teststim) == teststim;
% %               end
% %             end
% 
% 
%         end
% 
% class.acc=mean(class_corr(:));
% class.p = myBinomTest(sum(class_corr(:)),length(class_corr(:)),1/length(whichstim),'Greater');
% 
% class_o.acc=mean(class_corr_o(:));
% class_o.p = myBinomTest(sum(class_corr_o(:)),length(class_corr_o(:)),1/length(whichstim),'Greater');
% 
% if svm_flag
% svm.acc=mean(class_corr_svm(:));
% svm.p = myBinomTest(sum(class_corr_svm(:)),length(class_corr_svm(:)),1/length(whichstim),'Greater');
% else
% svm.acc=NaN;
% svm.p=NaN;
% end
% 
% 
% 
% 
% 
% 
% % 
% % figure, hist(class_acc_s);
% % hold on
% % svm_pr = prctile(class_acc_s, 95);
% % 
% % line([svm_pr svm_pr],[0 300],'color','k');
% % line([class_acc class_acc],[0 300],'color','r');
% % 
% % keyboard