% dispmotioncorrect.m
%
%      usage: dispmotioncorrect.m()
%         by: justin gardner
%       date: 03/03/05
%    purpose: check motion correction
%
function retval = dispmotioncorrect(d)

if (nargin ~= 1)
  help dispmotioncorrect;
  return
end

titlestr = d.expname;
  
% get motion correct
titlestr = sprintf('%s scan=%i',d.expname,d.scanNum);
view = newView('Volume');
motioncorrect = viewGet(view,'transforms',d.scanNum,d.groupNum);
if ~isempty(motioncorrect)
  for i = 1:length(motioncorrect)
    d.motioncorrect.mat(i,:,:) = motioncorrect{i};
  end
else
  disp(sprintf('No associated motion correction transforms'));
  return
end

selectGraphWin;
plots = 3;plotnum = 0;

% now plot the fields
plotnum = plotnum+1;
subplot(plots,1,plotnum);
plot(squeeze(d.motioncorrect.mat(:,1,4)),'k.-');hold on
plot(squeeze(d.motioncorrect.mat(:,2,4)),'r.-');
plot(squeeze(d.motioncorrect.mat(:,3,4)),'c.-');
if (isfield(d.motioncorrect,'allscansmat'))
  plot(squeeze(d.motioncorrect.allscansmat(1,4)),'ko');
  plot(squeeze(d.motioncorrect.allscansmat(2,4)),'ro');
  plot(squeeze(d.motioncorrect.allscansmat(3,4)),'co');
else
%  disp(sprintf('UHOH: Undefined allscan matrix'));
end
legend('x','y','z');
title(sprintf('%s Translation parameters',titlestr));
xlabel('Image number');
ylabel('Translation (voxels)');
xaxis(1,size(d.motioncorrect.mat,1));
hline(0);

plotnum = plotnum+1;
subplot(plots,1,plotnum);
plot(squeeze(d.motioncorrect.mat(:,1,2)),'k.-');hold on
plot(squeeze(d.motioncorrect.mat(:,1,3)),'r.-');
plot(squeeze(d.motioncorrect.mat(:,2,3)),'c.-');
if (isfield(d.motioncorrect,'allscansmat'))
  plot(squeeze(d.motioncorrect.allscansmat(1,2)),'ko');
  plot(squeeze(d.motioncorrect.allscansmat(1,3)),'ro');
  plot(squeeze(d.motioncorrect.allscansmat(2,3)),'co');
end
legend('xy','xz','yz');
title('Shear parameters');
xlabel('Image number');
ylabel('Shear');
xaxis(1,size(d.motioncorrect.mat,1));
hline(0);

% now convert to axis angle and plot the angle
for i = 1:size(d.motioncorrect.mat,1)
  R = squeeze(d.motioncorrect.mat(i,1:3,1:3));
  theta(i) = 2 * acosd(0.5 * sqrt(trace(R)+1));
end
plotnum = plotnum+1;
subplot(plots,1,plotnum);
plot(theta,'ko');
ylabel('Rotation (degrees)');
xlabel('Image number');