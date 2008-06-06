function dsgMRand = erConRandDsg1(dsgM)
%
%
%

  dsgMRand = zeros(size(dsgM,1),size(dsgM,2));

  for ii = 1:size(dsgM,2), countsCond(ii) = sum(dsgM(:,ii)); end

  aa = randperm(size(dsgMRand,1));
  dsgMRand(aa(1:countsCond(1)),1) = 1;
  for ii = 2:size(dsgMRand,2)

    dsgMRand(aa(countsCond(ii-1)+1:countsCond(ii-1)+1+countsCond(ii)),ii) = 1;
  end

return;

