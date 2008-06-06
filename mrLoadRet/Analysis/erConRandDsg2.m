function dsgMRand = erConRandDsg2(dsgM)
% This function randomises the design matrix as follows:
%   ones in the first column of the matrix stay put;
%   for all subsequent columns, ones are randomised: a one may
%   occur OFFSET_PROBE_TR TRs after the occurrence of a one in
%   the first column;
%   the number of ones in every column is preserved.

% For this regime, a probe for conditions >= 2 may only occur
% OFFSET_PROBE_TR TRs after the occurrence of a probe for
% condition 1.
OFFSET_PROBE_TR = 4;

dsgMRand = zeros(size(dsgM,1),size(dsgM,2));

for ii = 1:size(dsgM,2), countsCond(ii) = sum(dsgM(:,ii)); end

iiCandidateProbes = find(dsgM(:,1) == 1) + OFFSET_PROBE_TR;
iiCandidateProbes = iiCandidateProbes(randperm(length(iiCandidateProbes)));

dsgMRand(:,1) = dsgM(:,1);

countsCond(1) = 0; % fudge to make the following loop easy to code
for ii = 2:size(dsgMRand,2)

  dsgMRand(iiCandidateProbes((countsCond(ii-1)+1):(countsCond(ii-1)+countsCond(ii))),ii) = 1;
end

return;

