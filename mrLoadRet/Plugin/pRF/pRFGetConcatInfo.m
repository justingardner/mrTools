% pRFGetConcatInfo.m
%
%        $Id:$ 
%      usage: pRFGetConcatInfo()
%         by: justin gardner
%       date: 11/22/15
%    purpose: Gets concatInfo - making up one for non-concat scans
%
function concatInfo = pRFGetConcatInfo(v,scanNum)

concatInfo = viewGet(v,'concatInfo',scanNum);

% if there is no concatInfo, then make one that will
% treat the scan as a single scan
if isempty(concatInfo)
  nFrames = viewGet(v,'nFrames',scanNum);
  concatInfo.isConcat = false;
  concatInfo.n = 1;
  concatInfo.whichScan = ones(1,nFrames);
  concatInfo.whichVolume = 1:nFrames;
  concatInfo.runTransition = [1 nFrames];
  concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  if length(concatInfo.totalJunkedFrames > 1)
    % first check for consistency in totalJunkedFrames
    if length(unique(concatInfo.totalJunkedFrames)) > 1
      disp(sprintf('(pRFFit) totalJunkedFrames are different for different members of component scans - could be an average in which different scans with different number of junked frames were removed. This could cause a problem in computing what the stimulus was for the average. The total junked frames count was: %s, but we will use %i as the actual value for computing the stimulus',num2str(concatInfo.totalJunkedFrames),floor(median(concatInfo.totalJunkedFrames))));
    end
    concatInfo.totalJunkedFrames = floor(median(concatInfo.totalJunkedFrames));
  end
else
  concatInfo.isConcat = true;
  if ~isfield(concatInfo,'totalJunkedFrames')
    concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  end
end

