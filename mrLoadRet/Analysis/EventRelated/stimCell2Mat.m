% stimCell2Mat.m
%
%      usage: stimCell2Mat(stimOnsets, stimDurations{iStim}, framePeriod)
%         by: julien besle 
%       date: 05/12/2010
%    purpose: Transforms stimulus onsets/durations in cell array of vectors form into a stimulus matrix
%              $Id$

function [stimMatrix,runTransitions] = stimCell2Mat(stimOnsets, stimDurations,runTransitions)

if ~ieNotDefined('stimDurations')
  %apply stimDuration
  for iStim = 1:length(stimOnsets)
    stimOnsets{iStim} = reshape(stimOnsets{iStim},1,numel(stimOnsets{iStim})); %make sure stimOnsets{iStim} and stimDurations{iStim} 
    stimDurations{iStim} = reshape(stimDurations{iStim},1,numel(stimDurations{iStim})); %are row vectors
    maxDuration = max(stimDurations{iStim});
    stimNum = length(stimDurations{iStim});
    stimPresent = (cumsum(ones(maxDuration,stimNum)) ./ repmat(stimDurations{iStim},maxDuration,1))<=1;
    stimOnsets{iStim} = repmat(stimOnsets{iStim},maxDuration,1);
    stimOnsets{iStim} = stimOnsets{iStim} + cumsum(stimPresent) -1;
    stimOnsets{iStim} = stimOnsets{iStim}(stimPresent);
  end
end
if ieNotDefined('runTransitions')
  runTransitions = [1 1];
  for iStim = 1:length(stimOnsets)
    if ~isempty(stimOnsets{iStim})
      runTransitions(2) = max(runTransitions(2),stimOnsets{iStim}(end));
    end
  end
end
 

%fill the stim matrix
stimMatrix = zeros(runTransitions(end,2)-runTransitions(1,1)+1,length(stimOnsets));
for iRun = 1:size(runTransitions,1)
  thisStimMatrix = zeros(runTransitions(iRun,2)-runTransitions(iRun,1)+1,length(stimOnsets));
  for iStim = 1:length(stimOnsets)
    %only use stimOnsets{iStim}s that are within this runs volume numbers
    thisStimMatrix( stimOnsets{iStim}(stimOnsets{iStim}>=runTransitions(iRun,1) & stimOnsets{iStim}<=runTransitions(iRun,2) )- runTransitions(iRun,1)+1,iStim ) = 1;
  end
  stimMatrix(runTransitions(iRun,1):runTransitions(iRun,2),:)=thisStimMatrix;
end