% addSherlockJob.m
%
%   Adds jobs to Sherlock queue, given the split file handle.
%
%      by: akshay jagadeesh
%    date: 07/13/2017
function addSherlockJob(splitName, sherlockSessionPath, suid)

% get rid of .mat filename
if ~isempty(findstr(splitName, '.mat'))
  splitName = splitName(1:end-4);
end

pRFsavename = splitName(1:findstr(splitName, 'split')-2);
whichSplit = splitName(findstr(splitName, 'split')+5);

disp('Generating batch scripts');
system(sprintf('sh ~/proj/mrTools/mrLoadRet/Plugin/pRF/Sherlock/generateBatchScripts.sh "%s" "%s" "%s" "%s"', pRFsavename, sherlockSessionPath, suid, whichSplit))

disp('Transferring batch scripts to Sherlock and running');
system(sprintf('rsync -q Splits/Scripts/%s.sbatch %s@sherlock.stanford.edu:%s/Splits/Scripts/.', splitName, suid, sherlockSessionPath));
system(sprintf('ssh %s@sherlock.stanford.edu "cd %s/Splits/Scripts/; sbatch %s.sbatch"', suid, sherlockSessionPath, splitName));

