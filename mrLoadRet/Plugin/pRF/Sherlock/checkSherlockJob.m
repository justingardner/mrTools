% checkSherlockJob.m
% 
%

function isJobDone = checkSherlockJob(splitName)

curPath = pwd;
sherlockSessionPath = ['/share/PI/jlg/' curPath(findstr(curPath, 'data') : end)];
suid = mglGetParam('sunetID');

[~,out] = system(sprintf('ssh %s@sherlock.stanford.edu "if [ -f %s/Splits/Analysis/%s_Anal.mat ]; then echo exists; else echo doesNotExist; fi"', suid, sherlockSessionPath, splitName));

if strcmp(deblank(out), 'exists')
  disp('Analysis found. Job has successfully completed running! Pulling files to local machine now...');
  system(sprintf('rsync -q %s@sherlock.stanford.edu:%s/Splits/Analysis/%s_Anal.mat Splits/Analysis/.', suid, sherlockSessionPath, splitName));

end

keyboard
