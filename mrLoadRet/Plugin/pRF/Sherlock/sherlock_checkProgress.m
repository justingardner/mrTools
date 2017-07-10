function sherlock_CheckProgress

suid = 'akshayj';
system(sprintf('ssh %s@sherlock.stanford.edu "squeue -u %s"', suid, suid));
