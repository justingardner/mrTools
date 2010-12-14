% performs Threshold-Free Cluster Enhamcement on a 3-D array using flsmaths

function thresholded_data = applyFslTFCE(data, tempFilename,verbose)

if ieNotDefined('verbose')
  verbose = 1;
end
if ieNotDefined('tempFilename')
  tempFilename = 'temp.img';
end

if any(size(data)<3);
   mrErrorDlg('Volumes must have at least 3 voxels in each dimension to perform TFCE');
end

fslPath = mrGetPref('fslPath');
if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrErrorDlg('(applyFslTFCE) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
end

cbiWriteNifti(tempFilename, data,[],[],[],[],verbose);
tfce_H = 2.0;
tfce_E = 0.5;
tfce_connectivity = 6;
try
  [s,w] = unix(sprintf('%s/fslmaths %s -tfce %.2f %.2f %d %s',fslPath,  tempFilename, tfce_H, tfce_E, tfce_connectivity, tempFilename));
  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(applyFslTFCE) There was a problem running the TFCE unix command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the TFCE values
thresholded_data=cbiReadNifti(tempFilename,[],[],[],verbose);

if all(thresholded_data(:)==0)
  oneTimeWarning('tfceOutputsZeros','(applyFslTFCE) There is a problem with fslmaths -tfce (it outputs only zeros). try using another version of FSL');
end