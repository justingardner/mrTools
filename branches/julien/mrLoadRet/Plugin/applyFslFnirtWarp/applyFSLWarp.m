% warpedData = applyFSLWarp(data, warpCoefFilename, tempInputFilename, tempRefFilename)
%
%   applies non-linear FSL registration to 3/4D data
%
% jb 23/07/2010
%
% $Id$

function warpedData = applyFSLWarp(data, warpCoefFilename, tempFilename, hdr, interpMethod)


cbiWriteNifti(tempFilename, data, hdr);
if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end

%convert matlab flags to FNIRT flags
switch(interpMethod)
   case 'nearest' %Nearest neighbor interpolation
      interpMethod = 'nn';
   case 'linear'%Linear interpolation (default)
      interpMethod = 'trilinear';
   case {'spline', 'cubic'}
      mrWarnDlg(['(applyFSLWarp) ' interpMethod ' not implemented in FNIRT, using ''sinc'''])
      interpMethod = 'sinc';
end

try
  command =  sprintf('applywarp --ref=%s --in=%s --warp=%s --interp=%s --out=%s', tempFilename, tempFilename, warpCoefFilename, interpMethod, tempFilename);
  disp(command);
  [s,w] = unix(command);

  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(uhoh) looks like we have encountered a problem while running the applywarp FSL command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the TFCE values
warpedData=cbiReadNifti(tempFilename);
