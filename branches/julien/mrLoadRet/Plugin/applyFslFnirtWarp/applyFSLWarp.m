% warpedData = applyFSLWarp(data, warpCoefFilename, tempInputFilename, tempRefFilename)
%
%   applies non-linear FSL registration to 3/4D data
%
% jb 23/07/2010
%
% $Id$

function warpedData = applyFSLWarp(data, warpCoefFilename, tempFilename, hdr, interpMethod, verbose)


if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end
if ieNotDefined('verbose')
   verbose = true;
end

cbiWriteNifti(tempFilename, data, hdr,[],[],[],verbose);
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
  if verbose,disp(command);end;
  [s,w] = unix(command);

  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(applyFSLWarp) There was a problem running FSL applywarp command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the warped values
warpedData=cbiReadNifti(tempFilename);

if any(isnan(data(:))) %is there are NaNs, warp a mask of non-NaNs
  mask = ~isnan(data);
  warpedMask = logical(applyFSLWarp(mask, warpCoefFilename, tempFilename, hdr, 'nearest'));
  warpedData(~warpedMask) = NaN;
end

