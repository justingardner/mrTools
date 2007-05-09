function ext = niftiFileExtension()
% ext = niftiFileExtension()
%
% Uses 'niftiFileExtension' preference to determine *.img or *.nii
% Default: *.img if preference is not set.
%
% djh 7/2006

% Use niftiFileExtension preference or default to *.img
ext = mrGetPref('niftiFileExtension');
if isempty(ext)
    ext = '.img';
end

% Check valid nift file extension (either .nii or .img)
if ~strcmp(ext,'.img') & ~strcmp(ext,'.nii')
    mrErrorDlg('Invalid nifti file extension: ',ext);
end
