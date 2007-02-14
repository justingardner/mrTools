function [hdr,fid] = cbiWriteNiftiHeader(hdr,fname,no_overwrite,leave_open)
% [hdr,fid] = cbiWriteNiftiHeader(hdr,fname [,no_overwrite,leave_open])
% 
% Saves a NIFTI-1 file header using the hdr struct provided.
% See cbiReadNiftiHeader for details of header structure.
%
% WARNING: This function assumes that the header is valid - it does not check that
% all required fields are present or that the fields are congruent.
% A valid header struct can be generated with cbiCreateNiftiHeader(). 
%
% fname can be a file name or a file pointer.
% If fname is a file name: if the extension is:
%  - .hdr/.img: Will open or create the .hdr file and write to BOF.
%  - .nii: Will open or create a .nii file and write to BOF.
%  Unless the no_overwrite flag is nonzero, any existing file with the same name will
%  be deleted. (Use no_overwrite to replace the header data without destroying the
%  image data in a .nii file). This flag is ignored for dual file types (.img file is
%  not affected by this function).
%
% If fname is a file pointer: the data will be written to the file
% pointer at the current location (pointer should be at BOF)
%
% If fname is not specified or empty, will use the hdr.hdr_name field, if it exists.
% If this field is not set, fname must be specified.
%
% If leave_open is nonzero, will leave the file open for writing and return a 
% pointer to the open file (appropriate for writing image data to a .nii file)
% This flag is ignored for dual (hdr/img) file types.
% 

if (~exist('no_overwrite'))
  no_overwrite=0;
end
if (~exist('leave_open'))
  leave_open=0;
end
if (~exist('fname') | isempty(fname))
  if (~isfield(hdr,'hdr_name') | isempty(hdr.hdr_name))
    error('No file name specified!');
  else
    fname=hdr.hdr_name;
  end
end

if (isstr(fname))
  % Create a proper file name
  % Check name and file type.
  [pathstr,bname,ext,ver]=fileparts(fname);
  
  switch (ext)
   case '.nii'  
    hdr.single_file=1;
    hdr.hdr_name=fname;
    hdr.img_name=fname;
    % This hack is necessary as sprintf on Mac OS X is broken (refuses to print \0 character)
    hdr.magic=char(zeros(1,4));
    hdr.magic(1:3)=sprintf('%s','n+1');
   case {'.hdr','.img'}
    hdr.single_file=0;
    hdr.hdr_name=fullfile(pathstr,[bname '.hdr']);
    hdr.img_name=fullfile(pathstr,[bname '.img']);
    % This hack is necessary as sprintf on Mac OS X is broken (refuses to print \0 character)
    hdr.magic=char(zeros(1,4));
    hdr.magic(1:3)=sprintf('%s','ni1');
   case {'.gz','.Z'} % zipped
    error('No support for zipped NIFTI-1 format under Matlab.');
   otherwise
    error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img');
  end  
  if (~isfield(hdr,'endian'))
    hdr.endian='native';
  end
  if (hdr.single_file)
    fid=fopen(fname,'w',hdr.endian);
  else
    fid=fopen(hdr.hdr_name,'w',hdr.endian);
  end
else
  fid=fname;
end


% Write header with minimal checking

c=fwrite(fid,348,'int32');                   % sizeof_hdr
c=fwrite(fid,zeros(10,1),'char');           % data_type (unused)
c=fwrite(fid,zeros(18,1),'char');           % db_name (unused)
c=fwrite(fid,0,'int32');                    % extents (unused)
c=fwrite(fid,0,'short');                    % session_error (unused)
c=fwrite(fid,0,'char');                     % regular (unused)
c=fwrite(fid,hdr.dim_info,'char');          % dim_info (MRI slice ordering)
if (c~=1) error(['hdr.dim_info must be a single char (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.dim,'short');               % data array dimensions
if (c~=8) error(['hdr.dim must be 8 shorts (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.intent_ps,'float');        %
if (c~=3) error(['hdr.intent_ps must be 3 floats (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.intent_code,'short');
if (c~=1) error(['hdr.intent_code must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.datatype,'short');
if (c~=1) error(['hdr.datatype must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.bitpix,'short');
if (c~=1) error(['hdr.bitpix must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.slice_start,'short');
if (c~=1) error(['hdr.slice_start must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.pixdim,'float');
if (c~=8) error(['hdr.pixdim must be 8 floats (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.vox_offset,'float');
if (c~=1) error(['hdr.vox_offset must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.scl_slope,'float');
if (c~=1) error(['hdr.scl_slope must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.scl_inter,'float');
if (c~=1) error(['hdr.scl_inter must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.slice_end,'short');
if (c~=1) error(['hdr.slice_end must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.slice_code,'char');
if (c~=1) error(['hdr.slice_code must be a single char (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.xyzt_units,'char');
if (c~=1) error(['hdr.xyzt_units must be a single char (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.cal_max,'float');
if (c~=1) error(['hdr.cal_max must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.cal_min,'float');
if (c~=1) error(['hdr.cal_min must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.slice_duration,'float');
if (c~=1) error(['hdr.slice_duration must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.toffset,'float');
if (c~=1) error(['hdr.toffset must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,zeros(2,1),'int32');                % glmax, glmin (unused)
c=fwrite(fid,hdr.descrip,'char');
if (c~=80) error(['hdr.descrip must be 80 chars (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.aux_file,'char');
if (c~=24) error(['hdr.aux_file must be 24 chars (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.qform_code,'short');
if (c~=1) error(['hdr.qform_code must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.sform_code,'short');
if (c~=1) error(['hdr.sform_code must be a single short (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.quatern_b,'float');
if (c~=1) error(['hdr.quatern_b must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.quatern_c,'float');
if (c~=1) error(['hdr.quatern_c must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.quatern_d,'float');
if (c~=1) error(['hdr.quatern_d must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.qoffset_x,'float');
if (c~=1) error(['hdr.qoffset_x must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.qoffset_y,'float');
if (c~=1) error(['hdr.qoffset_y must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.qoffset_z,'float');
if (c~=1) error(['hdr.qoffset_z must be a single float (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.srow_x,'float');
if (c~=4) error(['hdr.srow_x must be 4 floats (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.srow_y,'float');
if (c~=4) error(['hdr.srow_y must be 4 floats (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.srow_z,'float');
if (c~=4) error(['hdr.srow_z must be 4 floats (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.intent_name,'char');
if (c~=16) error(['hdr.intent_name must be 16 chars (' num2str(c) ' written)']); end
c=fwrite(fid,hdr.magic,'char');
if (c~=4) error(['hdr.magic must be 4 chars (' num2str(c) ' written)']); end

if (~leave_open)
  fclose(fid);
  fid=[];
end
