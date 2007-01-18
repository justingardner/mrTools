function [byteswritten,hdr]=cbiWriteNifti(fname,data,hdr,prec,subset,short_nan);
% [byteswritten,hdr]=cbiWriteNifti(filename,data,hdr,prec,subset,short_nan) 
%  Uses user-defined header struct
% [byteswritten,hdr]=cbiWriteNifti(filename,data,[],prec);
%  Creates a new header
% 
% The header is always passed through cbiCreateNiftiHeader to ensure 
% consistency with the data.
%
%  prec:      Precision (dataype) of output. Default (if no header is specified) is 'float32'. Use [] for default.
%             Should be a recognized Matlab (or nifti) data type string.
%  subset:    4x1 cell array describing image subset to save. 1-offset (Matlab-style).
%             Only the following options are supported:
%             - to save a single z-slice (e.g. 4):  subset={[],[],4,[]}
%             - to save a single/multiple volumes of a time series (e.g. volumes 2-9):  subset={[],[],[],[2 9]}
%  short_nan: NaN handling for signed short (int16) data. If 1, will treat save NaN's as -32768 (smallest 
%             representable number) reserving -32767..32767 for scaled data; otherwise will save NaN as 0.
%             Default is 1 (use -32768 as NaN).
  
  

[pathstr,bname,ext,ver]=fileparts(fname);
  
switch (ext)
 case '.nii'  
  hdr.singleFile=1;
  hdr.hdr_name=fname;
  hdr.img_name=fname;
  hdr.magic=sprintf('%s\0','n+1');
 case {'.hdr','.img'}
  hdr.single_file=0;
  hdr.hdr_name=fullfile(pathstr,[bname '.hdr']);
  hdr.img_name=fullfile(pathstr,[bname '.img']);
  hdr.magic=sprintf('%s\0','ni1');
 case '.gz','.Z' % zipped
  error('No support for zipped NIFTI-1 format under Matlab.');
 otherwise
  error('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img');
end

% Ensure header matches data
if (exist('prec') & ~isempty(prec))
  hdr=cbiCreateNiftiHeader(hdr,data,'matlab_datatype',prec);
else
  hdr=cbiCreateNiftiHeader(hdr,data);
end

% Write header
no_overwrite=0;
[hdr,fid]=cbiWriteNiftiHeader(hdr,fname,no_overwrite,hdr.single_file);

% Prepare to write data
if (~hdr.single_file)
  fid=fopen(hdr.img_name,'wb',hdr.endian);
end

if (~exist('subset'))
  subset={[],[],[],[]};
end
if (~exist('short_nan'))
  short_nan=1;
end

headerdim=hdr.dim(2:5); % Matlab 1-offset - hdr.dim(1) is actually hdr.dim(0)
headerdim(headerdim==0)=1; % Force null dimensions to be 1
is5D=0;
if (hdr.dim(6)>1)
  if (hdr.dim(5)>1)
    error('No support for 5D data with multiple time points!');
  end
  is5D=1;
  headerdim(4)=hdr.dim(6);
  disp('5D data set detected');
end

loadSize=zeros(4,1);
for n=1:4
  if (length(subset{n})==1)
    subset{n}=[subset{n} subset{n}];    
  elseif (length(subset{n})>2)
    error('subset should be a scalar or 2-vector');
  end
  if (isempty(subset{n}))
    loadSize(n)=headerdim(n);
    subset{n}=[1 loadSize(n)];
  else
    loadSize(n)=subset{n}(2)-subset{n}(1)+1;  
  end
end

if (any(loadSize>headerdim(1:4)))
  error('subset index larger than image dimensions!');
elseif (any(loadSize(1:2)<headerdim(1:2)))
  error('no support for saving subvolumes of data; only entire z-slices may be saved.');
end

% Move to beginning of data
fseek(fid,hdr.vox_offset,'bof');

dataSize=prod(loadSize);
writeFormat=hdr.matlab_datatype;
switch (hdr.matlab_datatype)    
 case 'binary'
  error('No support for binary data');
 case 'complex64'
  writeFormat=float32;
 case 'complex128'
  writeFormat=float64;
 case 'RGB'
  error('No support for RGB data');
 case {'complex256','float128'}
  error('No support for 128-bit data on this platform!');
end

bytesPerElement=cbiSizeofNifti(writeFormat);
% For each time point:

% Move to correct location
readOrigin=sub2ind(headerdim(1:4)',subset{1}(1),subset{2}(1),subset{3}(1),subset{4}(1))-1; % now we're in C-land, hence 0-offset
											   
% Elements to read every time point
volSize=prod(headerdim(1:3));
readSize=prod(loadSize(1:3));
% Difference between volSize and readSize => offset to seek after every read
readOffset=volSize-readSize;
% Current position in data array
currPos=1;
if (strfind(hdr.matlab_datatype,'complex'))
  % Every voxel corresponds to two elements in the file
  readOrigin=readOrigin*2;
  readSize=readSize*2;
  readOffset=readOffset*2;
end
% Position file at first voxel
byteswritten=0;
fseek(fid,readOrigin*bytesPerElement,'cof');
for t=subset{4}(1):subset{4}(2)  
  % Extract and convert current subset of data
  saveData=convertData(data(currPos:currPos+readSize-1),hdr,short_nan);
  if (strfind(hdr.matlab_datatype,'complex'))
    % Separate complex data into real and imaginary parts
    realData=real(saveData);
    imagData=imag(saveData);
    saveData=zeros(2*prod(size(saveData)),1);
    saveData(1:2:readSize/2-1)=realData;
    saveData(2:2:readSize/2)=imagData;
  end

  % Write converted subset to file
  count=fwrite(fid,saveData,writeFormat);
  byteswritten=byteswritten+count;
  if (count~=readSize)
    fclose(fid);
    error(['Error writing to file ' hdr.img_name]);
  end
  if (readOffset>0)
    % fseek to next time point
    fseek(fid,readOffset*bytesPerElement,'cof');
  end
  currPos=currPos+readSize;
end

fclose(fid);

return


function data=convertData(data,hdr,short_nan);
% Scales and shifts data (using hdr.scl_slope and hdr.scl_inter)
% and changes NaN's to 0 or MAXINT for non-floating point formats
%
  % If data is the correct class and no scaling needed, just return
  if (strcmp(class(data),hdr.matlab_datatype) ... 
      & (isempty(hdr.scl_slope) | hdr.scl_slope==1) ...
      & (isempty(hdr.scl_inter) | hdr.scl_inter==0))
    return;
  end

  % Scale and shift data
  if (hdr.scl_slope~=1 | hdr.scl_inter~=0)
    data=double(data);
    data=(data-hdr.scl_inter)./hdr.scl_slope;
  end
    
  % Change NaNs for non-floating point datatypes
  switch (hdr.matlab_datatype)    
   case {'binary','uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
    data(isnan(data))=0;
   case 'int16'
    if (short_nan)
      data(isnan(data))=-32768;
    else
      data(isnan(data))=0;
    end
  end

