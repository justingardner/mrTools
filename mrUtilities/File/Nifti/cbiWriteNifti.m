function [byteswritten,hdr]=cbiWriteNifti(fname,data,hdr,prec,subset,short_nan,verbose)
% [byteswritten,hdr]=cbiWriteNifti(filename,data,hdr,prec,subset,short_nan,verbose) 
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
%             - to save a single z-slice, e.g. 4:  subset={[],[],4,[]}
%             - to save a single/multiple volumes of a time series, e.g. volumes 2-9:  subset={[],[],[],[2 9]}
%             - to append a single/multiple volumes to a time series, e.g. append 3 volumes to a 9-volume file: subset={[],[],[],[10 12]}
%  short_nan: NaN handling for signed short (int16) data. If 1, will treat save NaN's as -32768 (smallest 
%             representable number) reserving -32767..32767 for scaled data; otherwise will save NaN as 0.
%             Default is 1 (use -32768 as NaN).
% $Id$
  
if nargin == 0
  help cbiWriteNifti;
  return
end
if (ieNotDefined('subset'))
  subset={[],[],[],[]};
end
if (ieNotDefined('short_nan'))
  short_nan=1;
end
if ieNotDefined('verbose')
  verbose=1;
end

[pathstr,bname,ext]=fileparts(fname);

% 1D curvature file fix: 
% check to see if this is 1D image that has a large second
% dimension, if so, we need to reorder the data (this is 
% necessary for surface data in which we *want* to store
% a 1xnumVertices image with the surface coloring data,
% but nifti's idiotic two byte width limitation prevents
% us from doing so). Note that there is a complimentary 
% piece of code in cbiReadNifti that handles unpacking
% this data. Note that the realWidth is saved in the
% header as part of the descripiton (i.e. as 'realWidth: rest of description string')
% -j.
maxWidth = 2^15-1;realWidth = [];
if (size(data,1) == 1) && (size(data,2) > maxWidth)
  realWidth = size(data,2);
  % find out how many rows we need to store the data
  numRows = ceil(realWidth/maxWidth);
  % fill the data out with nans past the real data
  data(1,end+1:numRows*maxWidth) = nan;
  % and reshape
  data = reshape(data,numRows,maxWidth);
  % print out what happened
  if verbose, disp(sprintf('(cbiWriteNifti) Saving a 1x%i image, reshaping to %ix%i to fit nifti limitation of max width of %i',realWidth,numRows,maxWidth,maxWidth));end;
end
  
if isstruct(data)
  disp(sprintf('(cbiWriteNifti) UHOH! Data appears to be a structure rather than a matrix.'));
  return
end

switch (ext)
 case '.nii'  
  hdr.single_file=1;
  hdr.hdr_name=fname;
  hdr.img_name=fname;
  hdr.magic=sprintf('%s\0','n+1');
  hdr.vox_offset=352;
 case {'.hdr','.img'}
  hdr.single_file=0;
  hdr.hdr_name=fullfile(pathstr,[bname '.hdr']);
  hdr.img_name=fullfile(pathstr,[bname '.img']);
  hdr.magic=sprintf('%s\0','ni1');
  hdr.vox_offset=0; % important!
 case '.gz','.Z' % zipped
  mrErrorDlg('No support for zipped NIFTI-1 format under Matlab.');
 otherwise
  mrErrorDlg('Not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img');
end

% deal with subset and append options
[dataSize(1),dataSize(2),dataSize(3),dataSize(4)] = size(data); %size of the new data to save
subsetIndices = zeros(4,2);
for n=1:4
  if (length(subset{n})==1)
    subsetIndices(n,:)=[subset{n} subset{n}];
  elseif (length(subset{n})>2)
    mrErrorDlg('subset should be a scalar or 2-vector');
  elseif (isempty(subset{n}))
    subsetIndices(n,:)=[1 dataSize(n)];
  else
    subsetIndices(n,:) = subset{n};
  end
end

if ~isequal(subsetIndices,[ones(4,1) dataSize']) % if writing a subset or appending volumes, need to check header of existing file
  hdrDestination = cbiReadNiftiHeader(fname);
  if isempty(hdrDestination)
    mrErrorDlg(sprintf('(cbiWriteNifti) Cannot read %s',fname));
  end
  hdrDestDim =  hdrDestination.dim(2:5);
  hdrDestDim(hdrDestDim==0)=1;
  if ~isequal(subsetIndices(1:2,:),[ones(2,1) hdrDestDim(1:2)])
    mrErrorDlg('no support for saving subvolumes of data on x and y dimensions; only entire z-slices or volumes may be saved.');
  elseif subsetIndices(3,2)>hdrDestDim(3)
    mrErrorDlg('z subset index larger than file image dimensions!');
  elseif subsetIndices(4,2)>hdrDestDim(4) && subsetIndices(4,1)-hdrDestDim(4)~=1
    mrErrorDlg('(cbiWriteNifti) When appending data on the 4th dimension, the first frame index must be exactly adjacent to the last frame in the file');
  end
end

% Ensure header matches data
hdr.dim(2:5) = max(dataSize',subsetIndices(:,2));
if ~ieNotDefined('prec')
  hdr=cbiCreateNiftiHeader(hdr,'matlab_datatype',prec);
else
  hdr=cbiCreateNiftiHeader(hdr);
end
if (~strcmp(class(data),hdr.matlab_datatype))
  if verbose, disp(['(cbiWriteNifti) Scaling data from ' class(data) ' to ' hdr.matlab_datatype]);end
%  disp('To avoid this, cast data to desired format before calling cbiWriteNifti, e.g.')
%  disp('cbiWriteNifti(''myfilename'',int16(data),hdr,''int16'')')
end

% get hdr scaling factor and convert data if necessary
[data,hdr]=convertData(data,hdr,short_nan);  

% 1D curvature file fix: squirrel away the real width of the image, if we had to
% convert for a 1D image.
if ~isempty(realWidth )
  % write the real with into the description
  hdr.descrip = sprintf('%i:%s',realWidth,hdr.descrip);
  % clip the length of the description to the max of 80
  hdr.descrip = hdr.descrip(1:min(end,80));
  hdr.descrip(end+1:80) = 0;
end

% Write header
no_overwrite=0;
[hdr,fid]=cbiWriteNiftiHeader(hdr,fname,no_overwrite,hdr.single_file);

headerdim=hdr.dim(2:5); % Matlab 1-offset - hdr.dim(1) is actually hdr.dim(0)
headerdim(headerdim==0)=1; % Force null dimensions to be 1
if (hdr.dim(6)>1)
  if (hdr.dim(5)>1)
    mrErrorDlg('No support for 5D data with multiple time points!');
  end
  headerdim(4)=hdr.dim(6);
  if verbose, disp('5D data set detected');end
end

writeFormat=hdr.matlab_datatype;
switch (hdr.matlab_datatype)    
 case 'binary'
  mrErrorDlg('No support for binary data')
 case 'complex64'
  writeFormat=float32;
 case 'complex128'
  writeFormat=float64;
 case 'RGB'
  mrErrorDlg('No support for RGB data');
 case {'complex256','float128'}
  mrErrorDlg('No support for 128-bit data on this platform!');
end
bytesPerElement=cbiSizeofNifti(writeFormat);

% Prepare to write data
if (~hdr.single_file)
  if subsetIndices(4,1)>headerdim(4)
    permission = 'a';
  else
    permission = 'wb';
  end
  fid=fopen(hdr.img_name,permission,hdr.endian);
  if fid == -1,mrErrorDlg(sprintf('(cbiWriteNiftiHeader) Could not open file %s',fname));end
end

try
  %  Write emptiness so that we can move to the right offset. (This library does not currently support extensions, which would otherwise go here)
  % 11/10/05 PJ
  if ftell(fid) < hdr.vox_offset
    %  fwrite(fid,0,sprintf('integer*%d',(hdr.vox_offset-ftell(fid))));
    c=hdr.vox_offset-ftell(fid);
    if (fwrite(fid,zeros(c,1),'uint8')~=c)
      mrErrorDlg('error writing extension padding')
    end
  end

  % Move to beginning of data
  status = fseek(fid,hdr.vox_offset,'bof');
  if status
    ferror(fid)
  end

  % Move to correct location
  readOrigin=sub2ind([headerdim(1:3)' max(headerdim(4),subsetIndices(4,2))],subsetIndices(1,1),subsetIndices(2,1),subsetIndices(3,1),subsetIndices(4,1))-1; % now we're in C-land, hence 0-offset

  % Elements to write every time point
  volSize=prod(headerdim(1:3));
  writeSize=prod(diff(subsetIndices(1:3,:),1,2)+1);
  % Difference between volSize and readSize => offset to seek after every write
  readOffset=volSize-writeSize;
  % Current position in data array
  currPos=1;
  if (strfind(hdr.matlab_datatype,'complex'))
    % Every voxel corresponds to two elements in the file
    readOrigin=readOrigin*2;
    writeSize=writeSize*2;
    readOffset=readOffset*2;
  end
  % Position file at first voxel
  byteswritten=0;
  fseek(fid,readOrigin*bytesPerElement,'cof');
  for t=subsetIndices(4,1):subsetIndices(4,2)
    % Extract current subset of data
    saveData=data(currPos:currPos+writeSize-1);
    if (strfind(hdr.matlab_datatype,'complex'))
      % Separate complex data into real and imaginary parts
      realData=real(saveData);
      imagData=imag(saveData);
      saveData=zeros(2*numel(saveData),1);
      saveData(1:2:writeSize/2-1)=realData;
      saveData(2:2:writeSize/2)=imagData;
    end

    % Write converted subset to file
    count=fwrite(fid,saveData,writeFormat);
    byteswritten=byteswritten+count;
    if (count~=writeSize)
      fclose(fid);
      mrErrorDlg(['Error writing to file ' hdr.img_name]);
    end
    if (readOffset>0)
      % fseek to next time point
      fseek(fid,readOffset*bytesPerElement,'cof');
    end
    currPos=currPos+writeSize;
  end
catch exception
  %if anything went wrong, we still want to close the file
  fclose(fid);
  rethrow(exception);
end
fclose(fid);


return
  
function [data,hdr]=convertData(data,hdr,short_nan)
% Scales and shifts data (using hdr.scl_slope and hdr.scl_inter)
% and changes NaN's to 0 or MAXINT for non-floating point formats
% Returns hdr with scale factor changed (bug fixed 20060824)

% Calculate scale factor for non-floating point data
  switch (hdr.matlab_datatype)    
   case 'binary'
    mrErrorDlg('unsupported format')
   case 'uint8'
    MAXINT=2^8-1;
   case {'uint16','ushort'}
    MAXINT=2^16-1;
   case {'uint32','uint'}
    MAXINT=2^32-1;
   case 'uint64'
    MAXINT=2^64-1;
   case 'int8'
    MAXINT=2^7-1;
   case {'int16','short'}
    MAXINT=2^15-1;
   case {'int32','int'}
    MAXINT=2^31-1;
   case 'int64'
    MAXINT=2^63-1;
   otherwise
    MAXINT=0;
  end
  if (MAXINT)
    hdr.scl_slope=max(data(:))/MAXINT;
  end
  
  % Scale and shift data if scale factor is nonzero
  if (~isnan(hdr.scl_slope) && hdr.scl_slope~=0)
    if (hdr.scl_slope~=1 || hdr.scl_inter~=0)
      data=double(data);
      data=(data-hdr.scl_inter)./hdr.scl_slope;
      switch (hdr.matlab_datatype)    
       case {'binary','uint8','uint16','short','ushort','uint32','uint','int','uint64','int8','int16','int32','int64'}
	data=round(data);
       otherwise
	% nothing
      end
    end
  end
    
  % Change NaNs for non-floating point datatypes
  switch (hdr.matlab_datatype)    
   case {'binary','uint8','uint16','ushort','uint32','uint','int','uint64','int8','int16','int32','int64'}
    data(isnan(data))=0;
   case {'int16','short'}
    if (short_nan)
      data(isnan(data))=-32768;
    else
      data(isnan(data))=0;
    end
  end

