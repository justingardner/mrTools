% getfid.m
%
%      usage: getfid()
%         by: justin gardner
%       date: 05/08/03
%    purpose: reads k-space data from fid file
%             and transforms it into image. Dimensions are
%             set to x,y,slice,vol,coil (note that vol and coil
%             dimensions are swapped relative to what getfidk returns)
%
%             Also for 2D images with multiple receivers, note that Varian
%             encodes the image with the slices and receivers swapped. This
%             program (but not getfidk) fixes that swap.
%
function d = getfid(fidname,verbose,zeropad)

% check input arguments
if (nargin == 1)
  verbose = 0;
elseif ~ any(nargin == [2 3]) 
  help getfid;
  return
end

% default to no zeropad
if exist('zeropad') ~= 1,zeropad = 0;,end

t0 = clock;

% read the k-space data from the fid
if (verbose),disppercent(-inf,sprintf('(getfid) Reading %s...',fidname));end
d = getfidk(fidname,verbose);
if (verbose),disppercent(inf,sprintf('done.\n',fidname));end
% if it is empty then something has failed
if (isempty(d.data))
  return
end

d.dim = size(d.data);

% get fidinfo
[xform info] = fid2xform(fidname);

if info.receiversAndSlicesSwapped
  if size(d.data,4) > 1
    if verbose, disp(sprintf('(getfid) Swapping receivers and slices'));end
    d.data = reshape(d.data,size(d.data,1),size(d.data,2),size(d.data,4),size(d.data,3),size(d.data,5));
    d.data = permute(d.data,[1 2 4 3 5]);
  end
end

% everything is ok, then transform data
if(verbose),disppercent(-inf,'(getfid) Transforming data');end
for i = 1:size(d.data,3)
  for j = 1:size(d.data,4)
    for k = 1:size(d.data,5)
      % need to zeropad or not
      if zeropad
	% zeropad the data
	thisdata=zeros(zeropad,zeropad);
	% get this image
	thisdata(1:d.dim(1),1:d.dim(2)) = squeeze(d.data(:,:,i,j,k));
	% fft 
	data(:,:,i,k,j) = fftshift(abs(fft2(thisdata)))/(size(thisdata,1)*size(thisdata,2));
      else
	% simply fft data
	data(:,:,i,k,j) = fftshift(abs(fft2(d.data(:,:,i,j,k))))/(size(d.data,1)*size(d.data,2));
      end
    end
    if (verbose)
      disppercent(((i-1)*size(d.data,4)+j)/(size(d.data,3)*size(d.data,4)));
    end
  end
end

d.data = data;
d.dim = size(data);
% if zeropad fix up some stuff
if zeropad
  d.dim(1) = zeropad;d.dim(2) = zeropad;
end
d.zeropad = zeropad;


if (verbose), disppercent(inf); end
