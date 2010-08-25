% getfid.m
%
%      usage: getfid(fidname,<verbose>,<zeropad>,<movepro>,<kspace>)
%         by: justin gardner
%       date: 05/08/03
%    purpose: reads k-space data from fid file and transforms into an image
%             and transforms it into image. Dimensions are
%             set to x,y,slice,vol,coil (note that vol and coil
%             dimensions are swapped relative to what getfidk returns)
% 
%             optional arguments
%             verbose=0 to 1 for verbose info
%             zeropad=0 set to the size you want to zeropad out to (e.g. 256)
%             movepro=0 set to how much you want to movepro (default=0)
%
%             Also for 2D images with multiple receivers, note that Varian
%             encodes the image with the slices and receivers swapped. This
%             program (but not getfidk) fixes that swap.
%
%             If the kspace argument is set to 1, then this will return the
%             k-space data (this can be useful if you want to zeropad, moverpro etc.
%             which getfidk does not do.
%
function d = getfid(fidname,verbose,zeropad,movepro,kspace)

% check input arguments
if (nargin == 1)
  verbose = 0;
elseif ~ any(nargin == [2 3 4 5]) 
  help getfid;
  return
end

% default to no zeropad
if ieNotDefined('zeropad'),zeropad = 0;end
if ieNotDefined('movepro'),movepro = 0;end
if ieNotDefined('kspace'),kspace = 0;end
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
[xform info] = fid2xform(fidname,-1);

if info.receiversAndSlicesSwapped
  if size(d.data,4) > 1
    if verbose, disp(sprintf('(getfid) Swapping receivers and slices'));end
    d.data = reshape(d.data,size(d.data,1),size(d.data,2),size(d.data,4),size(d.data,3),size(d.data,5));
    d.data = permute(d.data,[1 2 4 3 5]);
  end
end

% find the needed phase shift to add to the image if we need to movepro
if movepro ~= 0
  % figure out how much we have to shift
  proshift = movepro/(info.procpar.lro/d.dim(1));
  % and create the shift
  phaseshift = proshift*(0:2*pi./d.dim(1):2*pi)';
  phaseshift = phaseshift(2:end);
  phaseshift = phaseshift*ones(1,d.dim(2));
  if (verbose),disp(sprintf('(getfid) Shifting pro by: %f',movepro));end
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
	if movepro == 0
	  % get this image
	  thisdata(1:d.dim(1),1:d.dim(2)) = squeeze(d.data(:,:,i,j,k));
	else
	  % get this image with pro shift
	  thisdata(1:d.dim(1),1:d.dim(2)) = squeeze(d.data(:,:,i,j,k)).*exp(sqrt(-1)*phaseshift);
	end
	% fft 
	data(:,:,i,k,j) = myfft(thisdata,kspace);
      else
	% simply fft data, moving pro if necessary
	if movepro == 0
	  data(:,:,i,k,j) = myfft(d.data(:,:,i,j,k),kspace);
	else
	  data(:,:,i,k,j) = myfft(d.data(:,:,i,j,k).*exp(sqrt(-1)*phaseshift),kspace);
	end
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

d.info = info;

if (verbose), disppercent(inf); end

%%%%%%%%%%%%%%%
%    myfft    %
%%%%%%%%%%%%%%%
function data = myfft(data,kspace)

% this just simply takes the 2D fft, shifts and gets the real part of the data.
% if kspace is set, it does nothing
if ~kspace
  data = fftshift(abs(fft2(data)))/(size(data,1)*size(data,2));
end