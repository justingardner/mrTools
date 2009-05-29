% getfid.m
%
%      usage: getfid()
%         by: justin gardner
%       date: 05/08/03
%    purpose: reads k-space data from fid file
%             and transforms it into image
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
if (verbose),disppercent(-inf,sprintf('Reading %s...',fidname));end
d = getfidk(fidname);
if (verbose),disppercent(inf,sprintf('done.\n',fidname));end

% if it is empty then something has failed
if (isempty(d.data))
  return
end

d.dim = size(d.data);

% everything is ok, then transform data
% everything is ok, then transform data
if(verbose),disppercent(-inf,'transforming data');end
for i=1:size(d.data,3)
  for j=1:size(d.data,4)
    if (verbose)
      disppercent(min(round(100*((i-1)*size(d.data,4)+j)/(size(d.data,3)*size(d.data,4))),99));
    end
    % need to zeropad or not
    if zeropad
      % zeropad the data
      thisdata=zeros(zeropad,zeropad);
      % get this image
      thisdata(1:d.dim(1),1:d.dim(2)) = squeeze(d.data(:,:,i,j));
      % fft 
      data(:,:,i,j) = fftshift(abs(fft2(thisdata)))/(size(thisdata,1)*size(thisdata,2));
    else
      % simply fft data
      d.data(:,:,i,j) = fftshift(abs(fft2(d.data(:,:,i,j))))/(size(d.data,1)*size(d.data,2));
    end
  end
end

% if zeropad fix up some stuff
if zeropad
  d.data = data;
  d.dim(1) = zeropad;d.dim(2) = zeropad;
end
d.zeropad = zeropad;


if (verbose), disppercent(inf); end
