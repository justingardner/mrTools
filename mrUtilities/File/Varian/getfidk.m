% getfidk.m
%
%        $Id:$ 
%      usage: getfidk(filename)
%         by: justin gardner
%       date: 08/10/11
%    purpose: returns raw fid data using getfidkraw.c. Reorders lines of k-space appropriately
%             so that images can be fourier transformed
%
function retval = getfidk(filename,varargin)

% check arguments
if nargin < 1
  help getfidk
  return
end

% parse arguments
verbose = 0;
getArgs(varargin,{'verbose=1'});

% check filename
if ~isdir(filename)
  filename = setext(filename,'fid');
  if ~isdir(filename)
    disp(sprintf('(getfidk) Could not find fid directory %s',filename));
    return
  end
end

% load the infromation from the procpar using fid2xform
[xform info] = fid2xform(filename);
if isempty(info),return,end

% look for fid file
fidFilename = fullfile(filename,'fid');
if ~isfile(fidFilename)
  disp(sprintf('(getfidk) Could not find fid file: %s', fidFilename));
  return
end

% load the fid file using getfidkraw
d = getfidkraw(fidFilename,verbose);

% make sure the dimensions match
if d.nblocks ~= prod(info.dim(2:4))
  disp(sprintf('(getfidk) Number of lines of k-space found in fid (%i) does not match expected number %ix%ix%i=%i',d.nblocks,info.dim(2),info.dim(3),info.dim(4),prod(info.dim(2:4))));
  return
end

% the data are stored such that we read one k-space line
% for each image, cyclying through all slices and volumes
% then we read the next k-space line etc.
kNum = 1;clear i;
for volNum = 1:info.dim(4)
  for sliceNum = 1:info.dim(3)
    for kLine = 1:info.dim(2)
      d.data(kLine,:,sliceNum,volNum) = d.real(kNum,:) + i*d.imag(kNum,:);
      kNum = kNum+1;
    end
  end
end
d2 = getfidk(filename);
keyboard

