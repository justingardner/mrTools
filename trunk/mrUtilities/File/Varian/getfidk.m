% getfidk.m
%
%        $Id:$ 
%      usage: getfidk(filename)
%         by: justin gardner
%       date: 08/10/11
%    purpose: returns raw fid data using getfidkraw.c. Reorders lines of k-space appropriately
%             so that images can be fourier transformed
%
function d = getfidk(filename,varargin)

d.data = [];
% check arguments
if nargin < 1
  help getfidk
  return
end

% parse arguments
verbose = [];
getArgs(varargin,{'verbose=0'});

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

% take transpose to make these easier to deal with
d.real = d.real';
d.imag = d.imag';

% compute ow any volumes, slices, phase encode lines and receivers we have
numVolumes = info.dim(4)+info.nRefVolumes;
numSlices = info.dim(3);
numPhaseEncodeLines = info.dim(1)/info.accFactor;
numReceivers = info.numReceivers;

% now compute how many lines we need to read the data from
if info.compressedFid
  % for compressedFid each line of data has all phase encode lines
  numLines = numSlices*numVolumes*numReceivers;
  % we will also need to know the line order as returned form petable
  if ~isfield(info.procpar,'pelist') || length(info.procpar.pelist) > 1
    lineorder = info.procpar.pelist-min(info.procpar.pelist)+1;
  else
    % pelist variable not set, read petable
    petable = readpetable(info.procpar.petable{1},'verbose',verbose);
    % petable not found, give up
    if isempty(petable)
      disp(sprintf('(getfidk) !!!! Unknown petable !!!!'));
      keyboard
    end
    % get the line ordering from petable
    lineorder = petable.t1';
    lineorder = lineorder-min(lineorder)+1;
  end
else
  numLines = numPhaseEncodeLines*numSlices*numVolumes*numReceivers;
end

% make sure the dimensions match
if d.nblocks ~= numLines
  disp(sprintf('(getfidk) !!! Number of lines of k-space found in fid (%i) does not match expected number %i (%ix%ix%ix%i %i receivers accFactor %i)',d.nblocks,numLines,info.dim(1),info.dim(2),info.dim(3),info.dim(4),info.numReceivers,info.accFactor));
  d.data = [];
  return
end

if verbose
  disp(sprintf('(getfidk) Numlines: %i (%ix%ix%ix%i %i receivers accFactor %i)',numLines,info.dim(1),info.dim(2),info.dim(3),info.dim(4),info.numReceivers,info.accFactor));
end

% read the data from fid block structure
kNum = 1;clear i;
d.data = nan(info.dim(2),numPhaseEncodeLines,numSlices,numReceivers,numVolumes);
if verbose,disppercent(-inf,'(getfidk) Reordering data');end
for sliceNum = 1:numSlices
  for receiverNum = 1:numReceivers
    for volNum = 1:numVolumes
	% compressed fids have all lines of k-space in one single block of data
      if info.compressedFid
	% note conversion here to double
	d.data(:,lineorder,sliceNum,receiverNum,volNum) = reshape(double(d.real(kNum,:)) + i*double(d.imag(kNum,:)),numPhaseEncodeLines,info.dim(2));
	kNum = kNum+1;
      else
	% uncompress fid contains one line of k-space per block
	for kLine = 1:numPhaseEncodeLines
	  % note conversion here to double
	  d.data(:,kLine,sliceNum,receiverNum,volNum) = double(d.real(kNum,:)) + i*double(d.imag(kNum,:));
	  kNum = kNum+1;
	end
      end
    end
    if verbose,disppercent(calcPercentDone(volNum,numVolumes,receiverNum,numReceivers));end
  end
end
if verbose,disppercent(inf);end

