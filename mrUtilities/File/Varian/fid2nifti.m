% fid2nifti.m
%
%        $Id:$ 
%      usage: fid2nifti(fidname)
%         by: justin gardner
%       date: 04/30/09
%    purpose: convert a fid 2D anatomy file into a nifti file. 
%             Preserves slice orientation infomation found in procpar
%             You can also use it with wild cards:
%
%       e.g.: fid2nifti 2d*.fid
% 
%             If you want to just load (not save a converted file)
%
%             [d h] = fid2nifti('fidname.fid');
%
%             You can get verbose information by doing
%
%             [d h] = fid2nifti('fidname.fid',1);
%
function [outdata outhdr] = fid2nifti.m(fidname,varargin)

outhdr = [];

% check arguments
if nargin < 1
  help fid2nifti
  return
end

% if the first argument is a numeric it means to set verbose setting
verbose = 0;
if (length(varargin) >= 1) && isscalar(varargin{1})
  verbose = varargin{1};
  temp = {};
  for i = 2:length(varargin)
    temp{i-1} = varargin{i};
  end
  varargin = temp;
end

% if this has no path, then check for search pattern
if ~iscell(fidname) && strcmp(getLastDir(fidname),fidname)
  % use regexp to serach for matching filenames
  dirlist = dir;
  % make sure fid ends in .fid
  fidname = setext(fidname,'fid',0);

  searchstr = sprintf('^%s$',fixBadChars(fidname,{'*','.*'}));
  fidname = {};
  for i = 1:length(dirlist)
    if regexp(dirlist(i).name,searchstr)
      fidname{end+1} = dirlist(i).name;
    end
  end
end

% cat with other arguments and make sure we have a cell array
if nargin > 1
  fidnames = cellcat(fidname,varargin);
else
  fidnames = cellArray(fidname);
end

% make sure we still have something to do
if isempty(fidnames)
  disp(sprintf('(fid2nifti) No matching files'));
  return
end

for i = 1:length(fidnames)
  % get this fidname and make sure it ends in fid
  fidname = fidnames{i};
  fidname = setext(fidname,'fid',0);
  
  if ~isdir(fidname)
    disp(sprintf('(fid2nifti) WARNING: Could not find file %s',fidname));
    continue
  end
  
  % get the fid
  fid = getfid(fidname,verbose);
  if isempty(fid.data)
    disp(sprintf('(fid2nifti) WARNING file %s could not be read',fidname));
    continue
  end

  % check to see if we have to merge coils
  keyboard
  if size(fid.data,5) > 1
    numReceivers = size(fid.data,5);
    % see if we have to merge coils
    if numReceivers > 1
      if numReceivers ~= fid.dim(end)
	disp(sprintf('(fid2nifti) Num receivers (%i) does not match data dim (%i)',numReceivers,fid.dim(end)));
      end
      % display what we are doing
      if verbose,disppercent(-inf,sprintf('(fid2nifti) Taking sum of squares of %i coils',numReceivers));end
      % merge the coils
      for volNum = 1:size(fid.data,4)
	sumOfSquares = zeros(fid.dim(1:3));
	for receiverNum = 1:numReceivers
	  sumOfSquares = sumOfSquares+fid.data(:,:,:,volNum,receiverNum).^2;
	end
	data(:,:,:,volNum) = sqrt(sumOfSquares);
	if verbose,disppercent(volNum/size(fid.data,4));end
      end
      if verbose,disppercent(inf);end
      fid.data = data;
      fid.dim = size(data);
    end
  end
  
  % create a header
  hdr = fid2niftihdr(fidname);
  
  % read the procpar
  procpar = readprocpar(fidname);
  
  % reorder slices if necessary
  [pss sliceIndex] = sort(procpar.pss);
  if ~isequal(sliceIndex,1:size(fid.data,3))
    fid.data = fid.data(:,:,sliceIndex,:);
  end

  % check the dimensions of the data versus the dimensions in the header
  if ~isequal([size(fid.data,1) size(fid.data,2) size(fid.data,3)],hdr.dim(2:4)')
    disp(sprintf('(fid2nifti) Header info from procpar does not match size %s with data read %s',mynum2str(hdr.dim(2:4)),mynum2str(size(fid.data))));
  end
  
  % make a filename 
  niftiname = setext(fixBadChars(stripext(fidname),{'.','_'}),'hdr');

  % write the file, but only if we aren't taking an output argument
  if nargout == 0
    disp(sprintf('(fid2nifti) Converting %s to %s',fidname,niftiname));
    cbiWriteNifti(niftiname,fid.data,hdr);
  else
    outhdr{i} = hdr;
    outdata{i} = fid.data;
  end
end

% return single header if that is all that is asked for
if length(outhdr) == 1,
  outdata = outdata{1};
  outhdr = outhdr{1};
end


% cellcat.m
%
%      usage: cellcat()
%         by: justin gardner
%       date: 06/24/07
%    purpose: cat two cell arrays. this preserves order and keeps duplicates (unlike union)
%
function c = cellcat(c1,c2)

% check arguments
if ~any(nargin == [2])
  c={};
  help cellcat
  return
end

% make sure that c1 and c2 are cell arrays
c1 = cellArray(c1);
c2 = cellArray(c2);

c = c1;
for i = 1:length(c2)
  c{end+1} = c2{i};
end
