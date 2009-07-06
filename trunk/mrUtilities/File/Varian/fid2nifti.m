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
%%
function fid2nifti(fidname,varargin)

% check arguments
if nargin < 1
  help fid2nifti
  return
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
  fid = getfid(fidname);
  if isempty(fid.data)
    disp(sprintf('(fid2nifti) WARNING file %s could not be read',fidname));
    continue
  end

  % get the qform
  qform44 = fid2xform(fidname);

  % create an empty header
  hdr = cbiCreateNiftiHeader(fid.data);

  % set the qform
  hdr = cbiSetNiftiQform(hdr,qform44);
  hdr = cbiSetNiftiSform(hdr,qform44);

  % get the volume TR (framePeriod) for EPI 
  procpar = readprocpar(fidname);
  tr = procpar.tr;
  % if we run mutliple shots, volume TR = slice TR * shots 
  % I do not know if shots equal navechoes or not, but as you read from
  % petable, they should be same
  if isfield(procpar,'navechoes')
    tr = tr * procpar.navechoes;
  end
  % if we run slice at not interleaved way, 
  %  volume TR = slice TR * shots * slice number
  if isfield(procpar,'intlv') && strcmp(procpar.intlv, 'n')
    tr = tr * length(procpar.pss);
  end

  % make a filename 
  niftiname = setext(fixBadChars(stripext(fidname),{'.','_'}),'hdr');

  % write the file
  disp(sprintf('(fid2nifti) Converting %s to %s',fidname,niftiname));
  cbiWriteNifti(niftiname,fid.data,hdr);
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
