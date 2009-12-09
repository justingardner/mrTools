% fid2niftihdr.m
%
%        $Id:$ 
%      usage: fid2niftihdr(fidname,<verbose>)
%         by: justin gardner
%       date: 12/09/09
%    purpose: strips out of fid2nifti the functionality to create nifti headers from procpar
%
function hdr = fid2niftihdr(fidname,verbose)

hdr = [];
if ieNotDefined('verbose'),verbose=1;end

% check arguments
if ~any(nargin == [1])
  help fid2niftihdr
  return
end

% create an empty header
hdr = cbiCreateNiftiHeader;

% read procpar
procpar = readprocpar(fidname);
if isempty(procpar),disp(sprintf('(fid2niftihdr) Could not find procpar in %s',fidname)),return,end

% get the number of navechoes
if isfield(procpar,'navechoes')
  navechoes = procpar.navechoes;
else
  navechoes = 0;
end

% get the dimensions of the scan
dim = [procpar.ni procpar.nv length(procpar.pss)];

% check to see if this is a sense reconstruction and what the sense acceleration factor is
senseFactor = [];
if isfield(procpar,'petable') && ~isempty(procpar.petable{1}) && (length(procpar.petable{1}>1))
  % petable name that ends in r means a sense protocol
  if isequal(procpar.petable{1}(end),'r')
    % check for the sense factor, (hopefully this number won't be greater than 9!!)
    senseFactor = str2num(procpar.petable{1}(end-1));
    % fix dimensions
    dim(1) = dim(1)*senseFactor;
    dim(2) = dim(2)*senseFactor;
    % print message
    if verbose>=0, disp(sprintf('(fid2niftihdr) Found a sense factor of %i. Dims are now [%i %i]',senseFactor,dim(1),dim(2)-navechoes));end
  end
end

% remove navigator echoes from k-space
dim(2) = dim(2) - navechoes;

% now set dimensions in header
hdr.dim(1:4) = [3 dim];

% get the qform
qform44 = fid2xform(fidname);

% set the qform
hdr = cbiSetNiftiQform(hdr,qform44);

% set voxel dimensions (note that this is done after cbiSetNiftiQform, since that
% function decides that the voxel spacing and the voxel size should be the same).
hdr.pixdim(2:4) = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) procpar.thk];

% get the volume TR (framePeriod) for EPI 
tr = procpar.tr;
% if we run mutliple shots, volume TR = slice TR * shots 
if isfield(procpar,'navechoes')
  tr = tr * procpar.numshots;
end
% if we run slice at not interleaved way, 
%  volume TR = slice TR * shots * slice number
if isfield(procpar,'intlv') && strcmp(procpar.intlv, 'n')
  tr = tr * length(procpar.pss);
end
hdr.pixdim(hdr.dim(1)+2) = tr*1000;



