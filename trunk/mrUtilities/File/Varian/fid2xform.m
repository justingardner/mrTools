% fid2xform.m
%
%        $Id:$ 
%      usage: [xform info] = fid2xform(fidname,<verbose>)
%         by: justin gardner
%       date: 05/05/09
%    purpose: Convert fields out of a named procpar into a rotation matrix
%             suitable for the nifti header. (fidname is either the name of 
%             a fid, or the structure returned by readprocpar). info is 
%             a sturcture that contains info about the scan like pixel dims
%             sense factors etc.
%
function [xform info] = fid2xform(procpar,verbose,varargin)

xform = [];

% check arguments
if ~any(nargin == [1 2 3])
  help fid2xform
  return
end

% get extra arguments
movepro=[];
getArgs(varargin,{'movepro=0'});

if ieNotDefined('verbose'),verbose = 0;end

% get the procpar
if isstr(procpar)
  procpar = readprocpar(procpar);
  if isempty(procpar),return,end;
elseif ~isstruct(procpar)
  help fid2xform;
  return
end

% move pro if called for
if movepro ~= 0
  procpar.pro = procpar.pro + movepro;
end

% get the number of navechoes
if isfield(procpar,'navechoes')
  navechoes = procpar.navechoes;
else
  navechoes = 0;
end

% check to see if this is an epi
info.isepi = 0;
if (isfield(procpar,'petable'))
  token = procpar.petable{1};
  % the petable name should be something like
  % "epi132alt8k". We want the second number
  if (strncmp(token,'epi',3)) 
    info.isepi = 1;
  end
end


% get the dimensions of the scan
% (procpar.ni is lines of k-space, procpar.nv is number of lines collected including navigator echoes)
% used to use procpar.nv and correct for navechoes, but seems more sensible to just use procpar.ni)
dim = [procpar.np/2 procpar.ni length(procpar.pss)];

% remove navigator echoes from k-space 
%dim(2) = dim(2) - navechoes;

% check to see if this is a sense reconstruction and what the sense acceleration factor is
if procpar.accfactor > 1
  % fix dimensions
  dim(2) = dim(2)*procpar.accfactor;
  % print message
  if verbose>0
    disp(sprintf('(fid2xform) Found a sense factor of %i. Dims are now [%i %i]',procpar.accfactor,dim(1),dim(2)));
  end
end

% make the rotation matrix from the procpar angles
rotmat = euler2rotmatrix(procpar.psi,-procpar.theta,-procpar.phi);

info.processed = 1;
if dim(2) == 0
 info.processed = 0;
 if verbose > 0
   disp(sprintf('(fid2xform) Unprocessed fid directory. Has dim(2)=0'));
 end
 voxsize = [nan nan nan 1];
 voxspacing = [nan nan nan 1];
else
  % get voxel sizes and put in diagonal to multiply against rotmat
  % so that we get the proper spacing
  voxsize = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) procpar.thk 1];

  % vox spacing can be *different* from voxsize, if you skip in your pss
  voxspacing = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) 10*median(diff(sort(procpar.pss))) 1];
  
  % if we have only one slice, then get vox spacing from thk (since there isn't really any vox spacing
  if length(procpar.pss) == 1
    voxspacing(3) = procpar.thk;
  end
end

% check for 3d acquisition
info.acq3d = 0;
if procpar.nv2 > 1
  if verbose>0,disp(sprintf('(fid2xform) 3D acquisition'));end
  % since the 3rd dimension is taken as a single slice with multiple
  % phase encodes, we have to get the voxel size and dimensions differently
  voxsize(3) = 10*procpar.lpe2/procpar.nv2;
  voxspacing(3) = 10*procpar.lpe2/procpar.nv2;
  dim(3) = procpar.nv2;
  % keep in structure that this is a 3d acquisition
  info.acq3d = 1;
  % check to see if it has been processed or not (i.e. pss should be of
  % correct length
  if length(procpar.pss) ~= procpar.nv2
      info.processed = 0;
  end
end

% Now get the offset in mm from the center of the bore that the center of the
% volume is. We can not change the phase encode center. Note that dimensions
% are given in cm, so we must convert to mm.
% 
% Note here about 3D images. The pss is in "reverse" order in the sense that it
% goes from positive numbers to negative numbers. This is fixed in fid2nifit
% since fid2nifti sorts the pss and reorders the slices accordingly. This reordering fixes
% interleaved acquisition as well. Knowing that this will be the case,
% means that the offset to the first slice should always be the min(procpar.pss)
offset = 10*[-procpar.pro 0 min(procpar.pss)];

% make into a translation matrix
offset = [eye(3) offset';0 0 0 1];

% get the distance to the image origin (ie voxel 0,0,0) in number of voxels
% (i.e. the image dimensions divided by 2 as we assume that the offset specifies
% where the center of the volume is)
originOffset = -(dim-1)/2;originOffset(3) = 0;
originOffset = [eye(3) originOffset';0 0 0 1];

% this swaps the dimensions to the coordinate frame that Nifti is expecting.
swapDim =[0 0 1 0;1 0 0 0;0 1 0 0;0 0 0 1];
swapDim2 =[0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1];

% now create the final shifted rotation matrix
xform = swapDim2*rotmat*swapDim*offset*diag(voxspacing)*originOffset;

% testing rotmat
%rotmatpsi = euler2rotmatrix(procpar.psi,0,0);
%rotmattheta = euler2rotmatrix(0,procpar.theta,0);
%rotmatphi = euler2rotmatrix(0,0,-procpar.phi);
%xform = rotmatpsi*rotmattheta*rotmatphi*sliceOffset*swapDim*offset*diag(voxsize)*originOffset;

% round-off to zero
xform((xform < 1e-10) & (xform > - 1e-10)) = 0;

% verbose display only
if verbose > 0
  % display some info.
  disp(sprintf('(fid2xform) psi=%0.2f phi=%0.2f theta=%0.2f',procpar.psi,procpar.phi,procpar.theta));
  disp(sprintf('(fid2xform) Scan dims=[%i %i %i]',dim(1),dim(2),dim(3)));
  disp(sprintf('(fid2xform) Voxel size: [%0.2f %0.2f %0.2f]',voxsize(1),voxsize(2),voxsize(3)));
  disp(sprintf('(fid2xform) Voxel spacing: [%0.2f %0.2f %0.2f]',voxspacing(1),voxspacing(2),voxspacing(3)));
  disp(sprintf('(fid2xform) First slice offset: [%0.2f %0.2f %0.2f]',offset(1,4),offset(2,4),offset(3,4)));
  disp(sprintf('(fid2xform) pss = %s',num2str(procpar.pss)));
  disp(sprintf('(fid2xform) offset to origin: [%0.2f %0.2f %0.2f]',originOffset(1),originOffset(2),originOffset(3)));
end

% work out what the tr is, this is just for setting in the info output field and
% is used by fid2niftihdr
% get the volume TR (framePeriod) for EPI 
tr = procpar.tr;
% if we run mutliple shots, volume TR = slice TR * shots 
if isfield(procpar,'navechoes')
  tr = tr * procpar.numshots/procpar.accfactor;
end
% if we run slice at not interleaved way, 
%  volume TR = slice TR * shots * slice number
if isfield(procpar,'intlv') && strcmp(procpar.intlv, 'n')
  tr = tr * length(procpar.pss);
end

% check for 2d anatomy, to tell getfid that this has the slices and receivers mixed up
if ~info.acq3d && ~info.isepi
  info.receiversAndSlicesSwapped = 1;
  if verbose>0,disp(sprintf('(fid2xform) Receviers and slices are swapped'));,end
else 
  info.receiversAndSlicesSwapped = 0;
end
  

% pack up some useful information that we have learned about the scan
% into an output structure that can be used by other programs.
info.dim = dim;
info.voxsize = voxsize;
info.voxspacing = voxspacing;
info.offset = offset;
info.originOffset = originOffset;
info.pss = procpar.pss;
info.psi = procpar.psi;
info.phi = procpar.phi;
info.theta = procpar.theta;
info.accFactor = procpar.accfactor;
if isfield(procpar,'cntr')
  info.dim(4) = length(procpar.cntr)-1;
else
  info.dim(4) = 1;
end
info.tr = tr;
% get date and time
[info.dateStr info.timeStr] = getDateAndTimeFromVarianField(procpar.time_run);
info.procpar = procpar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getDateAndTimeFromVarianField    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dateStr timeStr] = getDateAndTimeFromVarianField(f)

year = str2num(f{1}(1:4));
month = str2num(f{1}(5:6));
day = str2num(f{1}(7:8));

hour = str2num(f{1}(10:11));
min = str2num(f{1}(12:13));
sec = str2num(f{1}(14:15));

dateStr = sprintf('%04i/%02i/%02i',year,month,day);
timeStr = sprintf('%02i:%02i:%02i',hour,min,sec);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%   euler2rotmatrix   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function r = euler2rotmatrix(a,b,g)

% convert to radians
a = pi*a/180;b = pi*b/180;g = pi*g/180;

% get cos and sin of angles
ca = cos(a);sa = sin(a);
cb = cos(b);sb = sin(b);
cg = cos(g);sg = sin(g);

% convert each rotation into a rotation matrix
arot1 = [ca  sa 0;...
	-sa ca 0;...
	0   0  1];
arot2 = [ca  sa 0;...
	0 1  0;...
	-sa 0 ca];
arot3 = [1  0 0;...
	 0 ca  sa;...
	 0  -sa ca];
brot1 = [cb sb 0;...
	-sb cb  0;...
	0 0 1];
brot2 = [cb 0 sb;...
	0 1  0;...
	-sb 0 cb];
brot3 = [1 0   0;...
	 0 cb  sb;...
	 0  -sb cb];
grot1 = [cg  sg 0;...
	-sg cg 0;...
	0   0  1];
grot2 = [cg  0  sg;...
	 0   1   0;...
         -sg  0 cg];
grot3 = [1 0   0;...
	 0 cg  sg;...
	 0  -sg cg];

% composite the rotations together
r = arot3*brot1*grot3;

% make into a homogenized xform
r = [[r;0 0 0] [0 0 0 1]'];

