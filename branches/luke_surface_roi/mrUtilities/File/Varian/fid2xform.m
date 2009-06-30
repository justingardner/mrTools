% fid2xform.m
%
%        $Id:$ 
%      usage: rotmat = fid2xform(procpar)
%         by: justin gardner
%       date: 05/05/09
%    purpose: Convert fields out of a named procpar into a rotation matrix
%             suitable for the nifti header. (procpar is either the name of 
%             a file, or the structure returned by readprocpar
%
function xform = fid2xform(procpar,verbose)

xform = [];

% check arguments
if ~any(nargin == [1 2])
  help fid2xform
  return
end

if nargin == 1,verbose = 0;end

% get the procpar
if isstr(procpar)
  procpar = readprocpar(procpar);
  if isempty(procpar),return,end;
elseif ~isstruct(procpar)
  help fid2xform;
  return
end

% get the number of navechoes
if isfield(procpar,'navechoes')
  navechoes = procpar.navechoes;
else
  navechoes = 0;
end

% get the dimensions of the scan
dim = [procpar.ni procpar.nv-navechoes length(procpar.pss)];

% make the rotation matrix from the procpar angles
rotmat = euler2rotmatrix(procpar.psi,-procpar.theta,-procpar.phi);

% get voxel sizes and put in diagonal to multiply against rotmat
% so that we get the proper spacing
voxsize = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) procpar.thk 1];

% check for 3d acquisition
if procpar.nv2 > 1
  disp(sprintf('(fid2xform) 3D acquisition'));
  % since the 3rd dimension is taken as a single slice with multiple
  % phase encodes, we have to get the voxel size and dimensions differently
  voxsize(3) = 10*procpar.lpe2/procpar.nv2;
  dim(3) = procpar.nv2;
  dim(1) = dim(2);
  procpar.pss = fliplr(procpar.pss);
  % flip to confrom to 2d images
  rotmat = rotmat*[-1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
end

% Now get the offset in mm from the center of the bore that the center of the
% volume is. We can not change the phase encode center. Note that dimensions
% are given in cm, so we must convert to mm
offset = 10*[-procpar.pro 0 procpar.pss(1)];

% make into a translation matrix
offset = [eye(3) offset';0 0 0 1];

% get the distance to the image origin (ie voxel 0,0,0) in number of voxels
% (i.e. the image dimensions divided by 2 as we assume that the offset specifies
% where the center of the volume is)
originOffset = -(dim-1)/2;originOffset(3) = 0;
originOffset = [eye(3) originOffset';0 0 0 1];

% this swaps the dimensions to the coordinate framte that Nifti is expecting.
swapDim =[0 0 1 0;1 0 0 0;0 1 0 0;0 0 0 1];

% now create the final shifted rotation matrix
xform = rotmat*swapDim*offset*diag(voxsize)*originOffset;

% testing rotmat
%rotmatpsi = euler2rotmatrix(procpar.psi,0,0);
%rotmattheta = euler2rotmatrix(0,procpar.theta,0);
%rotmatphi = euler2rotmatrix(0,0,-procpar.phi);
%xform = rotmatpsi*rotmattheta*rotmatphi*sliceOffset*swapDim*offset*diag(voxsize)*originOffset;

% round-off to zero
xform((xform < 1e-10) & (xform > - 1e-10)) = 0;

% verbose display only
if verbose
  % display some info.
  disp(sprintf('(fid2xform) psi=%0.2f phi=%0.2f theta=%0.2f',procpar.psi,procpar.phi,procpar.theta));
  disp(sprintf('(fid2xform) Scan dims=[%i %i %i]',dim(1),dim(2),dim(3)));
  disp(sprintf('(fid2xform) Voxel size: [%0.2f %0.2f %0.2f]',voxsize(1),voxsize(2),voxsize(3)));
  disp(sprintf('(fid2xform) First slice offset: [%0.2f %0.2f %0.2f]',offset(1),offset(2),offset(3)));
  disp(sprintf('(fid2xform) pss = %s',num2str(procpar.pss)));
  disp(sprintf('(fid2xform) offset to origin: [%0.2f %0.2f %0.2f]',originOffset(1),originOffset(2),originOffset(3)));
end

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

