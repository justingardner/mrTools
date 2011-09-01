% fid2niftihdr.m
%
%        $Id:$ 
%      usage: hdr = fid2niftihdr(fidname,<verbose>,<movepro=0>)
%         by: justin gardner
%       date: 12/09/09
%    purpose: strips out of fid2nifti the functionality to create nifti headers from procpar. Fidname
%             is either the name of a fid directory *or* can be a procpar structure returned by readprocpar.
%             
%
function [hdr info] = fid2niftihdr(fidname,verbose,varargin)

hdr = [];info = [];
if ieNotDefined('verbose'),verbose=1;end

movepro=[];
getArgs(varargin,{'movepro=0'});

% check arguments
if ~any(nargin == [1 2 3])
  help fid2niftihdr
  return
end

% set fid extension
if ~isstruct(fidname)
  fidname = setext(fidname,'fid',0);
end

% create an empty header
hdr = cbiCreateNiftiHeader;

% get the qform
[qform44 info] = fid2xform(fidname,verbose,sprintf('movepro=%f',movepro));
if isempty(qform44),hdr = [];return;end
  
% set the qform
hdr = cbiSetNiftiQform(hdr,qform44);

% now set dimensions in header
nDims = length(info.dim);
hdr.dim(1:nDims+1) = [nDims info.dim];

% set voxel dimensions 
hdr.pixdim(2:4) = info.voxsize(1:3);

% set the tr
hdr.pixdim(5) = info.tr*1000;



