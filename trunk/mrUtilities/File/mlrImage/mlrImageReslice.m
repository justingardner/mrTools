% mlrImageReslice.m
%
%      usage: [d h] = mlrImageReslice(fromImage,toImage)
%         by: justin gardner
%       date: 09/07/11
%    purpose: Reslices the volume named by fromFilename into the
%             same dimensions as the toFilename. This is useful
%             for instance if you did not collect a high-resolution
%             inplane anatomy and you need it for some reason. Then
%             for the fromFilename use your hi-res canonical
%             anatomy and your toFilename put your epi. This will
%             make a new anatomy image with twice the resolution as
%             the epis. The images can be filenames or any of the
%             other formats that are understood by mlrImageParseArgs
%   
%
function [resliceData resliceHeader] = mlrImageReslice(varargin)

% default return arguments
resliceData = [];resliceHeader = [];

% check arguments
if nargin < 1
  help mlrImageReslice
  return
end

% parse args
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
verbose = [];
getArgs(otherArgs,{'verbose=0'});

% check the images
if length(imageArgs) ~= 2
  disp(sprintf('(mlrImageReslice) Need to specify two images, from and a to image'));
  return
end

% make sure we can reload image
if ~mlrImageIsImage(imageArgs{1}) || ~mlrImageIsImage(imageArgs{2})
  disp(sprintf('(mlrImageReslice) Could not load image'));
  return
end
  
% load the images
disppercent(-inf,'Loading images');
[fromData fromHeader] = mlrImageLoad(imageArgs{1});
if isempty(fromData),return,end
disppercent(0.5);
[toData toHeader] = mlrImageLoad(imageArgs{2});
if isempty(toData),return,end
disppercent(inf);

% xform according to sforms if both exist
if ~isempty(fromHeader.sform) && ~isempty(toHeader.sform)
  vol2vol = inv(fromHeader.sform)*toHeader.sform;
  if verbose
    disp(sprintf('(mlrImageReslice) Transforming use sforms\n%s',mlrnum2str(vol2vol,'compact=0')));
  end
elseif ~isempty(fromHeader.qform) && ~isempty(toHeader.qform)
  vol2vol = inv(fromHeader.qform)*toHeader.qform;
  if verbose
    disp(sprintf('(mlrImageReslice) Transforming use qforms\n%s',mlrnum2str(vol2vol,'compact=0')));
  end
else
  disp(sprintf('(mlrImagerReslice) Could not find valid xforms to do reslice'));
  return
end

resFactor = [2 2 2];

% make x, y and z of image
x = 1:1/resFactor(1):toHeader.dim(1)+1/resFactor(1);
y = 1:1/resFactor(2):toHeader.dim(2)+1/resFactor(2);
z = 1:1/resFactor(3):toHeader.dim(3)+1/resFactor(3);

% make the grid
[x y z] = meshgrid(x,y,z);

% get the reslice dimensions
resliceDim = size(x);

% make into homgenous coordinates
toCoord(1,:) = x(:);
toCoord(2,:) = y(:);
toCoord(3,:) = z(:);
toCoord(4,:) = 1;

% swap 
swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];

% take the xform
fromCoord = swapXY*inv(shiftOriginXform)*vol2vol*shiftOriginXform*swapXY*toCoord;

% get the resliced image
interpMethod = mrGetPref('interpMethod');
resliceData = interp3(fromData,fromCoord(1,:),fromCoord(2,:),fromCoord(3,:),interpMethod);
resliceData = reshape(resliceData,resliceDim(1),resliceDim(2),resliceDim(3));

% and make the header
resliceHeader = toHeader;
resliceHeader.dim = size(resliceData);
resliceHeader.pixdim(1:3) = toHeader.pixdim(1:3);

% create the xform
scaleXform = eye(4);
scaleXform(1:3,1:3) = diag(1./resFactor);
offsetXform = eye(4);
offsetXform(1,4) = y(1)-1;
offsetXform(2,4) = x(1)-1;
offsetXform(3,4) = z(1)-1;

if ~isempty(resliceHeader.qform)
  resliceHeader.qform = resliceHeader.qform*offsetXform*scaleXform;
end
if ~isempty(resliceHeader.sform)
  resliceHeader.sform = resliceHeader.sform*offsetXform*scaleXform;
end

