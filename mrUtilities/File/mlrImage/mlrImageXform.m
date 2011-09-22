% mlrImageXform.m
%
%      usage: [d h] = mlrImageXform(image,commands)
%         by: justin gardner
%       date: 09/10/11
%    purpose: transforms an image while also adjusting header info including
%             xformation matrices accordingly. command is any of
%             the following:
%
%             swapXY,swapXZ,swapYZ
%             flipX,flipY,flipZ
%
%             shiftX=n,shiftY=n,shiftZ=n
%             rotateXY=deg,rotateYZ=deg,rotateXZ=deg
%
%       e.g.: [d h] = mlrImageXform('image.hdr','swapXY=1','flipZ=1','shiftX=13.2',rotateXY=30');
%
%
function [d h] = mlrImageXform(varargin)

% check arguments
if nargin < 1
  help mlrImageXform
  return
end

% parse arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
verbose=[];
swapXY=[];swapXZ=[];swapYZ=[];
flipX=[];flipY=[];flipZ=[];
shiftX=[];shiftY=[];shiftZ=[];
rotateXY=[];rotateXZ=[];rotateYZ=[];
interpMethod=[];
applyToHeader=[];applyToData=[];
getArgs(otherArgs,{'verbose=1','swapXY=0','swapXZ=0','swapYZ=0','flipX=0','flipY=0','flipZ=0','shiftX=0','shiftY=0','shiftZ=0','rotateXY=0','rotateXZ=0','rotateYZ=0','interpMethod=[]','applyToHeader=1','applyToData=1'});

% check that we have an image to xform
if length(imageArgs) < 1
  disp(sprintf('(mlrImageXform) Must specify atleast one image to xform'));
  return
end

% get default interpMethod for xforms that use interp3
if isempty(interpMethod)
  interpMethod = mrGetPref('interpMethod');
end

% now cycle through images xforming them and putting in output
for iImage = 1:length(imageArgs)
  % load the image
  [d h] = mlrImageLoad(imageArgs{iImage},'returnRaw=1');
  if isempty(d)
    disp(sprintf('(mlrImageXform) Could not open image: %s',mlrImageArgFilename(imageArgs{iImage})));
  else
    % apply swapXY
    if swapXY
      if verbose,disp(sprintf('(mlrImageXform) Swapping XY'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(1:2) = [2 1];
	d = permute(d,permuteDims);
      end
      % swap the header
      h = applyXform([0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1],h,d,applyToHeader);
    end
    % apply swapXZ
    if swapXZ
      if verbose,disp(sprintf('(mlrImageXform) Swapping XZ'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(1:3) = [3 2 1];
	d = permute(d,permuteDims);
      end
      % swap the header
      h = applyXform([0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1],h,d,applyToHeader);
    end
    % apply swapYZ
    if swapYZ
      if verbose,disp(sprintf('(mlrImageXform) Swapping YZ'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(2:3) = [3 2];
	d = permute(d,permuteDims);
      end
      % apply to header
      h = applyXform([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1],h,d,applyToHeader);
    end
    % apply flipX
    if flipX
      if verbose,disp(sprintf('(mlrImageXform) Flipping X'));end
      % flip data
      if applyToData
	d = flipdim(d,1);
      end
      % apply to header
      h = applyXform([-1 0 0 h.dim(1)-1;0 1 0 0;0 0 1 0;0 0 0 1],h,d,applyToHeader);
    end
    % apply flipY
    if flipY
      if verbose,disp(sprintf('(mlrImageXform) Flipping Y'));end
      % flip data
      if applyToData
	d = flipdim(d,2);
      end
      % apply to header
      h = applyXform([1 0 0 0;0 -1 0 h.dim(2)-1;0 0 1 0;0 0 0 1],h,d,applyToHeader);
    end
    % apply flipZ
    if flipZ
      if verbose,disp(sprintf('(mlrImageXform) Flipping Z'));end
      % flip data
      if applyToData
	d = flipdim(d,3);
      end
      % apply to header
      h = applyXform([1 0 0 0;0 1 0 0;0 0 -1 h.dim(3)-1;0 0 0 1],h,d,applyToHeader);
    end
    % apply shifts
    if any([shiftX shiftY shiftZ])
      if verbose,disp(sprintf('(mlrImageXform) Using %s to shift by: [%s]',interpMethod,mlrnum2str([shiftX shiftY shiftZ])));end
      % flip data
      if applyToData
	x = 1+shiftX:h.dim(1)+shiftX;
	y = 1+shiftY:h.dim(2)+shiftY;
	z = 1+shiftZ:h.dim(3)+shiftZ;
	[x y z] = ndgrid(x,y,z);
	d = interpn(d,x,y,z,interpMethod);
      end
      % apply to header
      h = applyXform([1 0 0 shiftX;0 1 0 shiftY;0 0 1 shiftZ;0 0 0 1],h,d,applyToHeader);
    end
    % apply rotation
    if any([rotateXY rotateXZ rotateYZ])
      if verbose,disp(sprintf('(mlrImageXform) Using %s to rotate: %s',interpMethod,mlrnum2str([rotateXY rotateXZ rotateYZ])));end
      % flip data
      if applyToData
	% shift coordinates to center
	x = (1:h.dim(1))-h.dim(1)/2;
	y = (1:h.dim(2))-h.dim(2)/2;
	z = (1:h.dim(3))-h.dim(3)/2;
	% make into grid, note the swapping of x and y
	[x y z] = ndgrid(x,y,z);
	% make into homogenous coordinates
	coords(1,:) = x(:);
	coords(2,:) = y(:);
	coords(3,:) = z(:);
	coords(4,:) = 1;
	% now multiply against the rotation matrix
 	rotMatrix = makeRotMatrix3D(pi*rotateXZ/180,pi*rotateYZ/180,pi*rotateXY/180,[0 0 0]);
	coords = rotMatrix*coords;
	% and shift back from center
	coords(1,:) = coords(1,:) + h.dim(1)/2;
	coords(2,:) = coords(2,:) + h.dim(2)/2;
	coords(3,:) = coords(3,:) + h.dim(3)/2;
	% and interp
	d = interpn(d,coords(1,:),coords(2,:),coords(3,:),interpMethod);
	d = reshape(d,h.dim(1),h.dim(2),h.dim(3));
      end
      % apply to header
      shiftToCenter = [1 0 0 h.dim(1)/2;0 1 0 h.dim(2)/2;0 0 1 h.dim(3)/2;0 0 0 1];
      h = applyXform(shiftToCenter*rotMatrix*inv(shiftToCenter),h,d,applyToHeader);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%%   applyXform   %%
%%%%%%%%%%%%%%%%%%%%
function h = applyXform(xform,h,d,applyToHeader)

% nothing to do if we are not applying to header
if ~applyToHeader,return,end

% apply to qform if it is present
if ~isempty(h.qform)
  h.qform = h.qform*xform;
end

% apply to sform if it is present
if ~isempty(h.sform)
  h.sform = h.sform*xform;
end

% get the dimensions of the data
sized = size(d);
h.dim(1:length(sized)) = sized;
