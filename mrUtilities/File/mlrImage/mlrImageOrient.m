% mlrImageOrient.m
%
%      usage: [data h] = mlrImageOrient(orient,data,h)
%         by: justin gardner
%       date: 09/05/11
%    purpose: re-orient image to a standard orientation like LPI
%             This will fix the header to have the correct
%             qform/sform
%
%       e.g.: [d h] = mlrImageOrient('LPI',data,h);
%  
%
function [data h] = mlrImageOrient(orient,data,h)

% check arguments
if ~any(nargin == [3])
  help mlrImageOrient
  return
end

% for now just support LPI - should be fairly trivial
% to handle other orientations by swapping axisReverseMapping
% and axisDir below appropriately - but be careful of the getting
% the rows in the xform correct
if ~strcmp(orient,'LPI')
  disp(sprintf('(mlrImageOrient) Only LPI supported'));
  return
end

% check header
if ~mlrImageIsHeader(h)
  disp(sprintf('(mlrImageOrient) Header is not a standard mlrImage header'));
  return
end

% check for qform
if isempty(h.qform)
  disp(sprintf('(mlrImageOrient) No qform available in image header'));
  return
end

% get the axis directions
[axisLabels axisDirLabels axisMapping axisReverseMapping axisDir] = mlrImageGetAxisLabels(h.qform);

% flip each axis that goes in the opposite direction
for i = 1:3
  if axisDir(i) == -1
    data = flipdim(data,i);
  end
end

% get the permutation order
permutationOrder = 1:h.nDim;
permutationOrder(1:3) = axisReverseMapping;

% and permute
data = permute(data,permutationOrder);

% make the xform matrix
for i = 1:3
  xform(i,1:4) = 0;
  xform(i,axisReverseMapping(i)) = axisDir(i);
  if axisDir(i) == -1
    xform(i,4) = h.dim(i)-1;
  end
end
xform(4,:) = [0 0 0 1];

h.qform = h.qform*xform;
if ~isempty(h.sform)
  h.sform = h.sform*xform;
end
