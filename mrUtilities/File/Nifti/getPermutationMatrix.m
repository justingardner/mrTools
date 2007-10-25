% getPermutationMatrix.m
%
%      usage: permutationMatrix = getPermutationMatrix(hdr)
%         by: david heeger
%       date: 10/24/07
%    purpose: extracts permutation matrix form nifti hdr
%
function permutationMatrix = getPermutationMatrix(hdr)

% check arguments
if ~any(nargin == [1])
  help getPermutationMatrix
  return
end

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(hdr.qform44(1:3,1:3)));
permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);
