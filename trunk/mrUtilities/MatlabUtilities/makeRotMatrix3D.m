% makeRotMatrix3D.m
%
%        $Id$ 
%      usage: makeRotMatrix3D(alpha,beta,gamma,<offset>)
%         by: justin gardner
%       date: 08/14/08
%    purpose: create 3x3 rotation matrix, we use the formula for the "Euler angles"
%             alpha is an inital rotation around the y axis.
%             beta is the subsequent rotation around where the x-axis has been rotated to
%             gamma is the subsequent rotation around where the z-axis has been rotated to.
%             angles are in radians
%             see: http://en.wikipedia.org/wiki/Euler_angles
% 
%             Offset should be a 3x1 vector specifying a shift of coordinates.
%             If offset is defined then this function will return a 4x4 homogenous xform
%
function R = makeRotMatrix3D(alpha,beta,gamma,offset)

% check arguments
if ~any(nargin == [1 2 3 4])
  help makeRotMatrix3D
  return
end

% default to zero
if ieNotDefined('alpha'),alpha = 0;end
if ieNotDefined('beta'),beta = 0;end
if ieNotDefined('gamma'),gamma = 0;end

% short cuts for cos/sin of the three angles
ca = cos(alpha);sa = sin(alpha);
cb = cos(beta);sb = sin(beta);
cg = cos(gamma);sg = sin(gamma);

% the rotation matrix
R = [ca*cg-sa*cb*sg -ca*sg-sa*cb*cg sb*sa;
     sa*cg+ca*cb*sg -sa*sg+ca*cb*cg -sb*ca;
     sb*sg          sb*cg           cb];

% set offset
if ~ieNotDefined('offset')
  R(1,4) = offset(1);
  R(2,4) = offset(2);
  R(3,4) = offset(3);
  R(4,:) = [0 0 0 1];
end