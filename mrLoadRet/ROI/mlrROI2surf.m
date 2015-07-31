% mlrROI2surf.m
%
%        $Id:$ 
%      usage: surf = mlrROI2surf(roi)
%         by: justin gardner
%       date: 07/30/15
%    purpose: function that will convert an mlr ROI into a surf (triangulated mesh with
%             vertices and triangles)
%
%             v = newView;
%             v = loadROI(v,'l_mt.mat');
%             roi = viewGet(v,'roi','l_mt');
%             surf = mlrROI2surf(roi);
%
%
function retval = mlrROI2surf(roi)

% check arguments
if ~any(nargin == [1])
  help mlrROI2surf
  return
end


keyboard