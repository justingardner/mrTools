% weightedMeanStd.m: computes weighted sample average and std deviation of overlay indices (indices are weighted by the overlay value)
%             This is useful to estimate preferred condition and (gaussian) tuning width in a set of ordered stimulus conditions
%             Negative overlay values are set to 0 beofre computing meand and std deviation
%
%        $Id: weightedMeanStd.m 2733 2013-05-13 11:47:54Z julien $
%
%   example: [average,std] = weightedMeanStd(x,y,z) returns average = (1*x+2*y+3*z)/(x+y+z)
%                                                    and stddev = (x*(1-average)^2 + y*(2-average)^2 + z*(3-average)^2) / (x+y+z)
%           

function [average,stddev] = weightedMeanStd(varargin)


nDims = length(size(varargin{1}));
array = varargin{1};
for i = 2:nargin
  %first check that all inputs have the same size
  if ~isequal(size(varargin{i}),size(varargin{1}))
    error('All inputs must have the same size')
  else
    %concatenate
    array=cat(nDims+1,array,varargin{i});
  end
end
overlaySize = size(varargin{1});
nOverlays = size(array,4);
array(array<=0)=0;

%compute the weighted average
indices = repmat(permute(1:nOverlays,[1 3 4 2]),[overlaySize 1]);
average = sum( array .* indices, 4) ./ sum(array,4);
% and weighted stddev
stddev = sqrt( sum( array .* (indices - repmat(average,[1 1 1 nOverlays]) ).^2, 4) ./ sum(array,4) );


