% indexMax.m: returns index of maximum value 
%
%        $Id$
%
%   example: index = indexMax(x,y,z) returns [2 3] if y(1)>x(1), y(1)>z(1) and z(2)>x(2), z(2)>y(2)
%           

function index = indexMax(varargin)


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

[dump,index] = max(array,[],nDims+1);

%put NaNs if all values are NaNs
index(all(isnan(array),nDims+1))=NaN;

