function result = convXYsep(input,xfilter,yfilter,shape) 
%
% result = convXYsep(input,xfilter,yfilter,[shape]) 
%
% input is a volume, 3d array
% xfilter, yfilter are vectors containing filter taps
%
% result is 3d volume in which each slice of input has been
% separably convolved by the two filters.
% 
% shape: 'full', 'same', or 'valid' (default 'same'). See help conv2.  

if ieNotDefined('shape')
    shape = 'same';
end

% convolve one slice to get the dimensions right
tmp=conv2sep(input(:,:,1),xfilter,yfilter,shape);
newdims=size(tmp);
result=zeros(newdims(1),newdims(2),size(input,3));
result(:,:,1)=tmp;

% convolve the other slices
for sliceNum=2:size(input,3)
  tmp=conv2sep(input(:,:,sliceNum),xfilter,yfilter,shape);
  result(:,:,sliceNum)=tmp;
end

return;

% Debug:
input = ones(5,5,5);
filter = [1/4,1/2,1/4];
result = convXYsep(input,filter,filter,'same');
