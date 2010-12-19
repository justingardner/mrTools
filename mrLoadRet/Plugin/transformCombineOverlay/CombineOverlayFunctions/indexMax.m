% indexMax.m: returns index of maximum value 
%
%        $Id$
%
%   example: index = max(x,y,z) returns 2 if y>x and y>z

function index = indexMax(varargin)

numbers = cell2mat(varargin);
[dump,index] = max(numbers);

