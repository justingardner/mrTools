%        $Id$

function index = indexMax(varargin)

numbers = cell2mat(varargin);
[dump,index] = max(numbers);

