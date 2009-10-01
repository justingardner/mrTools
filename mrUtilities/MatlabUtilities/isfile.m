% isfile.m
%
%      usage: isfile(filename)
%         by: justin gardner
%       date: 08/20/03
%       e.g.: isfile('filename')
%    purpose: function to check whether file exists
%
function [isit permission] = isfile(filename)

isit = 0;permission = [];
if (nargin ~= 1)
  help isfile;
  return
end
if isempty(filename),isit = 0;,return,end

% open file
fid = fopen(filename,'r');

% check to see if there was an error
if (fid ~= -1)
  fclose(fid);
  [dummy permission] = fileattrib(filename);
  isit = 1;
else
  isit = 0;
end

