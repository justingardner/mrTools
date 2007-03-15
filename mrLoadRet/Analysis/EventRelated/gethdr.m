% gethdr.m
%
%      usage: gethdr(d,x,y,s)
%         by: justin gardner
%       date: 08/21/03
%    purpose: retrieve the hdrs in a matrix
%             form from the data structure
%
function [hdr,time] = gethdr(d,x,y,s)

if (nargin ~= 4)
  help gethdr;
  return
end

for i = 1:d.nhdr
  hdr(i,1:d.hdrlen) = d.ehdr(x,y,s,(i-1)*d.hdrlen+1:(i)*d.hdrlen);
end

time = d.tr/2:d.tr:(d.hdrlen*d.tr);

