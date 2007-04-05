% gethdr.m
%
%      usage: gethdr(d,x,y,s)
%         by: justin gardner
%       date: 08/21/03
%    purpose: retrieve the hdrs in a matrix
%             form from the data structure
%
function [hdr,time,hdrste] = gethdr(d,x,y,s)

if ~any(nargin == [2 4])
  help gethdr;
  return
end

% if only two arguments, it must be a 3 tuple
if nargin == 2
  if length(x) == 3
    y = x(2);
    s = x(3);
    x = x(1);
  else
    help gethdr;
    return
  end
end

hdr = squeeze(d.ehdr(x,y,s,:,:));
hdrste = squeeze(d.ehdrste(x,y,s,:,:));
time = d.tr/2:d.tr:(d.hdrlen*d.tr);


