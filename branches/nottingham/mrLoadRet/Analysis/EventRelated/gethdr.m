% gethdr.m
%
%        $Id$	
%      usage: gethdr(d,x,y,s)
%         by: justin gardner, modified by julien besle
%       date: 08/21/03
%    purpose: retrieve the hdrs in a matrix
%             form from the data structure
%             
function [hdr,time,hdrste, contrastHdr, contrastHdrSte] = gethdr(d,x,y,s,contrasts)

if ~any(nargin == [2 4 5])
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

hdr = shiftdim(d.ehdr(x,y,s,:,:), 3);
hdrlen = d.hdrlen;

% check if this is a glm analysis
% then we scale the hdr 
if isfield(d, 'hrf')
%    hdr = hdr*(d.hrf'-mean(d.hrf));
    hdr = hdr*d.hrf';
    hdrlen = size(d.hrf,1);
end

% return the time as well if asked for
if nargout >= 2
  if ~isfield(d,'estimationSupersampling')
    d.estimationSupersampling =1;
  end
  if ~isfield(d,'acquisitionSubsample')
    d.acquisitionSubsample=1;
  end
  tr = d.tr/d.estimationSupersampling;
  time = (d.acquisitionSubsample-.5)*tr:tr:(hdrlen+d.acquisitionSubsample-1.5)*tr;
end

% return the standard error
if nargout >=3
  hdrste = shiftdim(d.ehdrste(x,y,s,:,:), 3);
  if isfield(d, 'hrf')
     %hdrste = sqrt(hdrste.^2*abs(d.hrf)'); %JB: 
     hdrste = sqrt(hdrste.^2*(d.hrf.^2)'); %the variance (squared stddev) of a random variable multiplied by a constant 
     %is equal to the variance multiplied by the square of the constant. 
  end

end

if ieNotDefined('contrasts')
  contrasts = [];
end
if nargout ==4
  contrastHdr = getContrastEstimate(d,x,y,s,contrasts);
elseif nargout >=5
  [contrastHdr,contrastHdrSte] = getContrastEstimate(d,x,y,s,contrasts);
end

if ~ieNotDefined('contrastHdr') && isfield(d, 'hrf')
  contrastHdr = contrastHdr*d.hrf';
  if ~ieNotDefined('contrastHdrSte')
    contrastHdrSte = sqrt(contrastHdrSte.^2*(d.hrf.^2)'); %the variance (squared stddev) of a random variable multiplied by a constant 
     %is equal to the variance multiplied by the square of the constant. 
  end
end




