% gethdr.m
%
%        $Id: getContrastEstimate.m 1891 2010-11-25 08:06:27Z julien $	
%      usage: getContrastEstimate(d,x,y,s)
%         by: julien besle, modified from gethdr by Justin Gardner
%       date: 29/11/2010
%    purpose: retrieve the contrast value and standard error from d structure. 
%             If no contrast standard error is present, computes the OLS standard error
%             
function [contrastBetas, contrastBetasSte] = getContrastEstimate(d,x,y,s,contrasts)

if ~any(nargin == [2 4 5])
  help getContrastEstimate;
  return
end

% if only two arguments, it must be a 3 tuple
if nargin == 2
  if length(x) == 3
    y = x(2);
    s = x(3);
    x = x(1);
  else
    help getContrastEstimate;
    return
  end
end

if ieNotDefined('contrasts')
  if ~isfield(d,'contrasts')
    mrWarnDlg('(getContrastEstimate) No contrast is defined');
    return;
  else
    contrasts = d.contrasts;
  end
end

betas = shiftdim(d.ehdr(x,y,s,:,:), 3);
if size(contrasts,2)~= size(betas,1)
   mrErrorDlg('(getContrastEstimate) size(contrasts,2) should be equal to the number of EVs');
end

contrastBetas = contrasts*betas;
if nargout>1
  if isfield(d,'contrastSte') && size(d.contrastSte,4)==size(contrastBetas,1) && size(d.contrastSte,5)==size(contrastBetas,2)
    contrastBetasSte = shiftdim(d.contrastSte(x,y,s,:,:), 3);
  else
    %this only works for Ordinary Least Squares
    if ~isfield(d,'rdf')
       d.rdf = size(d.scm,1)-size(contrasts,2);
    end
    if ~isfield(d,'scm')
      contrastBetasSte = [];
      mrWarnDlg('(getContrastEstimate) No design matrix in the d structure');
      return
    end
    contrasts = kron(contrasts,eye(size(betas,2)));
    contrastBetasSte = sqrt(d.rss(x,y,s)/d.rdf*diag(contrasts/(d.scm'*d.scm)*contrasts'));
  end
  contrastBetasSte = reshape(contrastBetasSte,size(contrastBetas,2),size(contrastBetas,1))';
end

