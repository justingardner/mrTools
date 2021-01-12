% weightedMeanStd.m: computes weighted sample average and std deviation of overlay indices (indices are weighted by the overlay value)
%             This is useful to estimate preferred condition and (gaussian) tuning width in a set of ordered stimulus conditions
%             Negative overlay values are set to 0 before computing mean and std deviation
%             Sample average and std deviation are also corrected for limited stimulus sampling range
%             using the method described in Besle et al. (2019) Cerebral Cortex, and returned as 3rd and 4th outputs
%             Optional argument: population mean range [a b], where a is the lower bound in overlay index scale 
%             (usually something between 0  and 1, default = 0) and b is a multiplier applied to the total number of overlays (default = 1.3)
%
%        $Id: weightedMeanStd.m 2733 2013-05-13 11:47:54Z julien $
%
%   example: [average,std] = weightedMeanStd(x,y,z) returns average = (1*x+2*y+3*z)/(x+y+z)
%                                                    and stddev = (x*(1-average)^2 + y*(2-average)^2 + z*(3-average)^2) / (x+y+z)
%           

function [average,stddev,correctedAverage,correctedStddev] = weightedMeanStd(varargin)

nDims = length(size(varargin{1}));
array = varargin{1};
for i = 2:nargin
  %first check that all inputs have the same size
  if ~isequal(size(varargin{i}),size(varargin{1}))
    if length(varargin{i})==2 && i==nargin
      populationMeanRange = varargin{i};
      mrWarnDlg('(weightedMeanStd) Last input has two values. Assuming that this is the population mean range');
      break;
    else
      error('(weightedMeanStd) All inputs must have the same size');
    end
  else
    %concatenate
    array=cat(nDims+1,array,varargin{i});
  end
end
overlaySize = size(varargin{1});
nOverlays = size(array,4);
array(array<=0)=0;

if ieNotDefined('populationMeanRange')
  populationMeanRange = [0 1.3];
end
%compute the weighted average
indices = repmat(permute(1:nOverlays,[1 3 4 2]),[overlaySize 1]);
average = sum( array .* indices, 4) ./ sum(array,4);
% and weighted stddev
stddev = sqrt( sum( array .* (indices - repmat(average,[1 1 1 nOverlays]) ).^2, 4) ./ sum(array,4) );


%correct partial-sampling bias
if nargout == 4
  popMean = populationMeanRange(1):.1:nOverlays* populationMeanRange(2);
  popStddev = .1:.1:nOverlays;
  nMeans = length(popMean);
  nStddevs = length(popStddev);
  [sampleAverage,sampleStddev] = meanEstimationBias(popMean,popStddev,1:nOverlays);
  
  nonNaN = find(all(~isnan(array),4))';
  correctedAverage = NaN(overlaySize);
  correctedStddev = NaN(overlaySize);
  hWaitBar = mrWaitBar(-inf,'(weightedMeanStd) Correcting bias');
  c = 0;
  for i=nonNaN
    c = c+1;
    minimization = abs(sampleAverage-average(i))+abs(sampleStddev-stddev(i));
    [minValue,minIndex] = min(minimization(:));
    [whichMean,whichStddev] = ind2sub([nMeans nStddevs],minIndex);
    correctedAverage(i) = popMean(whichMean);
    correctedStddev(i) = popStddev(whichStddev);
    mrWaitBar( c/numel(nonNaN), hWaitBar);
  end
  mrCloseDlg(hWaitBar);
end
    