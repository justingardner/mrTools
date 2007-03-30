% getr2.m
%
%      usage: getr2(d,type)
%         by: justin gardner
%       date: 07/30/03
%    purpose: function to get r2 image by
%             finding hdr for each voxel
%             and seeing how much of the
%             variance the estimate accounts for.
%             also, can calculate impulse response function
%             for mseq stimulus sequence and do same.
%       e.g.: d = getr2(d)
%             d needs to have a stimulus convolution matrix
%             in d.scm
%             d = getr2(d,'mseq')
%             this will calculate the correlation of the mseq
%             with the voxel timeseries. the vairable d.scorrm
%             needs to be set to the correlation matrix. also
%             for mseq, need to specify which volumes to use
%             in the d.mseq.volumes variable. this needs to
%             be an array of length same as the mseq. also
%             need to have the stimulus convolution matrix scm set.
%
function d = getr2(d,type,justr2,roi)

if (nargin == 1)
  type = 'leastsq';
  justr2 = 0;
  roi = [];
elseif (nargin == 2)
  justr2 = 0;
  roi = [];
elseif (nargin == 3)
  roi = [];
elseif (nargin ~= 4)
  help getr2;
  return
end

% init some variables
ehdr=[];r2 = [];

if (strcmp(type,'leastsq'))
  % precalculate the normal equation (this dramatically speeds up things)
  precalcmatrix = ((d.scm'*d.scm)^-1)*d.scm';
  % if this don't work then do pinv
  if length(isnan(precalcmatrix(:))) == length(precalcmatrix(:))
    disp(sprintf('Using pseudo inverse to invert convolution matrix'));
    precalcmatrix = pinv(d.scm);
  end
elseif (strcmp(type,'mseq'))
  % for mseq we just need to multiply by the correlation matrix
  precalcmatrix = d.scorrm;
else
  help getr2;
  return
end

% check roi
if (isempty(roi))
  % no roi, so just do whole volume
  slices = 1:d.dim(3);slicen = length(slices);
  xvals = 1:d.dim(1);xvaln = length(xvals);
  yvals = 1:d.dim(2);yvaln = length(yvals);
else
  % make sure this is just a singleton roi since
  % we can only do a box at a time.
  if (roi.n > 1)
    disp(sprintf('roi must have only one box.'));
    return
  end
  roi = getroixys(roi,d.dim(1));
  slices = roi.slice(1);slicen = 1;
  xvals = min(roi.x):max(roi.x);xvaln = length(xvals);
  yvals = min(roi.y):max(roi.y);yvaln = length(yvals);
end  
  
% preallocate memory
if (~justr2)
  d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),size(precalcmatrix,1));
  d.r2 = zeros(d.dim(1),d.dim(2),d.dim(3));
else
  R2 = zeros(d.dim(1),d.dim(2),d.dim(3));
end  
% display string
disppercent(-inf,'Calculating r2');
% cycle through images calculating the estimated hdr and r^s of the 
% estimate.
%
% this following section has been optimized to run faster by
% eliminating one of the loops. Various different methods were
% tested eliminating all the loops and doing one big calculation
% which thrashed memory too much, or eliminating different
% dimensions and it was found that eliminating the first dimension
% was by far the faster by a factor of about 2-3. 
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  disppercent(max((j-min(yvals))/yvaln,0.1));
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % convert to percent signal change
    timeseries = 100*timeseries./(onesmatrix*colmeans);
    % get hdr for the each voxel
    ehdr{j,k} = precalcmatrix*timeseries;
    % calculate error bars, first get sum-of-squares of residual
    % (in percent signal change)
    sumOfSquaresResidual = sum((timeseries-d.scm*ehdr{j,k}).^2);
    % now calculate the sum-of-squares of that error
    % and divide by the degrees of freedom (n-k where n
    % is the number of timepoints in the scan and k is 
    % the number of timepoints in all the estimated hdr)
    S2 = sumOfSquaresResidual/(length(d.volumes)-size(d.scm,2));
    % now distribute that error to each one of the points
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    ehdrste{j,k} = sqrt(diag(pinv(d.scm'*d.scm))*S2);
    % calculate variance accounted for by the estimated hdr
    r2{j,k} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  end
end
disppercent(inf);
if (~justr2)
  % reshape matrix. this also seems the fastest way to do things. we
  % could have made a matrix in the above code and then reshaped here
  % but the reallocs needed to continually add space to the matrix
  % seems to be slower than the loops needed here to reconstruct
  % the matrix from the {} arrays.
  disppercent(-inf,'Reshaping matrices');
  for i = xvals
    disppercent((i-min(xvals))/xvaln);
    for j = yvals
      for k = slices
	% get the ehdr
	d.ehdr(i,j,k,:) = ehdr{j,k}(:,i);
	% and the stderror of that
	d.ehdrste(i,j,k,:) = ehdrste{j,k}(:,i);
	% now reshape into a matrix
	d.r2(i,j,k) = r2{j,k}(i);
      end
    end
  end
else
  % does same thing, but just returns r2 values,does set anything
  % in data structure
  disppercent(-inf,'Reshaping matrices');
  for i = xvals
    disppercent((i-min(xvals))/xvaln);
    for j = yvals
      for k = slices
	% now reshape into a matrix
	R2(i,j,k) = r2{j,k}(i);
      end
    end
  end
  d = R2;
end  


% display time took
disppercent(inf);

