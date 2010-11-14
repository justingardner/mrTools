% getGlmContrast.m
%
%        $Id$
%      usage: getGlmContrast(d, contrasts)
%         by: farshad moradi
%       date: 09/10/07
%    purpose: function to get r2 image for given combinations of
%             explanatory variables by removing the contribution 
%             of other explanatory variables and then finding the
%             correlation between the residual data and the ortho-
%             goalized combinations (contrasts) 
%       e.g.: d = getr2(d, [1 0 0 0 0])
%             d needs to have a stimulus convolution matrix
%             in d.scm, and contrast must have the same number
%             of columns as d.scm
%
% modified by julien besle xx/12/2009 to speed things up by pre allocating 
%             the whole ehdr and ehdrste matrices with the right dimensions
%             and to compute the residual and total sum of squares of the model
%
% THIS FUNCTION IS DEPRECATED AND HAS BEEN REPLACED BY GETGLMSTATISICS.M

function d = getGlmContrast(d, contrasts)

% number of explanatory variables in the model
nEv = size(d.scm, 2);


cmat = zeros(nEv, nEv);
cmat(1:size(contrasts,1),:) = contrasts;

% design matrix excluding the contrast of (non?) interest
d.scm_ev = d.scm*contrasts';

%---------------------------- find an orthogonal transformation that gives the contrast of interest
nullmat = null(cmat);
if isempty(nullmat)
    scm_no_interest = [];
else
    scm_no_interest = d.scm*nullmat;
    %orthogonalize d.scm_ev w.r.t. scm_no_interest (JB: shouldn't we be doing the reverse ?
    %indeed, if d.scm_ev represents the contrasts of interest, then we don't want to modify it
    %on the other hand, we don't want scm_no_interest to take interesting effects away from us)...
    %d.scm_ev = d.scm_ev - scm_no_interest*(pinv(scm_no_interest)*d.scm_ev);
    scm_no_interest = scm_no_interest - d.scm_ev*(pinv(d.scm_ev)*scm_no_interest);  %my version (JB)
end

% hemodynamic responses are estimated only for the specified contrasts
d.nhdr = size(d.scm_ev,2)/d.hdrlen;


%----------------------------- precalculate the normal equation (this dramatically speeds up things)

if ~isempty(scm_no_interest)
    precalcmatrix1 = ((scm_no_interest'*scm_no_interest)^-1)*scm_no_interest';
    % if this don't work then do pinv
    if sum(isnan(precalcmatrix1(:))) == length(precalcmatrix1(:))
        disp(sprintf('(getGlmContrast) Using pseudo inverse to invert convolution matrix'));
        precalcmatrix1 = pinv(scm_no_interest);
    end
end

precalcmatrix2 = ((d.scm_ev'*d.scm_ev)^-1)*d.scm_ev';
% if this don't work then do pinv
if sum(isnan(precalcmatrix2(:))) == length(precalcmatrix2(:))
  disp(sprintf('(getGlmContrast) Using pseudo inverse to invert convolution matrix'));
  precalcmatrix2 = pinv(d.scm_ev);
end

% check roi
slices = 1:d.dim(3);slicen = length(slices);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);

% preallocate memory
d.ehdr = zeros(d.nhdr*d.hdrlen,d.dim(1),d.dim(2),d.dim(3));    %JB replaces: d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.ehdrste = zeros(d.nhdr*d.hdrlen,d.dim(1),d.dim(2),d.dim(3)); %JB replaces: d.ehdrste = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.r2 = zeros(d.dim(1),d.dim(2),d.dim(3));
d.rss = zeros(d.dim(1),d.dim(2),d.dim(3)); %JB: this is to store the sum of square of the residual error term 
d.tss = zeros(d.dim(1),d.dim(2),d.dim(3)); %JB: this is to store the total sum of square

% turn off warnings to avoid divide by zero warning
warning('off','MATLAB:divideByZero');

% display string
disppercent(-inf,'(getGlmContrast) Calculating r2');
% cycle through images calculating the estimated hdr and r^2s of the 
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
    if isempty(scm_no_interest)
        resid1 = timeseries;
    else
        resid1 = timeseries - scm_no_interest*precalcmatrix1*timeseries;
    end
    d.ehdr(:,:,j,k) = precalcmatrix2*resid1;    %JB replaces: ehdr{j,k} = precalcmatrix2*resid1;
    % calculate error bars, first get sum-of-squares of residual
    % (in percent signal change)
    d.rss(:,j,k) = sum((resid1-d.scm_ev*d.ehdr(:,:,j,k)).^2);    %JB: replaces: sumOfSquaresResidual = sum((resid1-d.scm_ev*ehdr{j,k}).^2););
    d.tss(:,j,k) = sum(resid1.^2);    %JB replaces: 
    % now calculate the sum-of-squares of that error
    % and divide by the degrees of freedom (n-k where n
    % is the number of timepoints in the scan and k is 
    % the number of timepoints in all the estimated hdr)
    rms = d.rss(:,j,k)/(length(d.volumes)-size(d.scm_ev,2));   %JB: replaces: S2 = sumOfSquaresResidual/(length(d.volumes)-size(d.scm_ev,2));
    % now distribute that error to each one of the points
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    d.ehdrste(:,:,j,k) = sqrt(diag(pinv(d.scm_ev'*d.scm_ev))*rms');            %JB replaces: ehdrste{j,k} = sqrt(diag(pinv(d.scm_ev'*d.scm_ev))*S2);
    % calculate variance accounted for by the estimated hdr
    d.r2(:,j,k) = (1-d.rss(:,j,k)'./sum(timeseries.^2))';     %JB: replaces: r2{j,k} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  end
end
disppercent(inf);

d.ehdr = reshape(shiftdim(d.ehdr,1),d.dim(1),d.dim(2),d.dim(3), d.nhdr, d.hdrlen);        %JB
d.ehdrste = reshape(shiftdim(d.ehdrste,1),d.dim(1),d.dim(2),d.dim(3), d.nhdr, d.hdrlen);  %JB

% JB: remove all the following since the matrix has been created directly in the right shape

% reshape matrix. this also seems the fastest way to do things. we
% could have made a matrix in the above code and then reshaped here
% but the reallocs needed to continually add space to the matrix
% seems to be slower than the loops needed here to reconstruct
% the matrix from the {} arrays.
% disppercent(-inf,'(getGlmContrast) Reshaping matrices');
% for i = xvals
%   disppercent((i-min(xvals))/xvaln);
%   for j = yvals
%     for k = slices
%       % get the ehdr
%       d.ehdr(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdr{j,k}(:,i)),[d.hdrlen d.nhdr])';
%       % and the stderror of that
%       d.ehdrste(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdrste{j,k}(:,i)),[d.hdrlen d.nhdr])';
%       % now reshape r2 into a matrix
%       d.r2(i,j,k) = r2{j,k}(i);
%     end
%   end
% end

% display time took
disppercent(inf);

warning('on','MATLAB:divideByZero');
