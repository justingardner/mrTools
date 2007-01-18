% mydetrend.m
%
%      usage: mydetrend(data,dispfig)
%         by: justin gardner
%       date: 07/05/05
%    purpose: detrend a vector. if data is
%             a matrix then the detrending goes
%             along the columns.
%       e.g.: mydetrend(rand(1,51)+(0:1:50)*0.1+3,1)
%
function retval = mydetrend(data,dispfig)

% check command line arguments
if (nargin == 1)
  dispfig = 0;
elseif (nargin ~= 2)
  help mydetrend;
  return
end

% make into column vector
if (size(data,1)) == 1
  data = data';
end

% get length
n = size(data,1);

% make regression matrix
A = [(0:1/(n-1):1)' ones(n,1)];

% get the slope and offsets
regcoef = ((A'*A)^-1)*A'*data;
%regcoef = pinv(A)*data;

% take out the slope
retval = data-A(:,1)*regcoef(1,:);

% plot it
if dispfig
  smartfig('mydetrend');
  for i = 1:size(data,2)
    plot(data,'g');
    hold on
    plot(A*regcoef,'k-');
    plot(retval,'r');
    legend('raw data','regression line','detrended data');
  end
end
