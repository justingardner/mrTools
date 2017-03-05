% maskBwithA(A,B,test,<newValue>)
%
%   masks several overlays according to a test on the values of the last
%   selected overlay
%     test is an anonymous function that will be applied to the last overlay, for example @(x)x<=.05 
%     newValue is the value to replace values in the other overlays that DO NOT pass the test (default is nan)
%     
% jb 26/10/2016
%

function varargout = maskWithLast(varargin)
%(A,B,test,newValue)

if strcmp(class(varargin{end}),'function_handle')
  nOverlays = nargin-2;
  test = varargin{end};
  newValue = NaN;
elseif strcmp(class(varargin{end-1}),'function_handle')
  nOverlays = nargin-3;
  test = varargin{end-1};
  newValue = varargin{end};
else
  mrWarnDlg('(maskWithLast)');
  return
end

for o = 1:nOverlays
	varargout{o} = varargin{o};
  varargout{o}(~test(varargin{nOverlays+1})) = newValue;
end
