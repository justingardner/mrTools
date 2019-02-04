% mlrMake.m
%
%        $Id:$ 
%      usage: mlrMake()
%         by: justin gardner
%       date: 05/24/18
%    purpose: 
%
function retval = mlrMake()

% check arguments
if ~any(nargin == [0])
  help mlrMake
  return
end

% get the mrTools directory
mlrTop = fileparts(fileparts(which('mrLoadRet')));

% list of functions that need compiling
compiledFunctionList = {...
    '/mrUtilities/MatlabUtilities/mrDisp.c',...
    '/mrLoadRet/ROI/dijkstrap.cpp',...
	};	    

for iFile = 1:length(compiledFunctionList)
  % get filename
  filename = fullfile(mlrTop,compiledFunctionList{iFile});
  % check for file
  if mlrIsFile(filename)
    % display what we are doing
    disp(sprintf('(mlrMake) mex: %s',filename));
    % mex the file
    mex(filename);
  else 
    % display that we can't find file
    disp(sprintf('(mlrMake) Could not find file: %s', filename));
  end
end

