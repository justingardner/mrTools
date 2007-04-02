% reshapeTSeries
%
%      usage: tseries = reshapeTSeries(tseries)
%         by: eli merriam
%       date: 03/22/07
%    purpose: 
%	$Id$	
%
function tseries = reshapeTSeries(tseries)

% % check arguments
% if ~any(nargin == [0])
%   help reshapeTSeries
%   return
% end

view = newView('Volume');
mrGlobals

% Load single slice and reshape to nFrames X nVoxels
tseries = squeeze(tseries);
% Reformat to nFrames x nVoxels
dims = size(tseries);
nFrames = dims(3);
sliceDims = dims(1:2);
nVoxels = prod(sliceDims);
tseries = reshape(tseries,[nVoxels,nFrames])';

% % Check # voxels and # temporal frames
% if (nFrames ~= viewGet(view,'totalFrames',scan))
%     mrWarnDlg('loadTSeries: number of frames in tseries file does not match expected.');
% end
% if sliceDims ~= viewGet(view,'sliceDims',scan)
%     mrWarnDlg('loadTSeries: number of voxels in tseries does not match expected.');
% end
