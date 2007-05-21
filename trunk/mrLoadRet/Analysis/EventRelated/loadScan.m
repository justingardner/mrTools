% loadScan.m
%
%      usage: loadScan(view,scanNum,groupNum,sliceNum)
%         by: justin gardner
%       date: 03/20/07
%    purpose: loads a scan into a "d" structure
%
function d = loadScan(view,scanNum,groupNum,sliceNum)

if ~any(nargin == [1 2 3 4])
  help('loadScan');
  return
end

% default to loading all slices
if ~exist('groupNum','var'),groupNum = [];end
if ~exist('sliceNum','var'),sliceNum = [];end

% set the group number
if ~isempty(groupNum)
  view = viewSet(view,'curGroup',groupNum);
end
% load parameters
d.ver = 4.5;
d.tr = viewGet(view,'framePeriod',scanNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum);
d.nFrames = viewGet(view,'nFrames',scanNum);
d.dim(4) = d.nFrames;
d.filename = viewGet(view,'tseriesfile',scanNum);
d.filepath = viewGet(view,'tseriespathstr',scanNum);
d.expname = getLastDir(fileparts(fileparts(fileparts(d.filepath))));
d.fullpath = fileparts(fileparts(fileparts(fileparts(d.filepath))));

% dispay string to say what we are loading
mrDisp(sprintf('Loading scan %i from group: %s',scanNum,viewGet(view,'groupName')));
if isempty(sliceNum)
  mrDisp(sprintf('\n'));
elseif length(sliceNum) == 2
  mrDisp(sprintf(' slices=%i:%i of %i\n',sliceNum(1),sliceNum(2),viewGet(view,'nSlices',scanNum)));
else
  mrDisp(sprintf(' slices=%i of %i\n',sliceNum(1),viewGet(view,'nSlices',scanNum)));
end
% load the data
d.data = loadTSeries(view,scanNum,sliceNum);
	
% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum);
d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);
% junk frames total is used by getStimvol to adjust
% volumes according to how many volumes have been thrown out
d.junkFrames = viewGet(view,'junkFramesTotal',scanNum);

% load dicom header
d.dicom = viewGet(view,'dicom',scanNum);

% load stimfile and set traces
d.stimfile = viewGet(view,'stimfile',scanNum);

% get any concat info
d.concatInfo = viewGet(view,'concatInfo',scanNum);

% get dimensions
d.dim = size(d.data);

