% loadScan.m
%
%      usage: loadScan(view,scanNum,groupNum,<sliceNum>)
%         by: justin gardner
%       date: 03/20/07
%    purpose: loads a scan into a "d" structure
%             sliceNum is either [] for all slices,
%             slice numbers for slices you want or 0 to not load data
%
function d = loadScan(view,scanNum,groupNum,sliceNum)

if ~any(nargin == [1 2 3 4])
  help('loadScan');
  return
end

if ~isview(view)
  disp(sprintf('(loadScan) First argument is not a view'));
  return
end

% default to loading all slices
if ~exist('groupNum','var'),groupNum = [];end
if ~exist('sliceNum','var'),sliceNum = [];end

% set the group number
if ~isempty(groupNum)
  view = viewSet(view,'curGroup',groupNum);
end
groupNum = viewGet(view,'curGroup');

% load parameters
d.ver = 4.5;
d.scanNum = scanNum;
d.groupNum = groupNum;
d.description = viewGet(view,'description',scanNum);
d.tr = viewGet(view,'framePeriod',scanNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum);
d.xform = viewGet(view,'scanXform',scanNum);
d.nFrames = viewGet(view,'nFrames',scanNum);
d.dim = viewGet(view,'scanDims',scanNum);
d.dim(4) = d.nFrames;
d.filename = viewGet(view,'tseriesfile',scanNum);
d.filepath = viewGet(view,'tseriespathstr',scanNum);
d.expname = getLastDir(fileparts(fileparts(fileparts(d.filepath))));
d.fullpath = fileparts(fileparts(fileparts(fileparts(d.filepath))));

% print out tr for 3d scans to make sure it is right
if viewGet(view,'3D',scanNum)
  disp(sprintf('(loadScan) 3D sequence. TR is %0.2f',d.tr));
end

% dispay string to say what we are loading
mrDisp(sprintf('Loading scan %i from group: %s',scanNum,viewGet(view,'groupName')));
if isempty(sliceNum)
  mrDisp(sprintf('\n'));
elseif isequal(sliceNum,0)
  mrDisp(sprintf('\n'));
elseif length(sliceNum) == 2
  mrDisp(sprintf(' slices=%i:%i of %i\n',sliceNum(1),sliceNum(2),viewGet(view,'nSlices',scanNum)));
else
  mrDisp(sprintf(' slices=%i of %i\n',sliceNum(1),viewGet(view,'nSlices',scanNum)));
end

% load the data
if ~isequal(sliceNum,0)
  d.data = loadTSeries(view,scanNum,sliceNum);
else
  d.data = [];
end
	
% set the number of slices approriately
d.dim(3) = size(d.data,3);

% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum);
if ~isempty(d.data)
  d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);
end
% junk frames total is used by getStimvol to adjust
% volumes according to how many volumes have been thrown out
d.junkFrames = viewGet(view,'totalJunkedFrames',scanNum);
% we need to add the junk frames junked here to the first
% of the array 'totalJunkedFrames' (one for each scan), or
% if that array is empty then we only have to condiser the
% first junkFrames
if isempty(d.junkFrames)
  d.junkFrames = junkFrames;
else
  d.junkFrames(1) = d.junkFrames(1)+junkFrames;
end
% load dicom header
d.dicom = viewGet(view,'dicom',scanNum);

% load stimfile and set traces
d.stimfile = viewGet(view,'stimfile',scanNum);

if length(d.junkFrames) ~= length(d.stimfile)
  if (d.junkFrames == 0)
    disp(sprintf('(loadScan) Setting total junk frames for all scans to 0'));
    d.junkFrames = zeros(1,length(d.stimfile));
  end
end

% get any concat info
d.concatInfo = viewGet(view,'concatInfo',scanNum);

