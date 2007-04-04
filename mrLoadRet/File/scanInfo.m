% scanInfo.m
%
%      usage: scanInfo(scanNum,groupNum)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: groupInfo(1);
%             groupInfo('Raw');
function retval = scanInfo(scanNum,groupNum)

% check arguments
if ~any(nargin == [2])
  help scanInfo
  return
end

view = newView('Volume');

% if groupNum is a string, then user passed in a name rather than
% a number
if isstr(groupNum)
  groupName = groupNum;
  groupNum = viewGet(view,'groupNum',groupName);
  if isempty(groupNum)
    disp(sprintf('(groupInfo): No group %s',groupName));
    return
  end
end

if (groupNum < 1) | (groupNum > viewGet(view,'numberOfGroups'))
  disp(sprintf('(groupInfo): No group number %i',groupNum));
  return
end
groupName = viewGet(view,'groupName',groupNum);
% set the group and scan
view = viewSet(view,'curGroup',groupNum);

% grab info
description = viewGet(view,'description',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
tr = viewGet(view,'framePeriod',scanNum,groupNum);
totalFrames = viewGet(view,'totalFrames',scanNum,groupNum);
filename = viewGet(view,'tSeriesFile',scanNum,groupNum);
originalFilename = viewGet(view,'originalFilename',scanNum,groupNum);
originalGroupname = viewGet(view,'originalGroupname',scanNum,groupNum);
stimFilename = viewGet(view,'stimFilename',scanNum,groupNum);
scanHdr = viewGet(view,'niftiHdr',scanNum,groupNum);
scanDims = viewGet(view,'scanDims',scanNum,groupNum);

% display info
disp(sprintf('%s',description));
disp(sprintf('Filename: %s GroupName: %s',filename,groupName));
for i = 1:length(originalFilename)
  disp(sprintf('Original Filename: %s Group: %s',originalFilename{i},originalGroupname{i}));
end
for i = 1:length(stimFilename)
  disp(sprintf('StimFilename: %s',stimFilename{i}));
end
  
disp(sprintf('voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f Dims: [%i %i %i] Volumes=%i',scanVoxelSize(1),scanVoxelSize(2),scanVoxelSize(3),tr,scanDims(1),scanDims(2),scanDims(3),totalFrames));

% display qform and sform
disp(sprintf('++++++++++++++++++++++++++ qform ++++++++++++++++++++++++++'));
for rownum = 1:4
  disp(sprintf('%f\t%f\t%f\t%f',scanHdr.qform44(rownum,1),scanHdr.qform44(rownum,2),scanHdr.qform44(rownum,3),scanHdr.qform44(rownum,4)));
end

disp(sprintf('++++++++++++++++++++++++++ sform ++++++++++++++++++++++++++'));
for rownum = 1:4
  disp(sprintf('%f\t%f\t%f\t%f',scanHdr.sform44(rownum,1),scanHdr.sform44(rownum,2),scanHdr.sform44(rownum,3),scanHdr.sform44(rownum,4)));
end


