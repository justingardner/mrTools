% groupInfo.m
%
%      usage: groupInfo(groupNum)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: groupInfo(1);
%             groupInfo('Raw');
function retval = groupInfo(groupNum)

% check arguments
if ~any(nargin == [0 1])
  help groupInfo
  return
end

view = newView('Volume');

%%%%%%%%%%%%%%%%%%%%%%%%%
% get the correct group
%%%%%%%%%%%%%%%%%%%%%%%%%
% if groupNum was not given, then default to MotionComp then Raw
if ieNotDefined('groupNum')
  motionCompGroupNum = [];rawGroupNum = [];
  for g = 1:viewGet(view,'numberOfGroups')
    if strcmp('MotionComp',viewGet(view,'groupName',g))
      motionCompGroupNum = g;
    elseif strcmp('Raw',viewGet(view,'groupName',g))
      rawGroupNum = g;
    end
  end
  if ~isempty(motionCompGroupNum)
    groupNum = motionCompGroupNum;
  elseif ~isempty(rawGroupNum)
    groupNum = rawGroupNum;
  else
    disp(sprintf('(groupInfo): Could not find a Raw or MotionComp group'));
    return
  end
end

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

% now go through scans and print information
for s = 1:viewGet(view,'numberOfScans',groupNum)
  % grab info
  description = viewGet(view,'description',s,groupNum);
  scanVoxelSize = viewGet(view,'scanVoxelSize',s,groupNum);
  tr = viewGet(view,'framePeriod',s,groupNum);
  totalFrames = viewGet(view,'totalFrames',s,groupNum);
  filename = viewGet(view,'tSeriesFile',s,groupNum);
  originalFilename = viewGet(view,'originalFilename',s,groupNum);
  originalGroupname = viewGet(view,'originalGroupname',s,groupNum);
  stimFilename = viewGet(view,'stimFilename',s,groupNum);
  
  % display info
  disp(sprintf('%i: %s',s,description));
  disp(sprintf('   Filename: %s GroupName: %s',filename,groupName));
  for i = 1:length(originalFilename)
    disp(sprintf('   Original Filename: %s Group: %s',originalFilename{i},originalGroupname{i}));
  end
  for i = 1:length(stimFilename)
    disp(sprintf('   StimFilename: %s',stimFilename{i}));
  end
  
  disp(sprintf('   voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f Volumes=%i',scanVoxelSize(1),scanVoxelSize(2),scanVoxelSize(3),tr,totalFrames));
end
