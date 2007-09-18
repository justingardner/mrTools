% groupInfo.m
%
%      usage: groupInfo(groupNum)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: groupInfo(1);
%             groupInfo('Raw');
%
%             to get info on all groups:
%             groupInfo
function retval = groupInfo(groupNum,verbose)

% check arguments
if ~any(nargin == [0 1 2])
  help groupInfo
  return
end

if ieNotDefined('verbose'),verbose = 1;end
view = newView('Volume');

%%%%%%%%%%%%%%%%%%%%%%%%%
% get the correct group
%%%%%%%%%%%%%%%%%%%%%%%%%
% if groupNum was not given then display info about all groups
if ieNotDefined('groupNum')
  subject = viewGet(view,'subject');
  description = viewGet(view,'sessionDescription');
  magnet = viewGet(view,'magnet');
  operator = viewGet(view,'operator');
  coil = viewGet(view,'coil');
  protocol = viewGet(view,'protocol');
  homeDir = viewGet(view,'homeDir');
  disp(sprintf('%s',repmat('=',1,40)));
  disp(sprintf('homeDir: %s',homeDir));
  disp(sprintf('description: %s',description));
  disp(sprintf('operator: %s subject: %s',operator,subject));
  disp(sprintf('magnet: %s coil: %s protocol %s',magnet,coil,protocol));
  disp(sprintf('%s',repmat('=',1,40)));
  for g = 1:viewGet(view,'nGroups')
    groupName = viewGet(view,'groupName',g);
    numScans = viewGet(view,'nScans',g);
    disp(sprintf('%i: %s (%i scans)',g,groupName,numScans));
  end
  return
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
  scanDims = viewGet(view,'scanDims',s,groupNum);
  totalJunkedFrames = viewGet(view,'totalJunkedFrames',s,groupNum);
  junkFrames = viewGet(view,'junkFrames',s,groupNum);
  
  % display info
  disp(sprintf('%i: %s',s,description));
  disp(sprintf('   Filename: %s GroupName: %s',filename,groupName));
  for i = 1:length(originalFilename)
    disp(sprintf('   Original Filename: %s Group: %s',originalFilename{i},originalGroupname{i}));
  end
  for i = 1:length(stimFilename)
    disp(sprintf('   StimFilename: %s',stimFilename{i}));
  end

  disp(sprintf('   junkFrames=[%s] totalJunkedFrames=[%s]',num2str(junkFrames),num2str(totalJunkedFrames)));
  disp(sprintf('   voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f Dims: [%i %i %i] Volumes=%i',scanVoxelSize(1),scanVoxelSize(2),scanVoxelSize(3),tr,scanDims(1),scanDims(2),scanDims(3),totalFrames));

  % if verbose is set above 1, then show scan transform
  if verbose >1 
    scanXform = viewGet(view,'scanXform',s,groupNum);
    disp(sprintf('   scanXform:'));
    for i = 1:4
      disp(sprintf('   %0.3f %0.3f %0.3f %0.3f',scanXform(i,1),scanXform(i,2),scanXform(i,3),scanXform(i,4)));
    end
  end
end

