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
% 
%             to print info on groups with scanSforms
%             groupInfo('Raw',1)
function retval = groupInfo(groupNum,verbose)

% check arguments
if ~any(nargin == [0 1 2])
  help groupInfo
  return
end

if ieNotDefined('verbose'),verbose = 0;end
view = newView('Volume');

% check home dir
homeDir = viewGet(view,'homeDir');
if ~isequal(homeDir,pwd)
  disp(sprintf('(groupInfo) Current directory (%s) is not the home directory of the session',getLastDir(pwd)));
  deleteView(view);
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% groupInfo with no arguments
%%%%%%%%%%%%%%%%%%%%%%%%%
% if groupNum was not given then display info about all groups
if ieNotDefined('groupNum')
  % get session info
  subject = viewGet(view,'subject');
  description = viewGet(view,'sessionDescription');
  magnet = viewGet(view,'magnet');
  operator = viewGet(view,'operator');
  coil = viewGet(view,'coil');
  protocol = viewGet(view,'protocol');
  % display session info
  disp(sprintf('%s',repmat('=',1,40)));
  disp(sprintf('homeDir: %s',homeDir));
  disp(sprintf('description: %s',description));
  disp(sprintf('operator: %s subject: %s',operator,subject));
  disp(sprintf('magnet: %s coil: %s protocol: %s',magnet,coil,protocol));
  disp(sprintf('%s',repmat('=',1,40)));
  % display each group
  for g = 1:viewGet(view,'nGroups')
    % get group info
    groupName = viewGet(view,'groupName',g);
    numScans = viewGet(view,'nScans',g);
    % use du to get disk usage
    [status result] = system(sprintf('du -k -d 0 %s',viewGet(view,'datadir',g)));
    dirSize = str2num(strtok(result));
    if (dirSize > 1000000)
      dirSize = sprintf('%0.1fG',dirSize/1000000);
    elseif (dirSize > 1000)
      dirSize = sprintf('%0.1fM',dirSize/1000);
    else
      dirSize = sprintf('%iK',dirSize);
    end
    % display group info
    disp(sprintf('%i: %s (%i scans) %s',g,groupName,numScans,dirSize));
  end
  deleteView(view);
  return
end

% if groupNum is a string, then user passed in a name rather than
% a number
if isstr(groupNum)
  groupName = groupNum;
  groupNum = viewGet(view,'groupNum',groupName);
  if isempty(groupNum)
    disp(sprintf('(groupInfo): No group %s',groupName));
    deleteView(view);
    return
  end
end

if (groupNum < 1) | (groupNum > viewGet(view,'numberOfGroups'))
  disp(sprintf('(groupInfo): No group number %i',groupNum));
  deleteView(view);
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

  % if verbose is set to 1 or above, then show scan transform
  if verbose >= 1 
    scanSform = viewGet(view,'scanSform',s,groupNum);
    disp(sprintf('   scanSform:'));
    for i = 1:4
      disp(sprintf('   %0.3f %0.3f %0.3f %0.3f',scanSform(i,1),scanSform(i,2),scanSform(i,3),scanSform(i,4)));
    end
  end
end

deleteView(view);