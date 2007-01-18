% mrCleanDir.m
%
%      usage: mrCleanDir()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%
function retval = mrCleanDir()

% check arguments
if ~any(nargin == [0])
  help mrCleanDir
  return
end

view = newView('Volume');

% check for unlinked files
groups = viewGet(view,'groupNames');

for g = 1:length(groups)
  % set the current group
  groupNum = viewGet(view,'groupNum',groups{g});
  view = viewSet(view,'curGroup',groupNum);
  nScans = viewGet(view,'nScans');
  
  if nScans > 0
    % get the directory
    [tseriesDirName] = fileparts(viewGet(view,'tseriesPath',1));
    tseriesDir = dir(sprintf('%s/*.hdr',tseriesDirName));
    for i = 1:length(tseriesDir)
      tseriesDir(i).match = 0;
    end
  end
  
  % look for unmatched files
  for scanNum = 1:nScans
    tseriesFilename = viewGet(view,'tseriesPath',scanNum,groupNum);
    [thisDirName thisFilename] = fileparts(tseriesFilename);
    for filenum = 1:length(tseriesDir)
      [dirname,filename] = fileparts(tseriesDir(filenum).name);
      if strcmp(filename,thisFilename) && strcmp(tseriesDirName,thisDirName)
	tseriesDir(filenum).match = 1;
      end
    end
  end

  % count to see if we have all matches
  matched = 0;
  for i = 1:length(tseriesDir)
    if tseriesDir(i).match
      matched = matched+1;
    end
  end

  % if we have more files in directory than that are matched,...
  if (matched < length(tseriesDir))
    % display the names of the hdr/img/mat files that are not matched
    for i = 1:length(tseriesDir)
      if ~tseriesDir(i).match
	[path baseFilename] = fileparts(tseriesDir(i).name);
	filename = sprintf('%s/%s.hdr',tseriesDirName,baseFilename);
	if isfile(filename),disp(filename),end
	filename = sprintf('%s/%s.img',tseriesDirName,baseFilename);
	if isfile(filename),disp(filename),end
	filename = sprintf('%s/%s.mat',tseriesDirName,baseFilename);
	if isfile(filename),disp(filename),end
      end
    end
    % and ask user if they should be deleted
    if askuser(sprintf('Delete files from group %s',groups{g}))
      for i = 1:length(tseriesDir)
	if ~tseriesDir(i).match
	  [path baseFilename] = fileparts(tseriesDir(i).name);
	  filename = sprintf('%s/%s.hdr',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  filename = sprintf('%s/%s.img',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  filename = sprintf('%s/%s.mat',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  
	end
      end
    end
  else
    disp(sprintf('Group %s matches (%i:%i)',groups{g},length(tseriesDir),nScans));
  end

end
