% copyNiftiFile.m
%
%        $Id$ 
%      usage: success = copyNiftiFile(fromFilename,toFilename,<makeLink>)
%         by: justin gardner
%       date: 01/07/09
%    purpose: copy a nift file. Handles copying both .hdr and .img files
%             checks for file existence. If makeLink is set to 1, will
%             link the files rather than copy them. makeLink set to 2 will make a hard link.
%             If there is an associated .mat file (i.e. same name) that will be copied as well
%             Set overwrite to 1 to overwrite existing files without asking
%
function success = copyNiftiFile(fromFilename,toFilename,makeLink,overwrite)

% set initial return value
success = 0;

% check arguments
if ~any(nargin == [2 3 4])
  help copyNiftiFile
  return
end

% default to copy
if ieNotDefined('makeLink'),makeLink = 0;end
if ieNotDefined('overwrite'),overwrite = 0;end

% get calling function name
[st,i] = dbstack;
callingFunction = stripext(st(min(i+1,length(st))).file);

% check to see what extensions we should be expecting
if any(strcmp(getext(fromFilename),{'img','hdr'})) 
  extensions = {'hdr','img'};
else
  extensions = {getext(fromFilename)};
end

% check for an associated mat file
if mlrIsFile([stripext(fromFilename) '.mat'])
  extensions{end+1} = 'mat';
end

% check to make sure all files exist
for extensionNum = 1:length(extensions)
  ext = extensions{extensionNum};
  filename = sprintf('%s.%s',stripext(fromFilename),ext);
  if ~mlrIsFile(filename)
    mrWarnDlg(sprintf('(%s) Could not find nifti file %s',callingFunction,fromFilename));
    return
  end
end
  
% now copy files
for extensionNum = 1:length(extensions)
  ext = extensions{extensionNum};
  thisFromFilename = sprintf('%s.%s',stripext(fromFilename),ext);
  thisToFilename = sprintf('%s.%s',stripext(toFilename),ext);
  % check if toFile exists
  if overwrite
    r = inf;
  else
    r = 0;
  end
  if (extensionNum == 1) && mlrIsFile(thisToFilename)
    if ~isinf(r)
      r = askuser(sprintf('(%s) File %s already exists, overwrite',callingFunction,getLastDir(toFilename)),1);
    end
    if r==0
      return
    end
  end
  if ~isequal(thisFromFilename,thisToFilename)
    if makeLink
      disp(sprintf('(%s) Linking file %s to %s',callingFunction,thisFromFilename,thisToFilename));
      linkFile(thisFromFilename,thisToFilename,makeLink);
    else
      disp(sprintf('(%s) Copying file %s to %s',callingFunction,thisFromFilename,thisToFilename));
      copyfile(thisFromFilename,thisToFilename);
    end
  end
end
success = 1;

