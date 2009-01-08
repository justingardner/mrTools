% copyNiftiFile.m
%
%        $Id$ 
%      usage: success = copyNiftiFile(fromFilename,toFilename,<makeLink>)
%         by: justin gardner
%       date: 01/07/09
%    purpose: copy a nift file. Handles copying both .hdr and .img files
%             checks for file existence. If makeLink is set to 1, will
%             link the files rather than copy them.
%
function success = copyNiftiFile(fromFilename,toFilename,makeLink)

% set initial return value
success = 0;

% check arguments
if ~any(nargin == [2 3])
  help copyNiftiFile
  return
end

% default to copy
if ieNotDefined('makeLink'),makeLink = 0;end

% get calling function name
[st,i] = dbstack;
callingFunction = stripext(st(min(i+1,length(st))).file);

% check to see what extensions we should be expecting
if any(strcmp(getext(fromFilename),{'img','hdr'})) 
  extensions = {'hdr','img'};
else
  extensions = {getext(fromFilename)};
end
  
% check to make sure all files exist
for extensionNum = 1:length(extensions)
  ext = extensions{extensionNum};
  filename = sprintf('%s.%s',stripext(fromFilename),ext);
  if ~isfile(filename)
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
  r = 0;
  if (extensionNum == 1) && isfile(thisToFilename)
    if ~isinf(r)
      r = askuser(sprintf('(%s) File %s already exists, overwrite',callingFunction,getLastDir(toFilename)),1);
    end
    if r==0
      return
    end
  end
  disp(sprintf('(%s) Copying file %s to %s',callingFunction,thisFromFilename,thisToFilename));
  if makeLink
    linkFile(thisFromFilename,thisToFilename);
  else
    copyfile(thisFromFilename,thisToFilename);
  end
end
success = 1;

