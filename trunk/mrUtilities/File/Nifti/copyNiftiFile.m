% copyNiftiFile.m
%
%        $Id$ 
%      usage: success = copyNiftiFile(fromFilename,toFilename)
%         by: justin gardner
%       date: 01/07/09
%    purpose: copy a nift file. Handles copying both .hdr and .img files
%             checks for file existence.
%
function success = copyNiftiFile(fromFilename,toFilename)

success = 0;

% check arguments
if ~any(nargin == [2])
  help copyNiftiFile
  return
end

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
  copyfile(thisFromFilename,thisToFilename);
end
success = 1;

