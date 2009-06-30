% readprocpar.m
%
%      usage: readprocpar.m()
%         by: justin gardner
%       date: 07/03/03
%    purpose: reads procpar into matlab variable
%      usage: procpar = readprocpar;
%
function procpar = readprocpar(procdir,verbose)


% argument check
if (nargin == 0)
  procdir = './';
  verbose = 0;
elseif (nargin == 1)
  verbose = 0;
elseif (nargin ~= 2)
  help readprocpar;
  return
end
if (nargin >= 1)
  if (procdir(length(procdir)) ~= '/')
    procdir(length(procdir)+1) = '/';
  end
end

if (verbose)
  disp('(readprocpar) Reading procpar...');
end
procpar = [];

% open the procpar
fprocpar = fopen([procdir 'procpar'],'r');
if (fprocpar == -1)
  disp(sprintf('(readprocpar) ERROR: Could not open %s',[procdir 'procpar']));
  return
end

% scan through procpar looking for parameters needed for parsing
line = fgets(fprocpar);
while (line ~= -1) 
  % get the first token of the line
  token = strtok(line);
  % if it is a character and not doublequote then we have an array name
  if (~((token(1) >= '0') & (token(1) <= '9')) & (token(1) ~= '"'))
    fieldname = token;
    % get next line
    line = fgets(fprocpar);
    % first number is just length of array
    [len line] = strtok(line);
    len = str2num(len);
    % just an array of numbers
    if ((length(line) > 2) & (line(2) ~= '"'))
      eval(sprintf('procpar.%s = str2num(line);',fieldname));
      % confirm length
      arraylen = eval(sprintf('length(procpar.%s)',fieldname));
      if (arraylen ~= len)
	sprintf('WARNING: %s should be %i elements, but found only %i',fieldname,len,arraylen);
      end
    % character arrays
    else
      for j = 1:len
	% strip leading white space
	while (length(line) & (line(1) == ' '))
	  line = line(2:length(line));
	end
	token = line;
	% make sure we have a valid string token
	if (isempty(token) | (token(1) ~= '"'))
	  sprintf('WARNING: %s missing string element %i/%i',fieldname,j,len);
	  eval(sprintf('procpar.%s{%i} = '''';',fieldname,j));
	  j = len;
	else
	  str = strtok(line,sprintf('"\n'));
	  % convert single to double quotes
	  str(findstr(str,'''')) = '"';
	  eval(sprintf('procpar.%s{%i} = ''%s'';',fieldname,j,str));
	end
	line = fgets(fprocpar);
      end
    end
  end
  line = fgets(fprocpar);
end

fclose(fprocpar);

% for epi images, there are navigator echos, which
% should be subtracted from the number of lines.
% this can be known from the name of the petable
if (isfield(procpar,'petable'))
  token = procpar.petable{1};
  % the petable name should be something like
  % "epi132alt8k". We want the second number
  if (strncmp(token,'epi',3)) 
    j = 1;
    % go past any intial characters
    while((j < length(token)) && (isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % then skip the numbers
    while((j < length(token)) && (~isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % and skip the next characters
    while((j < length(token)) && (isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % find end of numbers
    k = j;
    while((k < length(token)) && (~isempty(strfind('0123456789',token(k))))),k=k+1;,end
    procpar.navechoes = str2num(token(j:k-1));
  end
end

if (verbose),disp('(readprocpar) Finished.');,end