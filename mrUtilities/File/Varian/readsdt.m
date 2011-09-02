% readsdt
%
%      usage: readsdt(filename,numvolumes)
%         by: justin gardner
%       date: 02/26/03
%    purpose: read sdt/spm format 
%    history:  11/26/2008 modified by Pei, add the endian part
%              endian:ieee-le or ieee-be don't have a space between
%              field and data
%              8/15/2011 jlg: reads edt files now. 
%
function d = readsdt(filename,numvolumes)

t0 = clock;

d = [];
% check arguments
if (nargin == 1)
  numvolumes = inf;
elseif (nargin ~= 2)  
  disp('USAGE: readstd(filename)');
  return
end

% get the extension
inputFilename = filename;
ext = getext(filename);
if ~any(strcmp({'sdt' 'spr' 'edt' 'epr'},ext)) 
  filename = setext(filename,'spr',0);
  if ~isfile(filename)
    filename = setext(filename,'epr');
    if ~isfile(filename)
      disp(sprintf('(readsdt) Could not find file %s',inputFilename));
      return
    end
  end
end
if ~isfile(filename)
  disp(sprintf('(readsdt) Could not find file %s',inputFilename));
  return
end

% see if this is an sdt or edt file
ext = getext(filename);
if any(strcmp(ext,{'sdt','spr'}))
  type = 'sdt';
  dataFilename = setext(filename,'sdt');
  paramFilename = setext(filename,'spr');
else
  type = 'edt';
  dataFilename = setext(filename,'edt');
  paramFilename = setext(filename,'epr');
end

% open the parameter file 
fspr = fopen(paramFilename,'r');
if (fspr == -1)
  disp(sprintf('(readsdt) Could not open parameter file %s',paramFilename));
  return
end

% cylce through line by line to get parameters
line = fgetl(fspr);
while (~isempty(line) & (line ~= -1))
  % find the colon
  colonpos = findstr(line,':');
  if (~isempty(colonpos))
    paramname = line(1:colonpos(1)-1);
    paramvalue = line(colonpos(1)+1:length(line));
    % see if it is a string or numeric values
    if (sum(double(paramvalue) > double('A')))
      % get rid of any spaces at beginning of paramvalue
      spacepos = strfind(paramvalue,' ');
      if (~isempty(spacepos))
	strstart = 1;
	while ((strstart <= length(spacepos)) && (spacepos(strstart) == strstart))
	  strstart = strstart+1;
	end
	paramvalue = paramvalue(strstart:length(paramvalue));
      end	
      eval(sprintf('d.%s = ''%s'';',paramname,paramvalue));
    else
      eval(sprintf('d.%s = [%s];',paramname,paramvalue));
    end
  end
  line = fgetl(fspr);
end
fclose(fspr);

% make sure we have necessary parameter values
needfields = {'dim','numDim'};
for i = 1:length(needfields)
  if (~isfield(d,needfields{i}))
    disp(sprintf('(readsdt) Missing %s field',needfields{i}));
    return
  end
end

if ~isfield(d,'dataType')
  if strcmp(type,'edt')
    d.dataType = 'COMPLEX';
  else
    disp(sprintf('(readsdt) !!! dataType not found in header, setting to float32'));
    d.dataType = 'REAL';
  end
end

% get wordisze
datatypes = strvcat('BYTE','WORD','LWORD','REAL','COMPLEX');
datasizes = [           1      2       4      4         4];
matlabtypes = {'int8','int16','int32','float32','float32'};
whichtype = strmatch(d.dataType,datatypes);
if (isempty(whichtype))
  disp(sprintf('(readsdt) Unrecogonized dataType %s',d.dataType));
  return
else
  d.wordsize = datasizes(whichtype);
end

% calculate data size
d.filelen = prod(d.dim)*d.wordsize;
if strcmp(d.dataType,'COMPLEX'),d.filelen = d.filelen*2;end

if numvolumes <= 0,return,end

% if we have the endian field, we read data in as it indicates
% if we do not have this field, read data as default (big endian)
if isfield(d, 'endian')
  if strcmp(d.endian,'ieee-le')
    fsdt = fopen(dataFilename,'r','ieee-le');
  elseif strcmp(d.endian,'ieee-be')
    fsdt = fopen(dataFilename,'r','ieee-be');
  end
else
  fsdt = fopen(dataFilename,'r', 'ieee-be');
end

% check file status
if (fsdt == -1)
  disp(sprintf('(readsdt) Could not open data file %s',dataFilename));
  return
end

% try to seek to the end
if (fseek(fsdt,0,'eof') == -1)
  disp(sprintf('(readsdt) Could not seek to end of %s',dataFilename));
  return
end
fpos = ftell(fsdt);
if (fpos == -1)
  disp(sprintf('(readsdt) Could not find position of end of %s',dataFilename));
  return
end

% see if we have enough bytes
if (fpos ~= d.filelen)
  disp(sprintf('(readsdt) %s has %i bytes, but should be %i',dataFilename,fpos,d.filelen));
  return
end  

% seek back to beginning
if (fseek(fsdt,0,'bof') == -1)
  disp(sprintf('(readsdt) Could not seek to beginning of %s',dataFilename));
  return
end

% reset the number of volumes if user calls for it
if ((numvolumes ~= inf) && (numvolumes < d.dim(d.numDim)))
  d.dim(d.numDim) = numvolumes;
end

% now read the data
% allocate memory for data
d.data = zeros(prod(d.dim),1);

% just read as one big chunk
if isfield(d, 'endian')
  if strcmp(d.endian,'ieee-le')
    [d.data, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-le');
  elseif strcmp(d.endian,'ieee-be')
    [d.data, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-be');
  end
else
  [d.data, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-be');
end
% make sure we have read enough bytes
if (total ~= prod(d.dim))
  disp(sprintf('(readsdt) Only read %i of %i data points',count,prod(d.dim)));
  return
end
% reshape, first generate command to reshape
command = 'd.data = reshape(d.data';
for i = 1:d.numDim
  command = strcat(command,sprintf(',%i',d.dim(i)));
end
command = strcat(command,');');
% then evaluate that command
eval(command);


% for complex data read the other part of data
if strcmp(d.dataType,'COMPLEX')
  if isfield(d, 'endian')
    if strcmp(d.endian,'ieee-le')
      [d.datai, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-le');
    elseif strcmp(d.endian,'ieee-be')
      [d.datai, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-be');
    end
  else
    [d.datai, total] = fread(fsdt,prod(d.dim),matlabtypes{whichtype}, 'ieee-be');
  end
  % make sure we have read enough bytes
  if (total ~= prod(d.dim))
    disp(sprintf('(readsdt) Only read %i of %i data points',count,prod(d.dim)));
    return
  end
  % reshape, first generate command to reshape
  command = 'd.datai = reshape(d.datai';
  for i = 1:d.numDim
    command = strcat(command,sprintf(',%i',d.dim(i)));
  end
  command = strcat(command,');');
  % then evaluate that command
  eval(command);
end

% close the file. we are done.
fclose(fsdt);

d.filename = filename;

% get full path of fid
if (isfield(d,'fidName'))
  d.fidName = [getpath(filename) getlastdir(d.fidName)];
else
  d.fidName = '';
end

% also write out the endian information
if ~isfield(d, 'endian')
  d.endian = 'ieee-be';
end

% convert to complex data
if strcmp(d.dataType,'COMPLEX')
  d.data = d.data + d.datai*sqrt(-1);
  d = rmfield(d,'datai');
end

disp(sprintf('(readsdt) Took %0.2f seconds',etime(clock,t0)));

