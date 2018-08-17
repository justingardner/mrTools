% mlrWriteJSON.m
%
%        $Id:$ 
%      usage: str = mlrWriteJSON(obj)
%         by: justin gardner
%       date: 07/15/18
%    purpose: Function to write simple json files,
%             so that we don't have to rely on matla2018
%
%       e.g.: 
%
%obj.duh = [];
%obj.uhm = [3 4 5];
%obj.yowsa = 'hello';
%obj.n = 675;
%obj.huh.duh = 3;
%obj.huh.uhm = [4 5 7];
%
%jsonStr = mlrWriteJSON(obj);
%
function str = mlrWriteJSON(obj,str)

% check arguments
if ~any(nargin == [1 2])
  help mlrWriteJSON
  return
end

% start str
if nargin == 1,str = '';end

% standard indent
indentStr = '    ';

% start file
str = sprintf('%s{\n',str);

% now go through obect fields
f = fieldnames(obj);
for iField = 1:length(f)
  % write out field name
  str = sprintf('%s%s\"%s\": ',str,indentStr,f{iField});
  % now write values
  if isempty(obj.(f{iField}))
    % null value
    str = sprintf('%s%s',str,'null');
  % numeric values
  elseif isnumeric(obj.(f{iField}))
    % if it is of length one
    if length(obj.(f{iField})) == 1
      % print number
      str = sprintf('%s%f',str,obj.(f{iField}));
    else
      % an array
      str = sprintf('%s[ ',str);
      for iArray = 1:length(obj.(f{iField}))
	str = sprintf('%s%f',str,obj.(f{iField})(iArray));
	if iArray < length(obj.(f{iField}))
	  % add a comma
	  str = sprintf('%s, ',str);
	end
      end
      str = sprintf('%s ]',str);
    end
  % string values
  elseif isstr(obj.(f{iField}))
    % print the string
    str = sprintf('%s\"%s\"',str,obj.(f{iField}));
  % struct
  elseif isstruct(obj.(f{iField}))
    % call recursively
    str = mlrWriteJSON(obj.(f{iField}),str);
    if isempty(str),return,end
  % unknown type
  else
    disp(sprintf('(mlrWriteJSON) Aborting. Do not know how to write out field: %s',f));
    str = [];
    return
  end
  % end of line
  if iField < length(f)
    str = sprintf('%s,\n',str);
  else
    str = sprintf('%s\n',str);
  end
end

% end file
str = sprintf('%s}\n',str);



