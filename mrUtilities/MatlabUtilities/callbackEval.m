function varargout = callbackEval(callback,firstArgs,lastArgs)

if ~ieNotDefined('firstArgs')
  arguments = firstArgs;
else
  arguments = [];
end

if iscell(callback)
  passedArgs = callback(2:end);
  callback = callback{1};
  for iArg = 1:length(passedArgs)
    arguments{end+1} = passedArgs{iArg};
  end 
end

if ~ieNotDefined('lastArgs')
  arguments = [arguments lastArgs];
end

% create the string to call the function
if nargout>=1
  funcall = '[';
  for i=1:nargout
    funcall = sprintf('%svarargout{%i},',funcall,i);
  end
  funcall(end:end+1)= ']=';
else
  funcall = '';
end
funcall = [funcall 'feval(callback'];
for i = 1:length(arguments)
  funcall = sprintf('%s,arguments{%i}',funcall,i);
end
funcall = sprintf('%s);',funcall);

% and call it
eval(funcall);