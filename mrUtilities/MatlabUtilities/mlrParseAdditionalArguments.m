function [arguments, nArgs] = mlrParseAdditionalArguments(argumentString, separator)
   
%parse string of arguments separated by separator and put them into a cell array of numerical and string arguments
%non-numerical values that are not between quotes are converted into strings
%
% Julien Besle, 08/07/2010
nArgs = 0;
arguments = cell(0);
remain = argumentString;
while ~isempty(remain)
  nArgs = nArgs+1;
  [token,remain] = strtok(remain, separator);
  try
    arguments{nArgs} = eval(token);
  catch
    arguments{nArgs} = token;
  end
      
end