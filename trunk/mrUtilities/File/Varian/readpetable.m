% readpetable.m
%
%      usage: readpetable(filename)
%         by: justin gardner
%       date: 08/11/11
%    purpose: 
%
function petable = readpetable(filename,varargin)

petable = [];
% check arguments
if nargin < 1
  help readpetable
  return
end

% petable
petablePath = [];
getArgs(varargin,{'petablePath=~/vnmrsys/tablib'});

% look for petable path
if ~isdir(petablePath)
  petablePath = '/usr4/justin/vnmrsys/tablib';
  if ~isdir(petablePath)
    disp(sprintf('(readpetable) Could not find petable directory: %s',petablePath));
    return
  end
end

% try to read petable file
if ~isfile(fullfile(petablePath,filename))
  disp(sprintf('(readpetable) Could not find petable: %s',fullfile(petablePath,filename)));
  return
end

% open petable
fPetable = fopen(fullfile(petablePath,filename),'r');

varname = textscan(fPetable,'%s =',1);
varval = textscan(fPetable,'%d');
while ~isempty(varname{1})
  % set the value in the petable variable
  petable.(varname{1}{1}) = double(varval{1});
  % and load the next part of the file
  varname = textscan(fPetable,'%s =',1);
  varval = textscan(fPetable,'%d');
end

% close petable
fclose(fPetable);
