% eventRelatedPreProcess.m
%
%      usage: d = eventRelatedPreProcess(d)
%         by: justin gardner
%       date: 03/29/07
%    purpose: make various changes to the d structure before
%             passing it on to r2 processing, this allows for
%             little modifications that might be necessary to
%             fix small things
%
function d = eventRelatedPreProcess(d,type)

% check arguments
if ~any(nargin == [2])
  help eventRelatedPreProcess
  return
end

% move all the stimvolumes by some amount
stimvolShift  = getfieldnum(type,'stimvolShift');
if ~isempty(stimvolShift)
  disp(sprintf('Shifting by %i volumes',stimvolShift));
  for stimnum = 1:length(d.stimvol)
    d.stimvol{stimnum} = d.stimvol{stimnum}+stimvolShift;
  end
end

% length of tr
tr = getfieldnum(type,'tr');
if ~isempty(tr)
  disp(sprintf('Setting TR to %0.2f',tr));
  d.tr = tr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if the type has something like fieldname=num
% it returns the num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldval = getfieldnum(type,fieldname)

fieldval = getfieldstr(type,fieldname);
fieldval = str2num(fieldval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if the type has something like fieldname=string
% it returns the string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldval = getfieldstr(type,fieldname)

% see if it is there
if ~strfind(type,sprintf('%s=',fieldname))
  fieldval = [];
end

% get the field
fieldval = type(first(strfind(type,sprintf('%s=',fieldname))):length(type));
% up to the end
fieldval = fieldval(first(strfind(fieldval,'='))+1:length(fieldval));
% remove everything after first space
if (strfind(fieldval,' '))
  fieldval = fieldval(1:first(strfind(fieldval,' '))-1);
end


