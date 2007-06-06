% evalargs.m
%
%      usage: evalargs(varargin,<alwaysUseNextArg>)
%         by: justin gardner
%       date: 12/13/05
%    purpose: passed in varargin, returns a string
%             that once evaluated sets the variables
%             called for. Thus, allows arguments like:
%
%             fun('var1','var2=3','var3',[3 4 5]);
%             will set var1=1, var2=3 and var3=[3 4 5];
%
%             should be run in the folling way:
%
%             function fun(varargin)
%             eval(evalargs(varargin));
%
%             if alwaysUseNextArg is set, then instead of
%             interpreting
%             fun('var1','str')
%               var1=1 str=1
%             will do
%               var1='str'
%
%             when eval is called, this will print out
%             what the variables are being set to. If
%             you want it to run quietly, do:
%             eval(evalargs({'gVerbose=0','var1',...}));
function evalstr = evalargs(args,alwaysUseNextArg)

if ~any(nargin == [1 2])
  help evalargs;
  return
end

if ~exist('alwaysUseNextArg','var'),alwaysUseNextArg=1;end

evalstr = 'global gVerbose;oldgVerbose=gVerbose;gVerbose=1;';
% check arguments in
skipnext = 0;
for i = 1:length(args)
  % skip if called for
  if skipnext
    skipnext = 0;
    continue
  end
  % evaluate anything that has an equal sign in it
  if isstr(args{i}) && ~isempty(strfind(args{i},'='))
    % if the argument is a numeric, than just set it
    % also if the argument happens to be a function, don't
    % even do str2num on that since it gives an annoyinw warning
    if ((exist(args{i}(strfind(args{i},'=')+1:end)) ~= 2) && ...
	~strcmp(args{i}(strfind(args{i},'=')+1:end),'string') && ...
	~isempty(str2num(args{i}(strfind(args{i},'=')+1:end))))
      evalstr = sprintf('%s%s;',evalstr,args{i});
    % same for a quoted string
    elseif args{i}(strfind(args{i},'=')+1)==''''
      evalstr = sprintf('%s%s;',evalstr,args{i});
    % otherwise, we got a unquoted string, so we need to set the quotes
    else      
      evalstr = sprintf('%s%s''%s'';',evalstr,args{i}(1:strfind(args{i},'=')),args{i}(strfind(args{i},'=')+1:end));
    end
    % any quote needs to be two single quotes
    args{i}(strfind(args{i},''''))='"';
    evalstr = sprintf('%sif gVerbose,disp(sprintf(''setting: %s''));,end,',evalstr,args{i});
  % if it is not evaluated then either it means to set the variable
  % or to set the variable to the next argument, we determine this
  % by whether the next argument is a string or not. If it is not
  % a string then it means to set the variable to that argument
  elseif isstr(args{i})
    if (length(args) >= (i+1)) && (~isstr(args{i+1}) || alwaysUseNextArg)
      % set the variable to the next argument
      if ~isstr(args{i+1})
	evalstr = sprintf('%s%s=varargin{%i};',evalstr,args{i},i+1);
	evalstr = sprintf('%sif gVerbose,disp(sprintf(''setting: %s=varargin{%i}''));,end,',evalstr,args{i},i+1);
      else
	evalstr = sprintf('%s%s=''%s'';',evalstr,args{i},args{i+1});
	evalstr = sprintf('%sif gVerbose,disp(sprintf(''setting: %s=''''%s''''''));,end,',evalstr,args{i},args{i+1});
      end
      skipnext = 1;
    else
      % just set the variable to one, since the next argument
      % does not contain a non string
      evalstr = sprintf('%s%s=1;',evalstr,args{i});
      evalstr = sprintf('%sif gVerbose,disp(sprintf(''setting: %s=1''));,end,',evalstr,args{i});
    end
  else
    % skip anythin we don't know about
    if ~skipnext
    else
      skipnext = 0;
    end
  end
end

evalstr = sprintf('%sgVerbose=oldgVerbose;',evalstr);


