% disppercent.m
%
%      usage: disppercent(percentdone)
%         by: justin gardner
%       date: 10/05/04
%    purpose: display percent done
%             call with -inf for init
%             call with percent to update
%             call with inf to end
%       e.g.:disppercent(-inf,'Doing stuff');
%            disppercent(.5);
%            disppercent(inf);
%
function retval = disppercent(percentdone,mesg)

if ~isunix
  if percentdone == -inf
    disp(mesg);
  end
  return
end

% check command line arguments
if ((nargin ~= 1) && (nargin ~= 2))
  help disppercent;
  return
end

% global to turn off printing
global gVerbose;
if ~gVerbose
  return
end

global gDisppercent;

% systems without mrDisp (print w/out return that flushes buffers)
if exist('mrDisp', 'file') ~= 3
  if (percentdone == -inf) && (nargin == 2)
    disp(mesg);
  end
  return
end

% if this is an init then remember time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (percentdone == -inf)
  % set starting time
  gDisppercent.t0 = clock;
  % default to no message
  if (nargin < 2)
    mrDisp(sprintf('00%% (00:00:00)'));
  else
    mrDisp(sprintf('%s 00%% (00:00:00)',mesg));
  end    
% display total time at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (percentdone == inf)
  % get elapsed time
  numsecs = etime(clock,gDisppercent.t0);
  % separate seconds and milliseconds
  numms = round(1000*(numsecs-floor(numsecs)));
  numsecs = floor(numsecs);
  % if over a minute then display minutes separately
  if numsecs>60
    nummin = floor(numsecs/60);
    numsecs = numsecs-nummin*60;
    % check if over an hour
    if nummin > 60
      numhours = floor(nummin/60);
      nummin = nummin-numhours*60;
      timestr = sprintf('%i hours %i min %i secs %i ms',numhours,nummin,numsecs,numms);
    else
      timestr = sprintf('%i min %i secs %i ms',nummin,numsecs,numms);
    end
  else
    timestr = sprintf('%i secs %i ms',numsecs,numms);
  end
  % display time string
  mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\btook %s\n',timestr));
% otherwise show update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  % avoid things that will end up dividing by 0
  if (percentdone >= 1)
    percentdone = .99;
  elseif (percentdone <= 0)
    percentdone = 0.01;
  end
  % display percent done and estimated time to end
  if (gDisppercent.percentdone ~= floor(100*percentdone))
    mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%02i%% (%s)',floor(100*percentdone),disptime(etime(clock,gDisppercent.t0)*(1/percentdone - 1))));
  end
end
% remember current percent done
gDisppercent.percentdone = floor(100*percentdone);

% display time
function retval = disptime(t)


hours = floor(t/(60*60));
minutes = floor((t-hours*60*60)/60);
seconds = floor(t-hours*60*60-minutes*60);

retval = sprintf('%02i:%02i:%02i',hours,minutes,seconds);
