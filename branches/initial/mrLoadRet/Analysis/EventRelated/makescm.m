% makescm.m
%
%      usage: makescm(d,hdrlen)
%         by: justin gardner
%       date: 07/28/04
%       e.g.: makescm(d)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it. if d.hdrlen is set then you can
%             just pass in the data structure
%
function d = makescm(d,hdrlen)

if (nargin == 1)
  if (isfield(d,'hdrlen'))
    hdrlen = d.hdrlen;
  else
    hdrlen = 25;
  end
elseif (nargin ~= 2)
  help makescm;
  return
end

% if run transitions are not set that means this 
% is just one experiment being computed. If it
% is set then it means that we have concatenated
% a bunch of runs together, so the convolution
% matrix should be calculated for each run and
% then stacked on top of each other--this way
% transitions between runs are respected
if ~isfield(d,'names')
  % make stimulus convolution matrix
  stimcmatrix = [];
  d.nhdr = 0;
  % ignore first stimulus length cause that
  % is just the task pulse
  if ((d.pulselens(1) ~= 1) && (d.pulselens(1) < 50)),startstim = 2;,else,startstim = 1;,end
  for i = startstim:length(d.stimvol)
    % also, if there is only one pulse that
    % means it is just a start pulse, so ignore
    % that as well
    if (length(d.stimvol{i}) == 1)
      disp(sprintf('Ignoring pulse of length %0.2f sec (only one pulse at %0.2f sec found)',d.pulselens(i)/1000,d.stimtimes{i}/100));
    else
      stimarray = zeros(1,d.dim(4));
      stimarray(d.stimvol{i}) = 1;
      % shorten, to remove timeseries before first volume
      stimarray = stimarray(d.firstvol:d.dim(4));
      % stack stimcmatrices horizontally
      stimcmatrix = [stimcmatrix stimconv(stimarray,hdrlen)];
      % update number of hdrs
      d.nhdr = d.nhdr+1;
    d.stimpulselens(d.nhdr) = d.pulselens(i);
    end
  end
  d.scm = stimcmatrix;
  d.hdrlen = hdrlen;
  d.volumes = d.firstvol:d.dim(4);
else
  % cycle through each run, constructing the convolution matrix
  allstimcmatrix = [];
  for runnum = 1:d.names.n
    % make stimulus convolution matrix
    stimcmatrix = [];
    d.nhdr = 0;
    % ignore first stimulus length cause that
    % is just the task pulse
    if ((d.pulselens(1) ~= 1) && (d.pulselens(1) < 50)),startstim = 2;,else,startstim = 1;,end
    for i = startstim:length(d.stimvol)
      % also, if there is only one pulse that
      % means it is just a start pulse, so ignore
      % that as well
      if (length(d.stimvol{i}) == 1)
	disp(sprintf('Ignoring pulse of length %0.2f sec (only one pulse at %0.2f sec found)',d.pulselens(i)/1000,d.stimtimes{i}/100));
      else
	stimarray = zeros(1,d.dim(4));
	stimarray(d.stimvol{i}) = 1;
	% select out only the part of the stim array for this run
	% for the first run, throw out firstvols
	if (runnum == 1)
	  stimarray = stimarray(d.firstvol:d.names.runtransition(runnum,2));
	else
	  stimarray = stimarray(d.names.runtransition(runnum,1):d.names.runtransition(runnum,2));
	end
	% stack stimcmatrices horizontally
	stimcmatrix = [stimcmatrix stimconv(stimarray,hdrlen)];
	% update number of hdrs
	d.nhdr = d.nhdr+1;
	d.stimpulselens(d.nhdr) = d.pulselens(i);
      end
    end
    % stack this run's stimcmatrix on to the last one
    allstimcmatrix = [allstimcmatrix;stimcmatrix];
  end
  
  d.scm = allstimcmatrix;
  d.hdrlen = hdrlen;
  d.volumes = d.firstvol:d.dim(4);
end  

