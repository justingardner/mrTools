% mlrExportForAnalysis.m
%
%        $Id:$ 
%      usage: output = mlrExportForAnalysis(v)
%         by: justin gardner
%       date: 12/06/18
%    purpose: Export a structure for data analysis by pyton etc
%
function output = mlrExportForAnalysis(varargin)

% check arguments
if ~any(nargin == [2])
  help mlrExportForAnalysis
  return
end

% check if first argument is not view, that
% means that we have been called from the gui
% so we got an hObject
if ishandle(varargin{1})
  % get view
  v = viewGet(getfield(guidata(varargin{1}),'viewNum'),'view');
elseif isview(varargin{1})
  v = varargin{1};
else
  disp(sprintf('(mlrExportForAnalysis) Call with a view'));
  return
end

% get group name and scan
output.groupName = viewGet(v,'groupName');
output.scanNum = viewGet(v,'curScan');

% load the time series
disppercent(-inf,'(mlrExportForAnalysis) Loading time series');
output.tSeries = loadTSeries(v);
disppercent(inf);
output.dims = size(output.tSeries);

% check what type of analysis is here
analysis = viewGet(v,'analysis');
if ~isempty(analysis)
  if strcmp(analysis.type,'erAnal')
    output = getEventRelatedData(v,analysis,output);
  end
end

% save file
[filename,pathname] = uiputfile({'*.mat','Matlab file'},'Save as');
if isequal(filename,0) || isequal(pathname,0)
  return
else
  save(fullfile(pathname,filename),output,'-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getEventRelatedData    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getEventRelatedData(v,analysis,output)


% get data (the one in gEventRelatedPlot has some fields removed to save space)
% so re get it here.
output.eventRelated.d = analysis.d{output.scanNum};
concatInfo = viewGet(v,'concatInfo');
keyboard
% get stimvol
output.eventRelated.stimvol = output.eventRelated.d.stimvol;

% resort time series as time from stimvol
% first pull out some variables from structures
hdrlen = output.eventRelated.d.hdrlen;
stimvol = output.eventRelated.stimvol;
% get run transitions
if ~isempty(concatInfo) && isfield(concatInfo,'runTransition') && ~isempty(concatInfo.runTransition)
  runTransition = concatInfo.runTransition;
else
  runTransition = [1 size(output.tSeries,4)];
end

disppercent(-inf,'(mlrExportForAnalysis) Exporting tSeries sorted by events');
for iStimType = 1:length(stimvol)
  % start trace event time
  transitionNum = 1;
  % init
  output.eventRelated.tSeriesByStim{iStimType} = nan(output.dims(1),output.dims(2),output.dims(3),hdrlen);
  for iStim = 1:length(stimvol{iStimType})
    % get the start of the event
    eventStart = stimvol{iStimType}(iStim);
    % check transition
    if eventStart > runTransition(transitionNum,2)
      % if we are passed a scan boundary go to next scan
      transitionNum = min(size(runTransition,1),transitionNum+1);
    end
    % get end of event 
    eventEnd = min(eventStart + hdrlen-1,runTransition(transitionNum,2));
    % now get the relevant part of the time series
    eventTSeries = output.tSeries(:,:,:,eventStart:eventEnd);
    % pack to end with nan if we are missing data
    eventTSeries(:,:,:,end+1:hdrlen) = nan;
    % and stick in the field for voxel
    output.eventRelated.tSeriesByStim{iStimType}(:,:,:,iStim,1:hdrlen) = eventTSeries;
    disppercent(calcPercentDone(iStimType,length(stimvol),iStim,length(stimvol{iStimType})));
  end
end
disppercent(inf);

