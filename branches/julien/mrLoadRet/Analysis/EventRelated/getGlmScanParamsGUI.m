% getScanParamsGUI.m
%
%        $Id$
%      usage: params = getScanParamsGUI(thisView,params,useDefault)
%         by: julien besle, extracted from eventRelatedGLMGUI so taht it can be called from eventRelatedGUI eventually
%       date: 23/03/2010
%    purpose: return scan specific parameters for GLM analysis
%

function [scanParams, params] = getGlmScanParamsGUI(thisView,params,useDefault)

keepAsking = 1;
groupNum = viewGet(thisView,'groupNum',params.groupName);
nScans = viewGet(thisView,'nScans',groupNum);
if isfield(params,'scanParams') && length(params.scanParams)==nScans
   scanParams = params.scanParams;
else
   % make the output as long as the number of scans
   scanParams = cell(1,nScans);
end

while keepAsking
  % check for stimfile, and if it is mgl/type then ask the
  % user which variable they want to do the anlysis on
  for scanNum = 1:nScans
  if ~ismember(scanNum,params.scanNum)
    scanParams{scanNum} = [];
  else
   % get scan info and description
    tr = viewGet(thisView,'framePeriod',scanNum,groupNum);
    if ~isfield(scanParams{scanNum},'scanInfo')
       scanParams{scanNum}.scanInfo = sprintf('%i: %s',scanNum,viewGet(thisView,'description',scanNum,groupNum));
    end
    if ~isfield(scanParams{scanNum},'description')
       scanParams{scanNum}.description = sprintf('Event related analysis of %s: %i',params.groupName,scanNum);
    end
    if ~isfield(scanParams{scanNum},'preprocess')
       scanParams{scanNum}.preprocess = '';
    end
    if ~isfield(scanParams{scanNum},'subsetBox') || ~strcmp(params.analysisVolume,'Subset box')
       scanDims = viewGet(thisView,'dims',scanNum,groupNum);
       scanParams{scanNum}.subsetBox = ['[1 ' num2str(scanDims(1)) ';1 ' num2str(scanDims(2)) ';1 ' num2str(scanDims(3)) ']'];
    end
    if ~isfield(scanParams{scanNum},'forceStimOnSampleOnset')
       scanParams{scanNum}.forceStimOnSampleOnset = 1;
    end
    if ~isfield(scanParams{scanNum},'estimationSupersampling') || isempty(scanParams{scanNum}.estimationSupersampling)
       scanParams{scanNum}.estimationSupersampling = 1;
    end
    if ~isfield(scanParams{scanNum},'acquisitionSubsample') || isempty(scanParams{scanNum}.acquisitionSubsample)
       scanParams{scanNum}.acquisitionSubsample = 1;
    end

 % Standard parameters to set
    
    subsetBoxVisibleOption = 'visible=1';
    subsetBoxEditableOption = 'editable=1';
    if strcmp(params.analysisVolume,'Loaded ROI(s)')
      subsetBoxVisibleOption = 'visible=0';
    elseif strcmp(params.analysisVolume,'Whole volume')
      subsetBoxEditableOption = 'editable=0';
    end

    taskVarParams = {...
      {'scan',scanParams{scanNum}.scanInfo,'type=statictext','Description of scan to set parameters for (not editable)'},...
      {'description',scanParams{scanNum}.description,'Event related analysis of [x...x]','Description of the analysis'}...
      {'preprocess',scanParams{scanNum}.preprocess,'String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'}...
      {'subsetBox', scanParams{scanNum}.subsetBox, subsetBoxVisibleOption, subsetBoxEditableOption, 'subset of voxels,  the form [X1 X2;Y1 Y2;Z1 Z2] (Zs are optional)'};
         };

  % Timing parameters
    % make sure we are running on a set with a stimfile
    stimfile = viewGet(thisView,'stimfile',scanNum,groupNum);

    if isempty(stimfile)
     mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',scanNum,params.groupName));
     scanParams = [];
     return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see if we have a stimfile from mgl, in which case we should
    % ask the user what the variable name is that they want ot use for the analysis
    if strfind(stimfile{1}.filetype,'mgl')

      % check to see what style this is, if the task variable does
      % not have a segmentTrace then it mus be an old style, in which
      % we used channels
      task = cellArray(stimfile{1}.task,2);
      if isfield(stimfile{1}.myscreen,'traces') && ~isfield(task{1}{1},'segmentTrace')
        % this is the old style, get the stimtrace number
        taskVarParams{end+1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]','incdecType=plusMinus',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
      else
        if exist('getTaskVarnames') ~= 2
          mrErrorDlg('(eventRelatedGUI) MGL function getTaskVarnames is not in path. You must have mgl in the path to extract stimulus timing information from an mgl stim file');
        end
        % this is the new tyle, ask for a variable name
        [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
        % if there is more than one task, then ask the user for that
        task = cellArray(stimfile{1}.task,2);
        if length(task)>1
          taskVarParams{end+1} = {'taskNum',num2cell(1:length(task)),'The task you want to use'};
        end
        % if there are multiple phases, then ask for that
        maxPhaseNum = 0;
        maxSegNum = 0;
        for tnum = 1:length(task)
          phaseNum{tnum} = num2cell(1:length(task{tnum}));
          maxPhaseNum = max(maxPhaseNum,length(task{tnum}));
          % if there are multiple _segments_, then ask for that
          for pnum = 1:length(task{tnum})
            segNum{tnum}{pnum} = num2cell(1:length(task{tnum}{pnum}.segmin));
            maxSegNum = max(maxSegNum,length(segNum{tnum}{pnum}));
          end
        end
        if maxPhaseNum > 1
          if length(task) == 1
            taskVarParams{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
          else
            taskVarParams{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
          end
        end

         % if there is more than one segment in any of the phases, ask the user to specify
         % should add some error checking.
        if maxSegNum > 1
          taskVarParams{end+1} = {'segmentNum',1,'The segment of the task you want to use','incdec=[-1 1]','incdecType=plusMinus'};
        end

         % set up to get the variable name from the user
        taskVarParams{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)};
      end
    elseif strfind(stimfile{1}.filetype,'eventtimes')  && ~isfield(scanParams{scanNum},'stimDuration') && isfield(stimfile{1}.mylog,'stimdurations_s')
          scanParams{scanNum}.stimDuration = 'fromFile';
    end

    if ~isfield(scanParams{scanNum},'stimDuration') || strcmp(params.hrfModel,'hrfDeconvolution')
       scanParams{scanNum}.stimDuration = tr;
    end

    taskVarParams{end+1} = {'stimDuration', scanParams{scanNum}.stimDuration, 'duration of stimulation/event (seconds, min=0.01s), a boxcar function that is convolved with hrf. If using deconvolution, should be equal to the frame period'};
    taskVarParams{end+1} = {'forceStimOnSampleOnset', scanParams{scanNum}.forceStimOnSampleOnset, 'type=checkbox','Forces stimulus onset to coincide with (sub)sample onsets'};
    taskVarParams{end+1} = {'estimationSupersampling',scanParams{scanNum}.estimationSupersampling,'incdec=[-1 1]','incdecType=plusMinus','minmax=[1 inf]','Supersampling factor of the HRF model. Set this to more than one in order to resolve the estimated HDR at a temporal resolution that is less than the frame rate. This is only required if both the design and the acquisition have been designed to achieve subsample HDR estimation'};
    taskVarParams{end+1} = {'acquisitionSubsample',scanParams{scanNum}.acquisitionSubsample, 'incdec=[-1 1]','incdecType=plusMinus','minmax=[1 inf]','If the subsample estimation factor is more than 1, specifies at which subsample of the frame period the signal is actually acquired.'};


    % give the option to use the same variable for all
    if (scanNum == params.scanNum(1)) && (length(params.scanNum)>1)
     taskVarParams{end+1} = {'sameForAll',1,'type=checkbox','Use the same parameters for all scans'};
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    % now we have all the dialog information, ask the user to set parameters
    if useDefault
       scanParams{scanNum} = mrParamsDefault(taskVarParams);
    else
       scanParams{scanNum} = mrParamsDialog(taskVarParams,'Set Scan Parameters');
    end

    % user hit cancel
    if isempty(scanParams{scanNum})
       scanParams = [];
       return
    end

    %get the number of events after running the pre-processing function
    d = loadScan(thisView, scanNum, [], 0);
    d = getStimvol(d,scanParams{scanNum});
    disp(sprintf('Getting number of conditions from scan %d', scanNum)); 
    subsetBox = eval(scanParams{scanNum}.subsetBox);
    
    % Try the pre-processing function
    preProcessFailure = 0;
    if ~isempty(scanParams{scanNum}.preprocess)
      [d, preProcessFailure] = eventRelatedPreProcess(d,scanParams{scanNum}.preprocess);
    end

    %various controls and variables settings
    if scanParams{scanNum}.estimationSupersampling<scanParams{scanNum}.acquisitionSubsample
      mrWarnDlg('(getScanParamsGUI) The acquisition subsample must be less than the subsample estimation factor','Yes');
      scanParams{scanNum} = [];
    elseif params.covCorrection && any(diff(subsetBox(1:2,:),1,2)<=params.covEstimationAreaSize)
      mrWarnDlg('(getScanParamsGUI) subset box size is too small for covariance estimation area size','Yes');
      scanParams{scanNum} = [];
    elseif (ischar(scanParams{scanNum}.stimDuration) || scanParams{scanNum}.stimDuration ~=tr) ...
        && strcmp(params.hrfModel,'hrfDeconvolution')...
        && (params.computeTtests || params.numberFtests)
      mrWarnDlg('(getScanParamsGUI) subTR sampling is not (yet) compatible with statistics on deconvolution weights','Yes');
      scanParams{scanNum} = [];
    elseif preProcessFailure
      mrWarnDlg(['(getScanParamsGUI) There was a problem running pre-processing function ' scanParams{scanNum}.preprocess],'Yes');
      scanParams{scanNum} = [];
    else
      keepAsking = 0;
      params.numberEvents = length(d.stimvol);
      params.stimNames = d.stimNames;
      disp(sprintf('%d conditions found', params.numberEvents));

      % check if the varname is a cell array, then convert to a cell array
      % instead of a string this is so that the user can specify a variable
      % name like {{'varname'}}
      if (isfield(scanParams{scanNum},'varname') &&...
         ichar(scanParams{scanNum}.varname) && ...
         (length(scanParams{scanNum}.varname) > 1) && ...
         (scanParams{scanNum}.varname(1) == '{'))
         scanParams{scanNum}.varname = eval(scanParams{scanNum}.varname);
      end

      % if sameForAll is set, copy all parameters into all scans and break out of loop
      if isfield(scanParams{scanNum},'sameForAll') && ...
         scanParams{scanNum}.sameForAll
         for i = 2:length(params.scanNum)
            % set the other scans params to the same as this one
            scanParams{params.scanNum(i)} = scanParams{params.scanNum(1)};
            % change the description field appropriately for this scan num
            description = scanParams{params.scanNum(1)}.description;
            groupNameLoc = strfind(description,params.groupName);
            if ~isempty(groupNameLoc)
               description = sprintf('%s%s: %i',description(1:groupNameLoc(1)),params.groupName,params.scanNum(i));
            end
            scanParams{params.scanNum(i)}.description = description;
         end
         break
      end
  %          taskVarParams = {};
    end
    if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
      scanParams = [];
      return;
    end
  end
  end
end