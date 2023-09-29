% pRFGetSomatoStimImageFromStimfile
%
%        $Id:$ 
%      usage: stim = pRFGetSomatoStimImageFromStimfile(stimfile,<timePoints>)
%
%         by: ds - based mostly on code by justin gardner
%       date: 2016/02
%    purpose: Pass in a stimfile (can be either a string filename, or a strucutre
%             with myscreen/task) created with mgl / task code 

%             implementation here is based on how information is obtained from 
%             mglRetinotopy stimfile.
%
%             Will
%             create a volume of dimensions x,y,t with the stimulus image (load
%             in mlrVol to view). stim.x and stim.y are the X and Y coordinates
%             (units??). stim.t is the array of times at which image is taken.
%
%             Optionally arguments:
%
%             timePoints: array for which the stim image should be computed. 
%
%             Note for developers - this function needs to keep up-to-date with
%             any changes in the display loop of mglRetinotopy to interpret
%             the stimfiles correctly
%
function stim = pRFGetSomatoStimImageFromStimfile(stimfile,varargin)

% set default return arguments
stim = [];

% check arguments
if nargin < 1
  help pRFGetSomatoStimImageFromStimfile
  return
end

% parse arguments
timePoints = [];screenWidth = [];screenHeight = [];volTrigRatio = [];
xFlip = [];yFlip = [];timeShift = [];verbose = [];
getArgs(varargin,{'timePoints=[]','screenWidth=[]','screenHeight=[]','volTrigRatio=[]','xFlip=0','yFlip=0','timeShift=0','verbose=1','saveStimImage=0','recomputeStimImage=0'});

% handle cell array
if iscell(stimfile) && ((length(stimfile)>1) || (length(stimfile{1})>1))
  for i = 1:length(stimfile)
    % get current volTrigRatio
    if isempty(volTrigRatio)
      thisVolTrigRatio = [];
    else
      thisVolTrigRatio = volTrigRatio{i};
    end
    stim{i} = pRFGetSomatoStimImageFromStimfile(stimfile{i},'timePoints',timePoints,'screenWidth',screenWidth,'screenHeight',screenHeight,'volTrigRatio',thisVolTrigRatio,'xFlip',xFlip,'yFlip',yFlip,'timeShift',timeShift,'verbose',verbose,'saveStimImage',saveStimImage,'recomputeStimImage',recomputeStimImage);
    if isempty(stim{i}),stim = [];return;end
  end
  return
end

% check volTrigRatio
if iscell(volTrigRatio)
  if length(volTrigRatio) > 1
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile) volTrigRatio should not be of length greater than one (length=%i) using only the first value of %i',length(volTrigRatio),volTrigRatio{1}));
  end
  volTrigRatio = volTrigRatio{1};
end

% get the stimfile
s = getStimfile(stimfile);
if isempty(s),return,end

% check that we have a stimfile that is interpretable
% by this program
% THIS IS DONE TO CHECK THE RETINOTOPY CODE maps onto analysis...
% [tf s taskNum] = checkStimfile(s);
% if ~tf,return,end

% check to see if a stimImage exists
% use s{1} here - [ma] 
if ~isfield(s,'pRFStimImage') 
  % somato stim image needs to be obtaine from task variables...
  % for now this is done in separate step.
  disp('your stim files need to contain pRFStimImage struct')
  keyboard
else
  % stim image was stored, just reclaim it
  disp(sprintf('(pRFGetSomatoStimImageFromStimfile) Loaded stim image from stimfile.'));
  stim = s.pRFStimImage;
end

if timeShift
  disp(sprintf('(pRFGetSomatoStimImageFromStimfile) Time shifting stimulus image by %i',timeShift));
  stim.im = circshift(stim.im,[0 0 timeShift]);
end

%%%%%%%%%%%%%%%%%%%%%
%    getStimfile    %
%%%%%%%%%%%%%%%%%%%%%
function s = getStimfile(stimfile)

s = [];

% deal with a cell array of stimfiles (like in an average)
if iscell(stimfile)
  for i = 1:length(stimfile)
    s{i} = getStimfile(stimfile{i});
    if isempty(s{i}),return;end
  end
  return
end

% load stimfile
if isstr(stimfile)
  stimfile = setext(stimfile,'mat');
  if ~mlrIsFile(stimfile)
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile) Could not open stimfile: %s',stimfile));
    return
  end
  s = load(stimfile);
elseif isstruct(stimfile)
  % see if this is a myscreen
  if isfield(stimfile,'imageWidth')
    % check for task field
    if isfield(stimfile,'task')
      s.task = stimfile.task;
      stimfile = rmfield(stimfile,'task');
    end
    % check for stimulus field
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
      stimfile = rmfield(stimfile,'stimulus');
    end
    % set myscreen field
    s.myscreen = stimfile;
  % else a variable with myscreen, task and stimulus or pRFStimImage
  elseif isfield(stimfile,'myscreen') || isfield(stimfile,'pRFStimImage')
    % copy fields over
    if isfield(stimfile,'myscreen')
      s.myscreen = stimfile.myscreen;
    end
    if isfield(stimfile,'task')
      s.task = stimfile.task;
    end
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
    end
    if isfield(stimfile,'pRFStimImage')
      s.pRFStimImage = stimfile.pRFStimImage;
    end
  end
end

% if you have a pRFStimImage then don't bother with the rest of the fields
if ~isfield(s,'pRFStimImage')
  % check fields
  checkFields = {'myscreen','task','stimulus'};
  for i = 1:length(checkFields)
    if ~isfield(s,checkFields{i})
      stimfileName = '';
      if isfield(s,'myscreen') && isfield(s.myscreen,'stimfile')
	stimfileName = getLastDir(s.myscreen.stimfile);
      end
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile) !!! Missing variable: %s in stimfile %s !!!',checkFields{i},stimfileName));
      s = [];
      return
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%
function [tf s taskNum] = checkStimfile(s)

tf = true;
s = cellArray(s);
taskNum = [];

stimulusType = [];
barAngle = [];
direction = [];

for i = 1:length(s)
  thiss = s{i};
  if isempty(thiss)
    disp(sprintf('(pRFGetsomatoStimImageFromStimfile) Missing stimfile'));
    tf = false;
    return
  end
  % if this has a pRFStimImage then we are ok
  if isfield(thiss,'pRFStimImage')
    continue;
  end
  dispstr = sprintf('%s: vols=%i',thiss.myscreen.stimfile,thiss.myscreen.volnum);
  % first check if this is a retinotpy stimfile - it should
  % have a task which is mglRetinotopy
  taskNum = [];
  for iTask = 1:2
    if (length(thiss.task) >= iTask) && (isequal(thiss.task{iTask}{1}.taskFilename,'mglRetinotopy.m') || isequal(thiss.task{iTask}{1}.taskFilename,'gruRetinotopy.m'))
      taskNum = iTask;
    end
  end
  if isempty(taskNum)
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by mglRetinotopy'));
    tf = false;
    return
  end

  % check for proper saved fields
  missing = '';
  if ~isfield(thiss.task{taskNum}{1},'randVars') missing = 'randVars';end
  if ~isfield(thiss.task{taskNum}{1},'parameter') missing = 'parameter';end
  if ~any(strcmp('maskPhase',thiss.myscreen.traceNames)) missing = 'maskPhase';end
  if ~any(strcmp('blank',thiss.myscreen.traceNames)) missing = 'blank';end
  if ~isempty(missing)
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the field %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one.',missing));
    tf = false;
    return
  end

  % check for necessary variables
  e = getTaskParameters(thiss.myscreen,thiss.task{taskNum}{1});

  % now check for each variable that we need
  varnames = {'blank'};
  for i = 1:length(varnames)
    varval = getVarFromParameters(varnames{i},e);
    if isempty(varval)
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the variable %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one',varnames{i}));
      tf = false;
      return
    end
  end
  
  % check for matching stimfiles
  if ~isempty(stimulusType) && (stimulusType ~= thiss.stimulusType)
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! Have you averaged together scans with different stimulus conditions?'));
  end
  if any(thiss.stimulus.stimulusType == [3 4])
    varval = getVarFromParameters('barAngle',e);
    if ~isempty(barAngle) && ~isequal(varval,barAngle)
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The barAngles are different! Have you averaged together scans with different stimulus conditions?'));
    end
    barAngle = varval;
  else
    if ~isempty(direction) && (thiss.stimulus.direction ~= direction)
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The directions are different! Have you averaged together scans with different stimulus conditions?'));
    end
    direction = thiss.stimulus.direction;
  end
end

s = s{end};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    saveStimImageToStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveStimImageToStimfile(stim,stimfile)

% make sure stimfile is a cell array
stimfile = cellArray(stimfile);

% first reload the stimfile
for iStimfile = 1:length(stimfile)
  if isfield(stimfile{iStimfile},'filename')
    s = load(stimfile{iStimfile}.filename);
    if isempty(s)
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:saveStimImageToStimfile) Could not load stimfile %s. Unable to save stim image back to stimfile',stimfile{iStimfile}.filename));
    else
      % append the stim image and save back
      s.pRFStimImage = stim;
      save(stimfile{iStimfile}.filename,'-struct','s');
      disp(sprintf('(pRFGetSomatoStimImageFromStimfile:saveStimImageToStimfile) Saved stimImage to %s.',stimfile{iStimfile}.filename));
    end
  else
    disp(sprintf('(pRFGetSomatoStimImageFromStimfile:saveStimImageToStimfile) Missing filename in stimfile structure, could not save stimImage back to stimfile'));
  end
end

