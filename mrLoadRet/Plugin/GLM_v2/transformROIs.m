function [thisView,params] = transformROIs(thisView,params,varargin)
% transformROIs(thisView)
%
%   transforms  ROI(s) using pre-defined or custom functions
%
% jb 11/01/2011
%
%             To just get a default parameter structure:
% 
%             v = newView;
%             [v params] = transformROIs(v,[],'justGetParams=1');
%             [v params] = transformROIs(v,[],'justGetParams=1','defaultParams=1');
%             [v params] = transformROIs(v,[],'justGetParams=1','defaultParams=1','roiList=[1 2]');
%
% $Id: transformROIs.m 1982 2010-12-20 21:12:20Z julien $

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

currentBaseString = ['Current Base (' viewGet(thisView,'basename') ')'];

if ieNotDefined('params')
  %default params
  %get names of combine Functions in transformFunctions directory
  functionsDirectory = [fileparts(which('transformROIs')) '/transformROIFunctions/'];
  transformFunctionFiles =  dir([functionsDirectory '*.m']);
  for iFile=1:length(transformFunctionFiles)
     transformFunctions{iFile} = stripext(transformFunctionFiles(iFile).name);
  end

  %get names of transform functions in additional folder(s)
  roiTransformPaths = commaDelimitedToCell(mrGetPref('roiTransformPaths'));
  for i = 1:length(roiTransformPaths)
    roiTransformFiles =  dir([roiTransformPaths{i} '/*.m']);
    for iFile=1:length(roiTransformFiles)
       transformFunctions{end+1} = stripext(roiTransformFiles(iFile).name);
    end
  end
  transformFunctions = sort(transformFunctions); %re-order in alphabetical order

  params.transformFunction = [{'User Defined'} transformFunctions];
  params.customTransformFunction = '';
  params.roiNameSuffix = '';
  params.newRoiName = '';

  roiSpaceMenu = {'Native','Current scan',currentBaseString};
  passRoiModeMenu = {'One ROI at a time','All ROIs at once'};

  askForParams = 1;
  while askForParams
    params = {...
     {'transformFunction',params.transformFunction,'type=popupmenu','name of the function to apply. This is a list of existing functions in the transformROIFunctions directory. To get help for a specific function, type ''help functionName''. To use another function, select ''User Defined'' and type the function name below'},...
     {'customTransformFunction',params.customTransformFunction,'name of the function to apply. You can use any custom matlab function on the path that accepts an ROI structure as argument and output a new ROI structure.'},...
     {'passRoiMode',passRoiModeMenu,'Specifies if the transform function accepts one or several ROIs as inputs'},...
     {'newRoiName',params.newRoiName,'transformed ROI name (leave blank if ROI name should stay the same).'},...
     {'roiNameSuffix',params.roiNameSuffix,'suffix that will be appended to the transformed ROI name (leave blank if not suffix should be appended).'},...
     {'roiSpace',roiSpaceMenu,'In which space should the coordinates be converted before being passed'},...
     {'additionalArgs','','Additional arguments to the transform function. These arguments will be input at the end of each function call. They must be separated by commas.'},...
     {'printHelp',0,'type=pushbutton','callback',@printHelp,'passParams=1','buttonString=Print transformFunction Help','Prints transformation function help in command window'},...
            };

    % Initialize analysis parameters with default values
    if defaultParams
      params = mrParamsDefault(params);
    else
      params = mrParamsDialog(params, 'ROI transformation parameters');
    end
    % Abort if params empty
    if ieNotDefined('params'),return,end

    if strcmp(params.transformFunction,'User Defined')
      params.transformFunction = params.customTransformFunction;
    end

    if 0
      %control here
    %elseif
      %other control here
    else
      askForParams = 0;
      if defaultParams
        params.roiList = viewGet(thisView,'curROI');
      else
        params.roiList = selectInList(thisView,'rois');
        if isempty(params.roiList)
           askForParams = 1;
        end
      end
    end
  end
end

if ~ieNotDefined('roiList')
  params.roiList = roiList;
end

% if just getting params then return
if justGetParams,return,end

switch(params.roiSpace)
  case currentBaseString
      baseNum = viewGet(thisView,'currentbase');
      newXform = viewGet(thisView,'baseXform',baseNum);
      newSformCode = viewGet(thisView,'baseSformCode',baseNum);
      newVol2mag = viewGet(thisView,'baseVol2mag',baseNum);
      newVol2tal = viewGet(thisView,'baseVol2tal',baseNum);
      newVoxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
      whichVolume = 0;
    
  case 'Current scan'
      newXform = viewGet(thisView,'scanXform');
      newSformCode = viewGet(thisView,'scanSformCode');
      newVol2mag = viewGet(thisView,'scanVol2mag');
      newVol2tal = viewGet(thisView,'scanVol2tal');
      newVoxelSize = viewGet(thisView,'scanVoxelSize');
      whichVolume = [];
      baseNum = [];
end

needToRefresh = 0;
cRoi = 0;
for iRoi=params.roiList
  cRoi = cRoi+1;
  roi = viewGet(thisView,'roi',iRoi);
  if ~strcmp(params.roiSpace,'Native')
    disp(sprintf('(convertROI) Converting ROI %i:%s to %s coordinate space',iRoi,roi.name,params.roiSpace));
    if ~isempty(roi.coords)
      roi.coords = getROICoordinates(thisView,iRoi,whichVolume,[],'baseNum',baseNum);
    else
      mrWarnDlg(sprintf('(convertROI) ROI %i:%s has empty coordinates in transformation, skipping conversion...',iRoi,roi.name));
    end
    roi.sformCode = newSformCode;
    roi.xform = newXform;
    roi.vol2mag = newVol2mag;
    roi.vol2tal = newVol2tal;
    roi.voxelSize = newVoxelSize;
  end
  switch(params.passRoiMode)
    case 'One ROI at a time'
      rois{cRoi} = roi;
      
    case 'All ROIs at once'
      rois{1}(cRoi) = roi;
  end
end

%parse other additional inputs
additionalArgs = parseArguments(params.additionalArgs,',');
%construct function call
functionString='';
for iArg = 1:length(additionalArgs)
  functionString = [functionString ',' additionalArgs{iArg}];
end
functionString(end+(1:2)) = ');';
  
for iCall = 1:length(rois)
  if strcmp(params.passRoiMode,'One ROI at a time')
    disp(sprintf('(transformROI) Calling %s for ROI %s', params.transformFunction, rois{iCall}.name));
  end
  try
    roi = eval([params.transformFunction '(rois{iCall}' functionString]);
  catch exception
     mrWarnDlg(sprintf('(transformROI) There was an error evaluating function %s:\n%s\n',functionString,getReport(exception)));
     return
  end
  for iRoi = 1:length(roi)
    if ~fieldIsNotDefined(params,'newRoiName')
      roi(iRoi).name = params.newRoiName;
      if length(roi)>1
        roi(iRoi).name = [roi(iRoi).name '_' num2str(iRoi)];
      end
    end
    if ~fieldIsNotDefined(params,'roiNameSuffix')
      roi(iRoi).name = [roi(iRoi).name params.roiNameSuffix];
    end
    thisView = viewSet(thisView,'newROI',roi(iRoi));
    needToRefresh = 1;
  end
end
if needToRefresh
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end

function [arguments, nArgs] = parseArguments(argumentString, separator)
   
%parse string of arguments separated by separator and put them into a cell array of 
%  - strings if numerical 
%  - strings with double quotes for non-numerical values
% so that it can be used with eval
% Julien Besle, 08/07/2010
nArgs = 0;
arguments = cell(0);
remain = argumentString;
while ~isempty(remain)
   nArgs = nArgs+1;
   [token,remain] = strtok(remain, separator);
   if ~isempty(str2num(token))
      arguments{nArgs} = token;
   else
      arguments{nArgs} = ['''' token ''''];
   end
end   

function printHelp(params)

if strcmp(params.transformFunction,'User Defined')
  mrWarnDlg('(transformROIs) Please select a combination function');
else
  helpString = help(params.transformFunction);
  if isempty(helpString)
    mrWarnDlg('(transformROIs) No help available');
  else
    disp(helpString)
  end
end


