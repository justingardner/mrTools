function transformROIs(thisView)
% transformROIs(thisView)
%
%   transforms  ROI(s) using pre-definde or custom functions
%
% jb 11/01/2011
%
% $Id: transformROIs.m 1982 2010-12-20 21:12:20Z julien $

%default params
%get names of combine Functions in transformFunctions directory
functionsDirectory = [fileparts(which('transformROIs')) '/TransformROIsFunctions/'];
transformFunctionFiles =  dir([functionsDirectory '*.m']);
for iFile=1:length(transformFunctionFiles)
   transformFunctions{iFile} = stripext(transformFunctionFiles(iFile).name);
end

params.transformFunction = [{'User Defined'} transformFunctions];
params.customTransformFunction = '';

currentBaseString = ['Current Base (' viewGet(thisView,'basename') ')'];
roiSpaceMenu = {'Native','Current scan',currentBaseString};

askForParams = 1;
while askForParams
  params = {...
   {'transformFunction',params.transformFunction,'type=popupmenu','name of the function to apply. This is a list of existing functions in the transformROIsFunctions directory. To get help for a specific function, type ''help functionName''. To use another function, select ''User Defined'' and type the function name below'},...
   {'customTransformFunction',params.customTransformFunction,'name of the function to apply. You can use any custom matlab function on the path that accepts an ROI structure as argument and output a new ROI structure.'},...
   {'roiSpace',roiSpaceMenu,'In which space should the coordinates be converted before being passed'},...
   {'additionalArgs','','Additional arguments to the transform function. These arguments will be input at the end of each function call. They must be separated by commas. '},...
          };
  params = mrParamsDialog(params, 'Choose an ROI transfromation function');
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
    roiList = selectInList(thisView,'rois');
    if isempty(roiList)
       askForParams = 1;
    end
  end
end

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
for roiNum=roiList
  roi = viewGet(thisView,'roi',roiNum);
  if ~strcmp(params.roiSpace,'Native')
    disp(sprintf('(convertROI) Converting ROI %i:%s to %s coordinate space',roiNum,roi.name,params.roiSpace));
    if ~isempty(roi.coords)
      roi.coords = getROICoordinates(thisView,roiNum,whichVolume,[],baseNum);
    else
      mrWarnDlg(sprintf('(convertROI) ROI %i:%s has empty coordinates in transformation, skipping conversion...',roiNum,roi.name));
    end
    roi.sformCode = newSformCode;
    roi.xform = newXform;
    roi.vol2mag = newVol2mag;
    roi.vol2tal = newVol2tal;
    roi.voxelSize = newVoxelSize;
  end
  
  %parse other additional inputs
  additionalArgs = parseArguments(params.additionalArgs,',');
  functionString = [' ' params.transformFunction '(roi'];
  for iArg = length(additionalArgs)
    functionString = [functionString ',' additionalArgs{iArg}];
  end
  functionString(end+(1:2)) = ');';
  roi = eval(functionString);
  if ~isempty(roi)
    thisView = viewSet(thisView,'newROI',roi);
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
   if isnumeric(str2num(token))
      arguments{nArgs} = token;
   elseif isempty(arguments{nArgs})
      arguments{nArgs} = ['''' token ''''];
   end
end        
