function combineOverlays(~,~,thisView)
% combineOverlays(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%   combines (masked) Overlays according to matlab or custom operators  in current view and current analysis
%
% jb 07/07/2010
%
% $Id$

inputOutputTypeMenu = {'3D Array','Scalar'};%,'Anonymous (Scalar)', 'Anonymous (Array)'};
combinationModeMenu = {'Apply function to all overlays','Apply function to each overlay','Recursively apply to overlay pairs'};

%default params
%get names of combine Functions in combineFunctions directory
functionsDirectory = [fileparts(which('combineOverlays')) '/CombineOverlayFunctions/'];
combineFunctionFiles =  dir([functionsDirectory '*.m']);
for iFile=1:length(combineFunctionFiles)
   combineFunctions{iFile} = stripext(combineFunctionFiles(iFile).name);
end

params.combineFunction = [{'User Defined'} combineFunctions];
params.customCombineFunction = '';%(''@(x)max(0,-norminv(x))';
params.nOutputOverlays = 1;
params.additionalArrayArgs = '';
params.additionalArgs = '';
params.useMask = 0;
params.outputName = '';

askForParams = 1;
while askForParams
  params = {...
   {'combineFunction',params.combineFunction,'type=popupmenu','name of the function to apply. This is a list of existing combine functions in the combineFunctions directory. To use another function, select ''User Defined'' and type the function name below'},...
   {'customCombineFunction',params.customCombineFunction,'name of the function to apply. You can use any type of matlab function (including custom) that accepts either scalars or multidimensional arrays. Any string beginning with an @ will be considered an anonymous function and shoulde be of the form @(x)func(x), @(x,y)func(x,y) ..., where the number of variables equals the number of overlay inputs and additional arguments. '},...
   {'inputOutputType',inputOutputTypeMenu,'type=popupmenu','Type of arguments accepted by the combination function. ''3D Array'' is faster but not all functions accept multidimensional arrays as inputs.'},...
   {'combinationMode',combinationModeMenu,'type=popupmenu', 'In the default mode,the number of inputs expected by the function must match the number of selected overlays.'},...
   {'nOutputOverlays',params.nOutputOverlays,'incdec==[-1 1]','round=1','minmax=[1 Inf]','Number of outputs of the combineFunctions'},...
   {'additionalArrayArgs',params.additionalArrayArgs,'constant arguments for functions that accept them. Arguments must be separated by commas. for Array input/output type, each argument will be repeated in a matrix of same dimensions of the overlay '},...
   {'additionalArgs',params.additionalArgs,'constant scalar arguments for functions that take both arrays and scalars. These arguments will be input at the end of each function call. They must be separated by commas. '},...
   {'useMask',params.useMask,'type=checkbox','use overlay masked by other overlays of the analysis according to clip values'}...
   {'outputName',params.outputName,'radical of the output overlay names'}...
          };
  params = mrParamsDialog(params, 'Overlay Combination parameters');
  % Abort if params empty
  if ieNotDefined('params'),return,end

  if strcmp(params.combineFunction,'User Defined')
    params.combineFunction = params.customCombineFunction;
  end

  inputOutputTypeMenu = putOnTopOfList(params.inputOutputType,inputOutputTypeMenu);
  combinationModeMenu = putOnTopOfList(params.combinationMode,combinationModeMenu);

  if strcmp(params.combinationMode,'Recursively apply to overlay pairs') && params.combineFunction(1)=='@'
    mrWarnDlg('Anonymous functions cannot be applied recursively');
  %elseif
    %other control here
  else
    askForParams = 0;
    overlayList = selectOverlays(thisView);
    if isempty(overlayList)
       askForParams = 1;
    end
  end
end
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;



%get and mask the overlay data
if params.useMask
   [mask,overlayData] = maskOverlay(thisView,overlayList);
   if ~iscell(mask)
      mask = num2cell(mask,1,1);
      overlayData = num2cell(overlayData,1,1);
   end
   for iScan = 1:length(mask)
      for iOverlay = 1:size(overlayData,2)
         overlayData{iScan,iOverlay}(~mask{iScan})=0;
      end
   end
else
   scanList = 1:viewGet(thisView,'nScans');   
   overlayData = cell(length(scanList),length(overlayList));
   for iScan = 1:length(scanList)
      for iOverlay = 1:length(overlayList)
         overlayData{iScan,iOverlay} = viewGet(thisView,'overlayData',scanList(iScan),overlayList(iOverlay));
      end
   end
end
%overlay names
for iInput = 1:size(overlayData,2)
    overlayNames{iInput} = viewGet(thisView,'overlayName',overlayList(iInput));
end

%parse additional array inputs
additionalArrayArgs = parseArguments(params.additionalArrayArgs,',');
if ~isempty(additionalArrayArgs)
   if all(cellfun(@isnumeric,additionalArrayArgs)) && strcmp(params.inputOutputType,'3D Array') %if all arguments are numeric and the input type is Array
      additionalArrayInputs = cellfun(@(x)repmat(x,[size(overlayData{1}) 1]),additionalArrayArgs,'UniformOutput',false); %convert additional arguments to arrays
   else %if any additional argument is not a number
      additionalArrayInputs = cellfun(@(x)num2cell(repmat(x,[size(overlayData{1}) 1])),additionalArrayArgs,'UniformOutput',false); %convert additional arguments to cell arrays
      params.inputOutputType = 'Scalar';  %and force scalar
   end
else
   additionalArrayInputs = {};
end

%parse other additional inputs
additionalArgs = parseArguments(params.additionalArgs,',');

%convert overlays to cell arrays if scalar function
if strcmp(params.inputOutputType,'Scalar')
   for iScan = 1:length(scanList)
      for iOverlay = 1:length(overlayList)
         overlayData{iScan,iOverlay} = num2cell(overlayData{iScan,iOverlay}); %convert overlays to cell arrays
      end
   end
end

%reshape the input if function run separately on each input
if strcmp(params.combinationMode,'Apply function to each overlay')
   overlayData = cellReshape(overlayData,size(overlayData,1),size(overlayData,3),size(overlayData,2));
end


if strcmp(params.combinationMode,'Recursively apply to overlay pairs') %should add additional non-array arguments also ?
   nTotalargs = size(overlayData,2)+length(additionalArrayArgs);
   combineFunctionString = '@(x1';
   for iInput = 2:nTotalargs
      combineFunctionString = [combineFunctionString ',x' num2str(iInput)];
   end
   combineFunctionString = [combineFunctionString ')' ];
   for iInput = nTotalargs:-1:2
      combineFunctionString = [combineFunctionString params.combineFunction '(x' num2str(iInput) ','];
   end
   combineFunctionString = [combineFunctionString 'x1'];
   for iInput = 1:nTotalargs-1
      combineFunctionString = [combineFunctionString ')'];
   end
else
   combineFunctionString = params.combineFunction;
end
combineFunctionHandle = str2func(combineFunctionString);


%-----------------------------construct function call
%output arguments
functionString = '[';
for iOutput = 1:params.nOutputOverlays
   functionString = [functionString 'outputData{iScan,' num2str(iOutput) ',iOperations},'];
end
functionString(end)=']';

%array operator
if strcmp(params.inputOutputType,'Scalar') %if function operates on scalar
   functionString = [functionString ' = cellfun(combineFunctionHandle,']; %we'll use cellfun to apply the function to each voxel
else                                       %if it operates on arrays
   functionString = [functionString ' = feval(combineFunctionHandle,'];
end

%input arguments
for iInput = 1:size(overlayData,2)
   functionString = [functionString 'overlayData{iScan,' num2str(iInput) ',iOperations},'];
end
   
%additional arguments
for iInput = 1:length(additionalArrayInputs)
   functionString = [functionString 'additionalArrayInputs{' num2str(iInput) '},'];
end

%additional scalar arguments
for iInput = 1:length(additionalArgs)
   functionString = [functionString 'additionalArgs{' num2str(iInput) '},'];
end

functionString(end:end+1) = ');'; %replace the last comma by a closing bracket to end the function



%evaluate the function 
outputData = cell(size(overlayData,1),params.nOutputOverlays, size(overlayData,3));
outputOverlayNames = cell(params.nOutputOverlays, size(overlayData,3));
for iOperations = 1:size(overlayData,3)
   for iScan = 1:size(overlayData,1)
      try
         eval(functionString);
      catch exception
         mrWarnDlg(sprintf('There was an error evaluating the combining function:\n%s',getReport(exception,'basic')));
         return
      end
      for iOuput = 1:params.nOutputOverlays
%          if  strcmp(params.inputOutputType,'Scalar')
%             %convert back to numerical value
%             outputData{iScan,iOverlay,iOperations} = cell2mat(outputData{iScan,iOverlay,iOperations}); 
%          end
         %add 0 to all the results to convert logical to doubles, because mrLoadRet doesn't like logical overlays
         outputData{iScan,iOuput,iOperations} = outputData{iScan,iOuput,iOperations}+0;
         %check that the size is compatible
         if any(size(outputData{iScan,iOuput,iOperations})~=size(overlayData{1}))
            mrWarnDlg(['(combineOverlays) Dimensions of result are not compatible with overlay ([' num2str(size(outputData{iScan,iOuput,iOperations})) '] vs [' num2str(size(overlayData{1})) '])']);
         end
      end
   end
end

%name of output overlays
for iOutput=1:params.nOutputOverlays
   if params.nOutputOverlays>1
      name = ['Ouput ' num2str(iOuput) ' - '];
   else
      name = '';
   end
   if ~isempty(params.outputName)
      name = [name params.outputName '('];
   else
      name = [params.combineFunction '('];
   end
   for iOperations = 1:size(overlayData,3)
      if size(overlayData,3)>1
         outputOverlayNames{iOutput,iOperations} = [name overlayNames{iOperations} ','];
      else
         outputOverlayNames{iOutput,iOperations} = name;
         for iInput = 1:size(overlayData,2)
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} overlayNames{iInput} ','];
         end
      end
      for iInput = 1:length(additionalArrayArgs)
         if isnumeric(additionalArrayArgs{iInput})
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} num2str(additionalArrayArgs{iInput}) ','];
         else
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} additionalArrayArgs{iInput} ','];
         end
      end
      for iInput = 1:length(additionalArgs)
         if isnumeric(additionalArgs{iInput})
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} num2str(additionalArgs{iInput}) ','];
         elseif isa(additionalArgs{iInput},'function_handle')
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} func2str(additionalArgs{iInput}) ','];
         else
            outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} additionalArgs{iInput} ','];
         end
      end
      outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations}(1:end-1) ')'];
   end
end

if size(overlayData,3)>1 %reshape the overlay cell array
   outputData = cellReshape(outputData,size(outputData,1),size(outputData,2)*size(outputData,3));
   outputOverlayNames = cellReshape(outputOverlayNames,numel(outputOverlayNames),1);
end
            
%saves the output
outputOverlay.groupName = viewGet(thisView,'groupName');
outputOverlay.function = 'combineOverlays';
outputOverlay.interrogator = 'combineOverlays';
outputOverlay.type = 'combination';
outputOverlay.params = params;
for iOverlay = 1:size(outputData,2)
%   outputOverlay = overlay;
   outputOverlay.data = outputData(:,iOverlay);
   allScansData = cell2mat(outputData);
   maxValue = max(allScansData(allScansData<inf));
   minValue = min(allScansData(allScansData>-inf));
   outputOverlay.clip = [minValue maxValue];
   outputOverlay.range = [minValue maxValue];
   outputOverlay.name = outputOverlayNames{iOverlay};
   thisView = viewSet(thisView,'newoverlay',outputOverlay);
end
   
refreshMLRDisplay(thisView.viewNum);
set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow;

function [arguments, nArgs] = parseArguments(argumentString, separator)
   
%parse string of arguments separated by separator and put them into a cell array of numerical and string arguments
%non-numerical values that are not between quotes are converted into strings
%
% Julien Besle, 08/07/2010
nArgs = 0;
arguments = cell(0);
remain = argumentString;
while ~isempty(remain)
   nArgs = nArgs+1;
   [token,remain] = strtok(remain, separator);
   try
      arguments{nArgs} = eval(token);
   catch exception
      if strcmp(exception.identifier,'MATLAB:UndefinedFunction')
         arguments{nArgs} = token;
      else
         mrErrorDlg(['(parseArguments) could not read argument: ' exception.message]);
      end
   end
      
end
