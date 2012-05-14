function combineTransformOverlays(thisView)
% combineTransformOverlays(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%   combines (masked) Overlays according to matlab or custom operators  in current view and current analysis
%
% jb 07/07/2010
%
% $Id$

inputOutputTypeMenu = {'3D Array','Scalar','Structure'};%,'Anonymous (Scalar)', 'Anonymous (Array)'};
combinationModeMenu = {'Apply function to all overlays','Apply function to each overlay','Recursively apply to overlay pairs'};

%default params
%get names of combine Functions in combineFunctions directory
functionsDirectory = [fileparts(which('combineTransformOverlays')) '/combineTransformOverlayFunctions/'];
combineFunctionFiles =  dir([functionsDirectory '*.m']);
for iFile=1:length(combineFunctionFiles)
  combineFunctions{iFile} = stripext(combineFunctionFiles(iFile).name);
end

combineFunctionsMenu = [{'User Defined'} combineFunctions];
params.customCombineFunction = '';%(''@(x)max(0,-norminv(x))';
params.nOutputOverlays = 1;
params.additionalArrayArgs = '';
params.additionalArgs = '';
params.useMask = 0;
params.passView = 0;
params.outputName = '';

askForParams = 1;
while askForParams
  params = {...
   {'combineFunction',combineFunctionsMenu,'type=popupmenu','name of the function to apply. This is a list of existing combine functions in the combineFunctions directory. To use another function, select ''User Defined'' and type the function name below'},...
   {'customCombineFunction',params.customCombineFunction,'name of the function to apply. You can use any type of matlab function (including custom) that accepts either scalars or multidimensional arrays. Any string beginning with an @ will be considered an anonymous function and shoulde be of the form @(x)func(x), @(x,y)func(x,y) ..., where the number of variables equals the number of overlay inputs and additional arguments. '},...
   {'inputOutputType',inputOutputTypeMenu,'type=popupmenu','Type of arguments accepted by the combination function. ''3D Array'' is faster but not all functions accept multidimensional arrays as inputs. Choose ''Structure'' to pass the whole overlay structure'},...
   {'combinationMode',combinationModeMenu,'type=popupmenu', 'In the default mode,the number of inputs expected by the function must match the number of selected overlays.'},...
   {'nOutputOverlays',params.nOutputOverlays,'incdec=[-1 1]','round=1','minmax=[1 Inf]','Number of outputs of the combineFunctions'},...
   {'additionalArrayArgs',params.additionalArrayArgs,'constant arguments for functions that accept them. Arguments must be separated by commas. for Array input/output type, each argument will be repeated in a matrix of same dimensions of the overlay '},...
   {'additionalArgs',params.additionalArgs,'constant scalar arguments for functions that take both arrays and scalars. These arguments will be input at the end of each function call. They must be separated by commas. '},...
   {'passView',params.passView,'type=checkbox','Check this if the function requires the current mrLoadRet view'},...
   {'useMask',params.useMask,'type=checkbox','use overlay masked by other overlays of the analysis according to clip values'},...
   {'outputName',params.outputName,'radical of the output overlay names'},...
   {'printHelp',0,'type=pushbutton','callback',@printHelp,'passParams=1','buttonString=Print combineFunction Help','Prints combination function help in command window'},...
          };
  params = mrParamsDialog(params, 'Overlay Combination parameters');
  % Abort if params empty
  if ieNotDefined('params'),return,end

  if strcmp(params.combineFunction,'User Defined')
    params.combineFunction = params.customCombineFunction;
  end

  inputOutputTypeMenu = putOnTopOfList(params.inputOutputType,inputOutputTypeMenu);
  combinationModeMenu = putOnTopOfList(params.combinationMode,combinationModeMenu);
  combineFunctionsMenu = putOnTopOfList(params.combineFunction,combineFunctionsMenu);

  if strcmp(params.combinationMode,'Recursively apply to overlay pairs') && params.combineFunction(1)=='@'
    mrWarnDlg('(combineTransformOverlays) Anonymous functions cannot be applied recursively');
  elseif isempty(params.combineFunction) && isempty(params.customCombineFunction)
    mrWarnDlg('(combineTransformOverlays) Please choose a combination/transformation function');
  %elseif
    %other controls here
  else
    askForParams = 0;
    overlayList = selectInList(thisView,'overlays');
    if isempty(overlayList)
       askForParams = 1;
    end
  end
end
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;

%get and mask the overlay data
nScans = viewGet(thisView,'nScans');   
overlayData = viewGet(thisView,'overlays');
overlayData = overlayData(overlayList);
if params.useMask
   mask = maskOverlay(thisView,overlayList);
   for iScan = 1:length(mask)
      for iOverlay = 1:length(overlayData)
         overlayData(iOverlay).data{iScan}(~mask{iScan}(:,:,:,iOverlay))=NaN;
      end
   end
end

%overlay names
for iOverlay = 1:length(overlayList)
  overlayNames{iOverlay} = overlayData(iOverlay).name;
end

%reformat input data
switch(params.inputOutputType)
  case 'Structure'
    overlayData = num2cell(overlayData);
  case {'3D Array','Scalar'}
    newOverlayData = cell(nScans,length(overlayList));
    for iOverlay = 1:length(overlayList)
      for iScan = 1:nScans
        newOverlayData{iScan,iOverlay} = overlayData(iOverlay).data{iScan};
      end
    end
    overlayData = newOverlayData;
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
   for iScan = 1:nScans
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

if params.passView
  functionString = [functionString 'thisView,'];
end

functionString(end:end+1) = ');'; %replace the last comma by a closing bracket to end the function



%evaluate the function 
outputData = cell(size(overlayData,1),params.nOutputOverlays, size(overlayData,3));
outputOverlayNames = cell(params.nOutputOverlays, size(overlayData,3));
for iOperations = 1:size(overlayData,3)
  for iScan = 1:size(overlayData,1)
    %check for empty overlays (only if 3D array or scalar)
    emptyInput=false;
    if ismember(params.inputOutputType,{'3D Array','Scalar'})
      for iInput = 1:size(overlayData,2)
        if isempty(overlayData{iScan,iInput,iOperations})
          emptyInput=true;
        end
      end
    end
    if emptyInput
      for iOutput = 1:params.nOutputOverlays
        outputData{iScan,iOutput,iOperations} = [];
      end
    else
      try
%         tic
         eval(functionString);
%          toc
      catch exception
         mrWarnDlg(sprintf('There was an error evaluating function %s:\n%s',combineFunctionString,getReport(exception,'basic')));
         return
      end
      for iOutput = 1:params.nOutputOverlays
    %          if  strcmp(params.inputOutputType,'Scalar')
    %             %convert back to numerical value
    %             outputData{iScan,iOverlay,iOperations} = cell2mat(outputData{iScan,iOverlay,iOperations}); 
    %          end
         %add 0 to all the results to convert logical to doubles, because mrLoadRet doesn't like logical overlays
        switch(params.inputOutputType)
          case 'Structure'
            for jScan = 1:nScans
              outputData{iScan,iOutput,iOperations}.data{jScan} = outputData{iScan,iOutput,iOperations}.data{jScan}+0;
            end
          case {'3D Array','Scalar'}
           outputData{iScan,iOutput,iOperations} = outputData{iScan,iOutput,iOperations}+0;
        end
        %check that the size is compatible
        if ~isequal(size(outputData{iScan,iOutput,iOperations}),size(overlayData{iScan}))
          mrWarnDlg(['(combineTransformOverlays) Dimensions of result are not compatible with overlay ([' num2str(size(outputData{iScan,iOutput,iOperations})) '] vs [' num2str(size(overlayData{iScan})) '])']);
        end
      end
    end
  end
end

%name of output overlays
for iOutput=1:params.nOutputOverlays
   if params.nOutputOverlays>1
      name = ['Ouput ' num2str(iOutput) ' - '];
   else
      name = '';
   end
   if ~isempty(params.outputName)
      name = [name params.outputName '('];
   else
      name = [name params.combineFunction '('];
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
            
%save the output
defaultOverlay.groupName = viewGet(thisView,'groupName');
defaultOverlay.function = 'combineTransformOverlays';
defaultOverlay.interrogator = 'combineTransformOverlays';
defaultOverlay.type = 'combination';
defaultOverlay.params = params;
defaultOverlay.clip = [];
defaultOverlay.range = [];
defaultOverlay.name = [];
defaultOverlay.data = [];
for iOverlay = 1:size(outputData,2)
  switch(params.inputOutputType)
    case {'3D Array','Scalar'}
      outputOverlay(iOverlay) = defaultOverlay;
      outputOverlay(iOverlay).data = outputData(:,iOverlay);
      for iOutput = 1:size(outputData,1)
        isNotEmpty(iOutput) = ~isempty(outputData{iOutput,iOverlay});
      end
      allScansData = cell2mat(outputData(isNotEmpty,iOverlay));
    case 'Structure'
      outputOverlay(iOverlay) = copyFields(defaultOverlay,outputData{iOverlay});
      allScansData = cell2mat(outputOverlay(iOverlay).data);
  end
  maxValue = max(allScansData(allScansData<inf));
  minValue = min(allScansData(allScansData>-inf));
  outputOverlay(iOverlay).clip = [minValue maxValue];
  outputOverlay(iOverlay).range = [minValue maxValue];
  outputOverlay(iOverlay).name = outputOverlayNames{iOverlay};
end
thisView = viewSet(thisView,'newoverlay',outputOverlay);

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

function printHelp(params)

if strcmp(params.combineFunction,'User Defined')
  mrWarnDlg('(combineTransformOverlays) Please select a combination function');
else
  helpString = help(params.combineFunction);
  if isempty(helpString)
    mrWarnDlg('(combineTransformOverlays) No help available');
  else
    fprintf('\n');
    disp(helpString);
  end
end

