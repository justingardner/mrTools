function saveAnalysis(view,analysisName,confirm)
%
%   saveAnalysis(view,[analysisName],[confirm])
%
% Saves an analysis to a file. The filename = analysis.name.
%
% analysisname: Can be either the name or the number.
%          Default: current analysis
% confirm: If filename already exists, prompt user to over-write. 
%          Default: uses mrLoadRet 'verbose' preference or 0 (if preference
%          not defined.
%
% djh, 7/2006

if ieNotDefined('analysisName')
    analysisNum = viewGet(view,'currentAnalysis');
    analysisName = viewGet(view,'analysisName',analysisNum);
end
if isstr(analysisName)
    analysisNum = viewGet(view,'analysisNum',analysisName);
elseif isnumeric(analysisName)
	analysisNum = analysisName;
	analysisName = viewGet(view,'analysisName',analysisNum);
else
	myErrorDlg(['Bad analysis name: ',analysisName]);
end

if ieNotDefined('confirm')
    if ispref('mrLoadRet','verbose')
        confirm = getpref('mrLoadRet','verbose');
    else
        confirm = 0;
    end
end

% Assign local variable with analysisName = analysis
analysis = viewGet(view,'analysis',analysisNum);
eval([analysisName,'=analysis;']);

% Path
analysisdir = viewGet(view,'analysisdir',[],analysisNum);
filename = analysisName;
pathStr = fullfile(analysisdir,filename);
        
% Write, though check for over-writing
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
    if confirm
        saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
            'Save Analysis?','Yes','No','No');
	end
end
if strcmp(saveFlag,'Yes')
    fprintf('Saving %s...',pathStr);
    saveString = ['save(pathStr,','''',analysisName,'''',');'];
    eval(saveString);
    fprintf('done\n');
else
    fprintf('Analysis not saved...');
end
return;
