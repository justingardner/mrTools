function saveROI(view,roiName,confirm)
%
% saveROI(view,[roiName],[confirm])
%
% Saves an ROI to a file. The filename = ROI.name.
%
% roiName: Either the ROI name or number.
%          Default: current overlay.
% confirm: If filename already exists, prompt user to over-write.
%          Default: uses mrLoadRet 'verbose' preference or 0 (if preference
%          not defined.
%
% djh, 9/2005

if ieNotDefined('roiName')
	roiNum = viewGet(view,'currentROI');
	roiName = viewGet(view,'roiName',roiNum);
end
if isstr(roiName)
    roiNum = viewGet(view,'roiNum',roiName);
elseif isnumeric(roiName)
	roiNum = roiName;
	roiName = viewGet(view,'roiName',roiNum);
else
	myErrorDlg(['Bad ROI name: ',roiName]);
end

if ieNotDefined('confirm')
    if ispref('mrLoadRet','verbose')
        confirm = getpref('mrLoadRet','verbose');
    else
        confirm = 0;
    end
end

% Assign local variable with roiName = roi
roi = viewGet(view,'roi',roiNum);
% now go through and fix any characters that are not allowed in
% variable names
swapchars = {{'-','_'},{' ','_'},{'*','s'},{'+','p'},{'%','p'},{'[','B'},{']','B'},{'(','P'},{')','P'}};
for i = 1:length(swapchars)
  swaplocs = strfind(roiName,swapchars{i}{1});
  roiName(swaplocs) = swapchars{i}{2};
end
roiName
eval([roiName,'=roi;']);

% path to file
filename = roiName;
pathStr = fullfile(viewGet(view,'roiDir'),filename);

% Write, though check for over-writing
%
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
	if confirm
		saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
			'Save ROI?','Yes','No','No');
	end
end
if strcmp(saveFlag,'Yes')
	fprintf('Saving %s...',pathStr);
	saveString = ['save(pathStr,','''',roiName,'''',');'];
	eval(saveString);
	fprintf('done\n');
else
	fprintf('ROI not saved...');
end
return;

% Test/debug
saveROI(MLR.views{1},'ROI1');
saveROI(MLR.views{1});
