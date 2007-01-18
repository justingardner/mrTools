function saveOverlay(view,overlayName,analysisName,confirm)
%
%   saveOverlay(view,[overlayName],[analysisName],[confirm])
%
% Saves an overlay to a file. The filename = overlay.name.
%
% overlayName: Can be either the name or the number. 
%          Default: current overlay.
% analysisname: Can be either the name or the number.
%          Default: current analysis
% confirm: If filename already exists, prompt user to over-write. 
%          Default: uses mrLoadRet 'verbose' preference or 0 (if preference
%          not defined.
%
% djh, 7/2004
% djh, 7/2006 updated with analysisName

if ieNotDefined('overlayName')
    overlayNum = viewGet(view,'currentOverlay');
    overlayName = viewGet(view,'overlayName',overlayNum);
end
if isstr(overlayName)
    overlayNum = viewGet(view,'overlayNum',overlayName);
elseif isnumeric(overlayName)
	overlayNum = overlayName;
	overlayName = viewGet(view,'overlayName',overlayNum);
else
	myErrorDlg(['Bad overlay name: ',overlayName]);
end

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

% Path
filename = [analysisName,'-',overlayName];
pathStr = fullfile(viewGet(view,'overlayDir'),filename);
        
% Assign local variable with overlayName = overlay
overlay = viewGet(view,'overlay',overlayNum,analysisNum);
eval([overlayName,'=overlay;']);

% Write, though check for over-writing
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
    if confirm
        saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
            'Save Overlay?','Yes','No','No');
	end
end
if strcmp(saveFlag,'Yes')
    fprintf('Saving %s...',pathStr);
    saveString = ['save(pathStr,','''',overlayName,'''',');'];
    eval(saveString);
    fprintf('done\n');
else
    fprintf('Overlays not saved...');
end
return;
