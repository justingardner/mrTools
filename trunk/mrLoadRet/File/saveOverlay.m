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
% confirm: If nonzero, prompts user to determine what to do if filename
%          already exists. Otherwise, relies of 'overwritePolicy'
%          preference to choose what to do. Default: 0.
%
% djh, 7/2004
% djh, 7/2006 updated with analysisName
% djh, 5/2007 added stuff to handle overwrite (copied from saveAnalysis)


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
    confirm = 0;
end

% Path
filename = [analysisName,'-',overlayName,'.mat'];
pathStr = viewGet(view,'overlayDir');

% Assign local variable with overlayName = overlay
overlay = viewGet(view,'overlay',overlayNum,analysisNum);
eval([overlayName,'=overlay;']);

% Check for over-writing
if isfile(fullfile(pathStr,filename))

    % get overwrite preference or ask user
    saveMethod = mrGetPref('overwritePolicy');
    saveMethodTypes = {'Merge','Rename','Overwrite'};
    saveMethod = find(strcmp(saveMethod,saveMethodTypes));
    if confirm || isempty(saveMethod) || (saveMethod==0)
        % ask the user what to do
        paramsInfo = {{'saveMethod',saveMethodTypes,'Choose how you want to save the variable'}};
        params = mrParamsDialog(paramsInfo,sprintf('%s already exists',filename));
        if isempty(params),return,end
        saveMethod = find(strcmp(params.saveMethod,saveMethodTypes));
    end

    if saveMethod == 1
        disp('(saveOverlay) Merging with old analysis');
        % load the old overlay
        s = load(fullfile(pathStr,filename));
        varNames = fieldnames(s);
        oldOverlay = s.(varNames{1});
        oldOverlay.name = varNames{1};
        % get the new overlay
        newOverlay = eval(overlayName);
        % check if they have the same name and merge them
        if strcmp(oldOverlay.name,newOverlay.name)
            [mergedParams,mergedData] = feval(newOverlay.mergeFunction,newOverlay.groupName,...
                oldOverlay.params,newOverlay.params,oldOverlay.data,newOverlay.data);
            newOverlay.params = mergedParams;
            newOverlay.data = mergedData;
            % set overlayName variable so that it will be saved below
            eval(sprintf('%s = newOverlay;',overlayName));
            % replace overlay with the newly merged one
            view = viewSet(view,'deleteoverlay',overlayNum,analysisNum);
            view = viewSet(view,'newoverlay',newOverlay,analysisNum);
        else
            mrWarnDlg('(saveOverlay) Merge failed. Save overlay aborted.');
            return
        end

    elseif saveMethod == 2
        % put up a dialog to get new filename
        [filename pathStr] = uiputfile({'*.mat'},'Enter new name to save overlay',fullfile(pathStr,filename));
        % user hit cancel
        if isequal(filename,0)
            return
        end
        % otherwise accept the filename, make sure it has a .mat extension
        filename = [stripext(filename),'.mat'];

    elseif saveMethod == 3
        % this is the easiest, just overwrite
        disp('(saveAnalysis) Overwriting old analysis');
    end
end

% Finally, write the file
pathStr = fullfile(pathStr,filename);
fprintf('Saving %s...',pathStr);
saveString = ['save(pathStr,','''',overlayName,'''',');'];
eval(saveString);
fprintf('done\n');

return
