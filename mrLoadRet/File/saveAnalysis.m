function saveAnalysis(view,analysisName,confirm)
%
%   saveAnalysis(view,[analysisName],[confirm])
%
% Saves an analysis to a file. The filename = analysis.name.
%
% analysisname: Can be either the name or the number.
%          Default: current analysis
% confirm: If nonzero, prompts user to determine what to do if filename
%          already exists. Otherwise, relies of 'overwritePolicy'
%          preference to choose what to do. Default: 0.
%          
%
% djh, 7/2006
% jlg, 4/2007 added stuff to handle overwrite

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
analysisdir = viewGet(view,'analysisdir',[],analysisNum);
filename = [analysisName,'.mat'];
pathStr = analysisdir;

% Assign local variable with analysisName = analysis
analysis = viewGet(view,'analysis',analysisNum);
eval([analysisName,'=analysis;']);

% Write, though check for over-writing
if isfile(fullfile(pathStr,filename))
    % get the preference how to deal with what to do with over-writing
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
        disp(sprintf('(saveAnalysis) Merging with old analysis'));
        % load the old analysis
        s = load(fullfile(pathStr,filename));
        varNames = fieldnames(s);
        oldAnal = s.(varNames{1});
        oldAnal.name = varNames{1};
        % get the new analysis
        newAnal = eval(analysisName);    
        % check if they have the same name and merge them
        if strcmp(oldAnal.name,newAnal.name)
            [mergedParams,mergedData] = feval(newAnal.mergeFunction,newAnal.groupName,...
                oldAnal.params,newAnal.params);            
            newAnal.params = mergedParams;
            % loop through overlays and merge the params and data
            for oldNum = 1:length(oldAnal.overlays)
                oldOverlay = oldAnal.overlays(oldNum);
                matchedOverlay = 0;
                for newNum = 1:length(newAnal.overlays)
                    newOverlay = newAnal.overlays(newNum);
                    % see if there is a matching overlay
                    if strcmp(oldOverlay.name,newOverlay.name)
                        matchedOverlay = 1;
                        % now try to combine them
                        [mergedParams,mergedData] = feval(newOverlay.mergeFunction,newOverlay.groupName,...
                            oldOverlay.params,newOverlay.params,oldOverlay.data,newOverlay.data);
                        newOverlay.params = mergedParams;
                        newOverlay.data = mergedData;
                        newAnal.overlays(newNum) = newOverlay;
                        disp(sprintf('(saveAnalysis) Merged overlay %s from old analysis',oldOverlay.name));
                    end
                end
                % if there was no match, then simply add the overlay to the analysis struct
                if ~matchedOverlay
                    newAnal.overlays(end+1) = oldOverlay;
                    disp(sprintf('(saveAnalysis) Added overlay %s from old analysis',oldOverlay.name));
                end
            end                
            % set analysisName variable so that it will be saved below
            eval(sprintf('%s = newAnal;',analysisName));
            % replace analysis with the newly merged one
            view = viewSet(view,'deleteAnalysis',analysisNum);
            view = viewSet(view,'newAnalysis',newAnal);            
        else
            mrWarnDlg('(saveAnalysis) Merge failed. Save analysis aborted.');
            return
        end
 
    elseif saveMethod == 2
        % put up a dialog to get new save name
        [filename pathStr] = uiputfile({'*.mat'},'Enter new name to save analysis as',fullfile(pathStr,filename));
        % user hit cancel
        if isequal(filename,0)
            return
        end
        % otherwise accept the filename, make sure it has a .mat extension
        filename = sprintf('%s.mat',stripext(filename));
        
    elseif saveMethod == 3
        % this is the easiest, just overwrite
        disp(sprintf('(saveAnalysis) Overwriting old analysis'));
    end
end

% Finally, write the file
pathStr = fullfile(pathStr,filename);
fprintf('Saving %s...',pathStr);
saveString = ['save(pathStr,','''',analysisName,'''',');'];
eval(saveString);
fprintf('done\n');

return;
