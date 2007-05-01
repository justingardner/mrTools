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
    confirm = viewGet([],'pref','verbose');
    if isempty(confirm)
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
  
[pathStr filename] = fileparts(pathStr);
filename = sprintf('%s.mat',filename);
if isfile(fullfile(pathStr,filename))
  % get the preference how to deal with what to do with over-writing
  saveMethod = viewGet(view,'prefs','overwritePolicy');
  saveMethodTypes = {'Merge','Rename','Overwrite'};
  saveMethod = find(strcmp(saveMethod,saveMethodTypes));
  if confirm || isempty(saveMethod) || (saveMethod==0)
    % ask the user what to do
    paramsInfo = {{'saveMethod',saveMethodTypes,'Choose how you want to save the variable'}};

    params = mrParamsDialog(paramsInfo,sprintf('%s already exists',filename));
    if isempty(params),return,end
    saveMethod = find(strcmp(params.saveMethod,saveMethodTypes));
  end
  % now we know what to do, do the different options
  if saveMethod == 1
    disp(sprintf('(saveAnalysis) Merging with old analysis'));
    % load up old analysis
    oldAnal = load(fullfile(pathStr,filename));
    oldAnalFieldnames = fieldnames(oldAnal);
    oldAnal = oldAnal.(oldAnalFieldnames{1});
    % get the new analysis
    newAnal = eval(analysisName);
    % ok, now look at overlays field
    for oldNum = 1:length(oldAnal.overlays)
      oldOverlay = oldAnal.overlays(oldNum);
      matchedOverlay = 0;
      for newNum = 1:length(newAnal.overlays)
	% see if there is a matching overlay
	if strcmp(oldOverlay.name,newAnal.overlays(newNum).name)
	  matchedOverlay = 1;
	  % now try to combine them
	  for scanNum = 1:length(oldOverlay.data)
	    % if it is not empty in the new structure, but
	    % is in there in the old strucutre, copy it over
	    if ((scanNum >= length(newAnal.overlays(newNum).data)) || isempty(newAnal.overlays(newNum).data{scanNum})) && ~isempty(oldOverlay.data{scanNum})
	      newAnal.overlays(newNum).data{scanNum} = oldOverlay.data{scanNum};
	      if (isfield(newAnal.overlays(newNum),'params') && ...
		  iscell(newAnal.overlays(newNum).params) && ...
		  isempty(newAnal.overlays(newNum).params{newNum}))
		% copy over overlay params
		newAnal.overlays(newNum).params{newNum} = oldOverlay.params{scanNum};
	      end
	      disp(sprintf('(saveAnalysis) Merged overlay %s:%i from old analysis',oldOverlay.name,scanNum));
	    end
	  end
	end
      end
      % if there was no match, then simply add the overlay to the analysis struct
      if ~matchedOverlay
	newAnal.overlays(end+1) = oldAnal.overlays(oldNum);
	disp(sprintf('(saveAnalysis) Added overlay %s from old analysis',oldOverlay.name));
      end
    end
    % if there are any d structures, then do the same thing
    if isfield(oldAnal,'d') && isfield(newAnal,'d')
      for i = 1:length(oldAnal.d)
	% if the new analysis doesn't have it but the old one does
	if (((length(newAnal.d)<i) || isempty(newAnal.d{i})) &&...
	    (length(oldAnal.d)>=i) && ~isempty(oldAnal.d{i}))
	  newAnal.d{i} = oldAnal.d{i};
	  disp(sprintf('(saveAnalysis) Merged d{%i} from old analysis',i));
	end
      end
    end
    % if the params has fields called scanParams than those should
    % be merged
    if (isfield(oldAnal,'params') && isfield(oldAnal.params,'scanParams') &&...
	isfield(newAnal,'params') && isfield(newAnal.params,'scanParams'))
      for i = 1:length(oldAnal.params.scanParams)
	% if the new analysis doesn't have it but the old one does
	if (((length(newAnal.params.scanParams)<i) || ...
	    isempty(newAnal.params.scanParams{i})) && ...
	    (length(oldAnal.params.scanParams)>=i) && ...
	    ~isempty(oldAnal.params.scanParams{i}))
	  newAnal.params.scanParams{i} = oldAnal.params.scanParams{i};
	  disp(sprintf('(saveAnalysis) Merged params.scanParams{%i} from old analysis',i));
	end
      end
    end
    eval(sprintf('%s = newAnal;',analysisName));
    % now the save name gets the same and it can overwrite the old one
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

if strcmp(saveFlag,'Yes')
  pathStr = fullfile(pathStr,filename);
  fprintf('Saving %s...',pathStr);
  saveString = ['save(pathStr,','''',analysisName,'''',');'];
  eval(saveString);
  fprintf('done\n');
else
  fprintf('Analysis not saved...');
end
return;
