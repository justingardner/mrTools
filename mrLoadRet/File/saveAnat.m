function saveAnat(view,anatomyName,confirm)
%
% saveOverlay(view,[anatomyName],[confirm])
%
% Saves a base anatomy as a nifti file: anatomyName.img
%
% anatomyName: Can be either the name or the number (in which case it is
%          converted to the name). Default: current base volume.
% confirm: If filename already exists, prompt user to over-write. 
%          Default: uses mrLoadRet 'verbose' preference or 0 (if preference
%          not defined.
%
% djh, 7/2006

if ieNotDefined('anatomyName')
    anatomyNum = viewGet(view,'currentBase');
    anatomyName = viewGet(view,'baseName',anatomyNum);
end
if isstr(anatomyName)
    anatomyNum = viewGet(view,'baseNum',anatomyName);
elseif isnumeric(anatomyName)
	anatomyNum = anatomyName;
	anatomyName = viewGet(view,'baseName',anatomyNum);
else
	myErrorDlg(['Invalide base volume (anatomy) name: ',anatomyName]);
end

if ieNotDefined('confirm')
    confirm = viewGet([],'pref','verbose');
    if isempty(confirm)
        confirm = 0;
    end
end

% Extract data and hdr
baseVolume = viewGet(view,'baseVolume',anatomyNum);
data = baseVolume.data;
hdr = baseVolume.hdr;

% Path
pathStr = fullfile(viewGet(view,'anatomydir'),[anatomyName,niftiFileExtension]);
        
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
    [byteswritten,hdr] = cbiWriteNifti(pathStr,data,hdr);
    fprintf('done\n');
else
    fprintf('Anatomy not saved...');
end
return;
