% importOverlay.m
%
%        $Id$ 
%      usage: importOverlay(thisView)
%         by: julien besle
%       date: 20/01/2010
%    purpose: import overlay in an analysis of mrLoadRet GUI
%
function retval = importOverlay(thisView)

% check arguments
if ~any(nargin == [1])
  help importOverlay
  return
end

if ~isanalysis(viewGet(thisView,'analysis'))
    mrWarnDlg('(importOverlay) Overlays must be imported into an analysis. Use Edit -> Analysis -> New Analysis.')
    return
end

% go find the group that user wants to load here
[filename pathname] = uigetfile({sprintf('*%s',mrGetPref('niftiFileExtension')),'Nifti files'},'Select nifti file that you want to import');
if (filename==0)
  return
end

% get the full file name
fullFilename = fullfile(pathname,filename);

% read the nifti header
[data,hdr] = cbiReadNifti(fullFilename);
if isempty(hdr),return,end

% make sure it has only 1 frame
if hdr.dim(1)==3
   hdr.dim(5)=1;
end
if hdr.dim(5) ~= 1
  mrWarnDlg(sprintf('(importOverlay) Could not import image because it has %d frames',hdr.dim(5)));
  return
end

scan_number = viewGet(thisView,'nScans');
if scan_number>1
   scanlist = selectInList(thisView,'scans');
else
   scanlist = 1;
end

for i_scan = scanlist
scan_dims = viewGet(thisView,'datasize',i_scan);
   if any(hdr.dim([2 3 4])'~= scan_dims)
      mrWarnDlg(sprintf('(importOverlay) Could not import image because it is not compatible with scan (%s vs %s)', num2str(hdr.dim([2 3 4])'),num2str(scan_dims)) );
      return
   end
end
   
paramsInfo = {{'filename',filename,'The name of the nifti file that you are importing','editable=0'}};
paramsInfo{end+1} = {'description','','A description for the overlay you are importing'};

params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

max_overlay = max(data(:));
min_overlay = min(data(:));
overlay.name = [params.description ' (' params.filename ')']; 
overlay.groupName = viewGet(thisView,'groupName');
overlay.function = 'importOverlay';
overlay.reconcileFunction = 'defaultReconcileParams';
overlay.data = cell(scan_number,1);
overlay.data(scanlist) = squeeze(num2cell(repmat(data,[1 1 1 length(scanlist)]),[1 2 3]))';
overlay.date = datestr(now);
overlay.params = params;
% colormap is made with a little bit less on the dark end
overlay.colormap = hot(312);
overlay.colormap = overlay.colormap(end-255:end,:);
overlay.alpha = 1;
overlay.interrogator = '';
overlay.mergeFunction = 'defaultMergeParams';
overlay.colormapType = 'normal';
overlay.range = [min_overlay max_overlay];
overlay.clip = [min_overlay max_overlay];

thisView = viewSet(thisView,'newoverlay',overlay);
refreshMLRDisplay(thisView.viewNum);





