% importOverlay.m
%
%        $Id$ 
%      usage: [thisView,params] = importOverlay(thisView,params,<options>)
%         by: julien besle
%       date: 20/01/2010
%    purpose: import overlay in an analysis
%
%             to just get a default parameter structure:
%             [thisView,params] = importOverlay(thisView,[],'justGetParams=1','defaultParams=1','pathname=fileToImport');

function [thisView,params] = importOverlay(thisView,params,varargin)

% check arguments
if nargin<1
  help importOverlay
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

if ~isanalysis(viewGet(thisView,'analysis'))
    mrWarnDlg('(importOverlay) Overlays must be imported into an analysis. Use Edit -> Analysis -> New Analysis.')
    return
end


nScans = viewGet(thisView,'nScans');
if ieNotDefined('params')
  %first check that there is at least one scan and display scan selection dialog if more than one
  switch(nScans)
    case 0
      mrWarnDlg('(importOverlay) Group must have at least one scan.')
      return
    case 1
      params.scanlist = 1;
    otherwise
      if defaultParams
        params.scanlist = 1:nScans;
      else
        params.scanlist = selectInList(thisView,'scans');
      end
  end


  if ieNotDefined('pathname')
    
    if defaultParams
      mrWarnDlg('(importOverlay) Specify a filename when using option defaultParams');
      return
    end
    [filename, pathname] = uigetfile({sprintf('*%s',mrGetPref('niftiFileExtension')),'Nifti files'},'Select nifti file that you want to import');
    if (filename==0)
      return
    end

    % get the full file name
    params.pathname = fullfile(pathname,filename);
  else
    params.pathname = pathname;
  end
end

% read the nifti header
[data,hdr] = cbiReadNifti(params.pathname);
if isempty(hdr),return,end

% make sure it has only 1 frame
if hdr.dim(1)==3
   hdr.dim(5)=1;
end
if hdr.dim(5) < 1
  mrWarnDlg(sprintf('(importOverlay) Could not import image because it has %d frames',hdr.dim(5)));
  return
elseif hdr.dim(5) >= 1
  if fieldIsNotDefined(params,'frameList')
    if defaultParams || hdr.dim(5)==1
      params.frameList=1:hdr.dim(5);
    else
      params.frameList = buttondlg('Choose the frames to import', [cellstr(int2str((1:hdr.dim(5))'));{'All'}]);
      if isempty(params.frameList)
        return
      end
      if params.frameList(end)
        params.frameList = 1:hdr.dim(5);
      else
        params.frameList = find(params.frameList);
      end
    end
  end
end
nFrames = length(params.frameList);


for i_scan = params.scanlist
scan_dims = viewGet(thisView,'datasize',i_scan);
   if any(hdr.dim([2 3 4])'~= scan_dims)
      mrWarnDlg(sprintf('(importOverlay) Could not import image because it is not compatible with scan (%s vs %s)', num2str(hdr.dim([2 3 4])'),num2str(scan_dims)) );
      return
   end
end

maxFrames=20;
if fieldIsNotDefined(params,'nameFrame') && fieldIsNotDefined(params,'nameFrame1')
  [~,filename,ext]=fileparts(params.pathname);
  filename = [filename ext];
  paramsInfo = {{'filename',filename,'The name of the nifti file that you are importing','editable=0'}};
  if nFrames<=maxFrames
    for iFrame = 1:nFrames
      numFrame = num2str(params.frameList(iFrame));
      paramsInfo{end+1} = {['nameFrame' numFrame],['Frame ' numFrame],['A description for the overlay corresponding to frame' numFrame]};
    end
  else
    paramsInfo{end+1} = {'nameFrame','','A description for the overlay you are importing. the frame number will be added at the end for each overlay'};
  end
  if defaultParams
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo);
  end
  if isempty(tempParams),return,end
  params = copyFields(tempParams,params);
  params = rmfield(params,{'paramInfo'});
end

if justGetParams,return,end

if nFrames>maxFrames
  for iFrame = 1:nFrames
    numFrame = num2str(params.frameList(iFrame));
    params.(['nameFrame' numFrame]) = [params.nameFrame '(Frame ' numFrame ')'];
  end
end
  
max_overlay = max(data(:));
min_overlay = min(data(:));
defaultOverlay.name = '';
defaultOverlay.groupName = viewGet(thisView,'groupName');
defaultOverlay.function = 'importOverlay';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.data = cell(nScans,1);
defaultOverlay.date = datestr(now);
defaultOverlay.params = params;
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = '';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [min_overlay max_overlay];
defaultOverlay.clip = [min_overlay max_overlay];

for iFrame=1:nFrames
  numFrame = params.frameList(iFrame);
  overlays(iFrame) = defaultOverlay;
  overlays(iFrame).name = [params.(['nameFrame' num2str(numFrame)]) ' (' params.filename ')']; 
  overlays(iFrame).data(params.scanlist) = squeeze(num2cell(repmat(data(:,:,:,numFrame),[1 1 1 length(params.scanlist)]),[1 2 3]))';
end

thisView = viewSet(thisView,'newoverlay',overlays);
refreshMLRDisplay(thisView.viewNum);





