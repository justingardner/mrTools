% flexibleCopyPastePlugin.m
%
%        $Id$ 
%      usage: flexibleCopyPastePlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%    purpose: replaces behaviour of copy/paste overlay in Edit menu
%
function retval = flexibleCopyPastePlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help flexibleCopyPasteOverlaysPlugin
  return
end

switch action
 % return a help string
 case {'help','h','?'}
   retval = 'Allows copy/paste of several scans/overlays and automatic interpolation of overlays to the destination space';
   
 case {'install','i'}
  % check for a valid thisView
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(flexibleCopyPasteOverlaysPlugin) Need a valid thisView to install plugin'));
  else
    % if the thisView is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    mlrAdjustGUI(thisView,'set','/Edit/Overlay/Copy Overlay','Callback',@copyOverlayCallback);
    mlrAdjustGUI(thisView,'set','/Edit/Overlay/Paste Overlay','Callback',@pasteOverlayCallback);

    mlrAdjustGUI(thisView,'set','copyScanMenuItem','Callback',@copyScanCallback);
    mlrAdjustGUI(thisView,'set','pasteScanMenuItem','Callback',@pasteScanCallback);
    
    
    % return true to indicate successful plugin
    retval = true;
  end
   
 otherwise
   disp(sprintf('(flexibleCopyPasteOverlaysPlugin) Unknown command %s',action));
end

% --------------------------------------------------------------------
function copyOverlayCallback(hObject, dump)
mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
clipboard = copyOverlay(thisView); %calls copyOverlay which asks which overlay and for which scans to copy an returns all the copied overlays
if ~isempty(clipboard)
  MLR.clipboard = clipboard;
end

% --------------------------------------------------------------------
function pasteOverlayCallback(hObject, dump)
mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
pasteOverlay(thisView, MLR.clipboard);
refreshMLRDisplay(thisView.viewNum);


% --------------------------------------------------------------------
function copyScanCallback(hObject, dump)
mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
if viewGet(thisView,'nScans')
  scanList = selectInList(thisView,'scans','Select Scans to copy');
else
  scanList = 1;
end
cScan=0;
for iScan = scanList
  cScan = cScan+1;
  scan(cScan).scanParams = viewGet(thisView,'scanParams',iScan);
  scan(cScan).scanParams.fileName = [];
  % just save the filename instead of loading the whole tSeries
  % since some scans are very large
  %scan.tseries = loadTSeries(thisView,cScan,'all');
  scan(cScan).tseries = viewGet(thisView,'tseriesPathStr',iScan);
end
MLR.clipboard = scan;

% --------------------------------------------------------------------
function pasteScanCallback(hObject, dump)
mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
if isfield(MLR.clipboard,'tseries') && isfield(MLR.clipboard,'scanParams') 
  for iScan = 1:length(MLR.clipboard)
    if isscan(MLR.clipboard(iScan).scanParams)
    	saveNewTSeries(thisView,MLR.clipboard(iScan).tseries,MLR.clipboard(iScan).scanParams,MLR.clipboard(iScan).scanParams.niftiHdr);
    else
      mrWarnDlg(['(paste scan) Could not paste scan ' MLR.clipboard(iScan).tseries ' because its parameters are not valid'])
    end
  end
else
    mrErrorDlg('(paste scan) Cannot paste. Clipboard does not contain valid scans. Use Edit -> Scan -> Copy Scan.')
end
