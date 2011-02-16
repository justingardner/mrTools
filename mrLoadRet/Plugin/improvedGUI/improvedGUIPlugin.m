% improvedGUIPlugin.m
%
%        $Id: improvedGUIPlugin.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: improvedGUIPlugin(action,<thisView>)
%         by: julien besle
%       date: 13/02/11
%    purpose: 
%

function retval = improvedGUIPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help improvedGUIPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(improvedGUIPlugin) Need a valid view to install plugin'));
  else
    
    
    mlrAdjustGUI(thisView,'set','group','position',               [0.01    0.955   0.06   0.025]);
    mlrAdjustGUI(thisView,'set','groupPopup','position',          [0.08    0.95    0.2    0.035]);
    mlrAdjustGUI(thisView,'set','ROI','position',                 [0.01    0.91    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','roiPopup','position',            [0.06    0.905   0.22   0.035]);
    mlrAdjustGUI(thisView,'set','scan','position',                [0.01    0.85    0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','scanSlider','position',          [0.07    0.855   0.16   0.02 ]);
    mlrAdjustGUI(thisView,'set','scanText','position',            [0.24    0.855   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','slice','position',               [0.01    0.81    0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','sliceSlider','position',         [0.07    0.815   0.16   0.02 ]);
    mlrAdjustGUI(thisView,'set','sliceText','position',           [0.24    0.815   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','corticalDepth','position',       [0.01    0.78    0.07   0.06 ]);
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','position', [0.09    0.815   0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','corticalDepthText','position',   [0.24    0.815   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','corticalMaxDepthText','units', 'normalized');
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthSlider',...
        'SliderStep',min(1/mrGetPref('corticalDepthBins')*[1 3],1),...
        'Callback',@corticalMaxDepthSlider_Callback,...
                                     'style','slider','position', [0.09    0.79    0.14   0.02 ]);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthText',...
        'Callback',@corticalMaxDepthText_Callback,...
              'BackgroundColor','white','style','edit','position',[0.24    0.785   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','position', [0.01    0.78    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','axialRadioButton','position',    [0.11    0.78    0.07   0.025]);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','position',  [0.19    0.78    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','baseImage','position',           [0.01    0.725   0.06   0.025]);
    mlrAdjustGUI(thisView,'set','basePopup','position',           [0.07    0.72    0.21   0.035]);
    mlrAdjustGUI(thisView,'set','baseGamma','string','Gamma');
    mlrAdjustGUI(thisView,'set','baseGamma','position',           [0.02    0.68    0.07   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseGammaSlider','position',     [0.10    0.685   0.13   0.02 ]);
    mlrAdjustGUI(thisView,'set','baseGammaText','position',       [0.24    0.685   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','rotate','position',              [0.02    0.64    0.06   0.03 ]);
    mlrAdjustGUI(thisView,'set','rotateSlider','position',        [0.09    0.645   0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','rotateText','position',          [0.24    0.645   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','baseTilt','position',            [0.02    0.60    0.03   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseTiltSlider','position',      [0.055   0.605   0.175  0.02 ]);
    mlrAdjustGUI(thisView,'set','baseTiltText','position',        [0.24    0.605   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','analysis','position',            [0.01    0.55    0.08   0.025]);
    mlrAdjustGUI(thisView,'set','analysisPopup','position',       [0.09    0.545   0.19   0.035]);
    mlrAdjustGUI(thisView,'set','overlay','position',             [0.01    0.51    0.07   0.025]);
    mlrAdjustGUI(thisView,'set','overlayPopup','position',        [0.02    0.28    0.26   0.23 ]);
    mlrAdjustGUI(thisView,'set','overlayPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','overlayPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','overlayMin','position',          [0.02    0.235   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','overlayMin','string','Min');
    mlrAdjustGUI(thisView,'set','overlayMinText','position',      [0.055   0.24    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','mapMax','position',              [0.15    0.235   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','mapMax','string','Max');
    mlrAdjustGUI(thisView,'set','overlayMaxText','position',      [0.19    0.24    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','overlayMinSlider','position',    [0.02    0.21    0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','overlayMaxSlider','position',    [0.02    0.19    0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','alpha','position',               [0.02    0.15    0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','alphaSlider','position',         [0.08    0.155   0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','alphaText','position',           [0.23    0.155   0.05   0.025]);
    mlrAdjustGUI(thisView,'set','colorbar','position',            [0.35    0.12    0.58   0.07]);
    mlrAdjustGUI(thisView,'add','axes','colorbarRightBorder',...
      'YaxisLocation','right','XTick',[],'box','off','position',  [0.929   0.12    0.001  0.07]);

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Improvements include: superimposition of several overlays, ...';
 otherwise
   disp(sprintf('(improvedGUIPlugin) Unknown command %s',action));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function corticalMaxDepthSlider_Callback(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
value = get(hObject,'Value');
mlrGuiSet(viewNum,'corticalMaxDepth',value);
refreshMLRDisplay(viewNum);


function corticalMaxDepthText_Callback(hObject,dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if isempty(value) %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.corticalMaxDepthSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  mlrGuiSet(viewNum,'corticalMaxDepth',value);
  refreshMLRDisplay(viewNum);
end


