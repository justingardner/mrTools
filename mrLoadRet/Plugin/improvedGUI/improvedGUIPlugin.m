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
    
    % new uicontrols and reposition old ones
    if viewGet(thisView,'baseType')>0
      corticalDepthVisibility = 'on';
    else
      corticalDepthVisibility = 'off';
    end    
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
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','SliderStep',min(1/viewGet(thisView,'corticalDepthBins')*[1 3],1));
    mlrAdjustGUI(thisView,'set','corticalDepthText','position',   [0.24    0.815   0.04   0.025]);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthSlider',...
        'SliderStep',min(1/(viewGet(thisView,'corticalDepthBins')-1)*[1 3],1),...
        'Callback',@corticalMaxDepthSlider_Callback,'visible',corticalDepthVisibility,...
                          'value',.5,'style','slider','position', [0.09    0.79    0.14   0.02 ]);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthText',...
        'Callback',@corticalMaxDepthText_Callback,'visible',corticalDepthVisibility,'String','0.5',...
              'BackgroundColor','white','style','edit','position',[0.24    0.78    0.04   0.03 ]);
    mlrAdjustGUI(thisView,'add','control','linkMinMaxDepthCheck','style','checkbox','value',1,...
        'fontSize', 12,'visible',corticalDepthVisibility,...
        'String','Fix Depth Range','position',                    [0.04    0.755   0.23   0.025]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','position', [0.01    0.78    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','axialRadioButton','position',    [0.11    0.78    0.07   0.025]);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','position',  [0.19    0.78    0.1    0.025]);
    mlrAdjustGUI(thisView,'set','baseImage','position',           [0.01    0.71    0.06   0.025]);
    mlrAdjustGUI(thisView,'set','basePopup','position',           [0.07    0.705   0.21   0.035]);
    mlrAdjustGUI(thisView,'set','baseGamma','string','Gamma');
    mlrAdjustGUI(thisView,'set','baseGamma','position',           [0.02    0.665   0.07   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseGammaSlider','position',     [0.10    0.67    0.13   0.02 ]);
    mlrAdjustGUI(thisView,'set','baseGammaText','position',       [0.24    0.67    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','rotate','position',              [0.02    0.63    0.06   0.03 ]);
    mlrAdjustGUI(thisView,'set','rotateSlider','position',        [0.09    0.635   0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','rotateText','position',          [0.24    0.635   0.04   0.025]);
    mlrAdjustGUI(thisView,'set','baseTilt','position',            [0.02    0.595   0.03   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseTiltSlider','position',      [0.055   0.60    0.175  0.02 ]);
    mlrAdjustGUI(thisView,'set','baseTiltText','position',        [0.24    0.60    0.04   0.025]);
    mlrAdjustGUI(thisView,'set','analysis','position',            [0.01    0.545   0.08   0.025]);
    mlrAdjustGUI(thisView,'set','analysisPopup','position',       [0.09    0.54    0.19   0.035]);
    mlrAdjustGUI(thisView,'set','overlay','position',             [0.01    0.51    0.07   0.025]);
    mlrAdjustGUI(thisView,'add','control','clipAcrossOverlays','style','checkbox','value',0,...
        'callback',@clipAcrossOverlays_Callback,...
        'fontSize', 11,'String','Clip across overlays','position',[0.09    0.505   0.18   0.025]);
    mlrAdjustGUI(thisView,'add','control','colorBlending','style','text',...
        'fontSize', 12,'String','Blending mode','position',       [0.02    0.475   0.12   0.025]);
    mlrAdjustGUI(thisView,'add','control','colorBlendingPopup','style','popupmenu','value',2,...
        'callback',@colorBlendingPopup_Callback,'String',{'Alpha blend' 'Additive'},...
        'fontSize', 12,'position',                                [0.14    0.47    0.14   0.035]);
    mlrAdjustGUI(thisView,'set','overlayPopup','position',        [0.02    0.24    0.26   0.23 ]);
    mlrAdjustGUI(thisView,'set','overlayPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','overlayPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','overlayMin','position',          [0.02    0.205   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','overlayMin','string','Min');
    mlrAdjustGUI(thisView,'set','overlayMinText','position',      [0.055   0.21    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','mapMax','position',              [0.15    0.205   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','mapMax','string','Max');
    mlrAdjustGUI(thisView,'set','overlayMaxText','position',      [0.19    0.21    0.09   0.025]);
    mlrAdjustGUI(thisView,'set','overlayMinSlider','position',    [0.02    0.185   0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','overlayMaxSlider','position',    [0.02    0.165   0.26   0.02 ]);
    mlrAdjustGUI(thisView,'set','alpha','position',               [0.02    0.125   0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','alphaSlider','position',         [0.08    0.13    0.14   0.02 ]);
    mlrAdjustGUI(thisView,'set','alphaText','position',           [0.23    0.125   0.05   0.025]);
    mlrAdjustGUI(thisView,'set','colorbar','position',            [0.35    0.12    0.58   0.07 ]);
    mlrAdjustGUI(thisView,'set','colorbar','fontSize',10);
    mlrAdjustGUI(thisView,'add','axes','colorbarRightBorder',...
      'YaxisLocation','right','XTick',[],'box','off','position',  [0.929   0.12    0.001  0.07 ]);

    % Add some useful menus
    mlrAdjustGUI(thisView,'add','menu','Unlink Stimfile','/Edit/Scan/Link Stimfile','callback',@unlinkStimfileMenuItem_Callback);
    mlrAdjustGUI(thisView,'set','/Edit/Scan/Link Stimfile','separator','on');
    
    
    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Improvements include: superimposition of several overlays, average over a range of cortical depths, ...';
 otherwise
   disp(sprintf('(improvedGUIPlugin) Unknown command %s',action));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function corticalMaxDepthSlider_Callback(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = get(hObject,'Value');
if get(handles.linkMinMaxDepthCheck,'value')
  minDepth = viewGet(viewNum,'corticalMinDepth');
  maxDepth = str2num(get(handles.corticalMaxDepthText,'string'));
  newMinValue = minDepth+newMaxValue-maxDepth;
  if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
    mlrGuiSet(viewNum,'corticalMinDepth',newMinValue);
    mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
    drawnow;
    refreshMLRDisplay(viewNum);
  else
    set(hObject,'Value',maxDepth);
  end
else
  mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
  refreshMLRDisplay(viewNum);
end


function corticalMaxDepthText_Callback(hObject,dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = str2num(get(hObject,'String'));
if isempty(newMaxValue) %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.corticalMaxDepthSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  if get(handles.linkMinMaxDepthCheck,'value')
    minDepth = viewGet(viewNum,'corticalMinDepth');
    maxDepth = get(handles.corticalMaxDepthSlider,'value');
    newMinValue = minDepth+newMaxValue-maxDepth;
    if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
      mlrGuiSet(viewNum,'corticalMinDepth',newMinValue);
      mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
      refreshMLRDisplay(viewNum);
    else
      set(hObject,'string',num2str(maxDepth));
    end
  else
    mlrGuiSet(viewNum,'corticalMaxDepth',newMaxValue);
    refreshMLRDisplay(viewNum);
  end
end

function clipAcrossOverlays_Callback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
refreshMLRDisplay(viewNum);

function colorBlendingPopup_Callback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function unlinkStimfileMenuItem_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum, 'stimfilename', '');
    